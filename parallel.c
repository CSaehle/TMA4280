/*
  C-program to solve the two-dimensional Poisson equation on
  a unit square using one-dimensional eigenvalue decompositions
  and fast sine transforms

  einar m. ronquist
  ntnu, october 2000
  revised, october 2001
*/

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <mpi.h>

typedef double Real;

/* function prototypes */
Real *createRealArray (int n);
Real **createReal2DArray (int m, int n);
void transpose (Real **bt, Real **b, int m);
void fst_(Real *v, int *n, Real *w, int *nn);
void fstinv_(Real *v, int *n, Real *w, int *nn);


main(int argc, char **argv)
{
  Real *diag, **b, **bt, *z;
  Real pi, h, local_max, global_max;
  int i, j, n, m, nn;
  int mpi_size, mpi_rank, mpi_work;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  /* the total number of grid points in each spatial direction is (n+1) */
  /* the total number of degrees-of-freedom in each spatial direction is (n-1) */
  /* this version requires n to be a power of 2 */

  if (argc < 2) {
    if (mpi_rank == 0){
      printf("need a problem size\n");
    }
    goto end;
  }

  n  = atoi(argv[1]);
  m  = n-1;
  // mpi_work is the amount of work needed to be done by each mpi node. The last
  // mpi node may do slightly less work than the others, but that's the closest
  // we'll get to proper load balancing.
  mpi_work = 1 + ((m - 1) / mpi_size);
  nn = 4*n;

  diag = createRealArray (m);
  b    = createReal2DArray (m,m);
  bt   = createReal2DArray (m,m);
  z    = createRealArray (nn);

  h    = 1./(Real)n;
  pi   = 4.*atan(1.);

  #pragma omp parallel for
  for (i=0; i < m; i++) { // Everyone generate this one
    diag[i] = 2.*(1.-cos((i+1)*pi/(Real)n));
  }

  #pragma omp parallel for
  for (j=0; j < m; j++) { // MPI
    for (i=0; i < m; i++) { // OMP
      b[j][i] = h*h; // Or should this be calculated on node 0 and distributed?
    }
  }

  #pragma omp parallel for
  for (j=0; j < m; j++) { // MPI cut + OMP
    fst_(b[j], &n, z, &nn);
  }

  transpose (bt,b,m);

  #pragma omp parallel for
  for (i=0; i < m; i++) { // MPI cut + OMP
    fstinv_(bt[i], &n, z, &nn);
  }

  #pragma omp parallel for
  for (j=0; j < m; j++) { // MPI
    for (i=0; i < m; i++) {
      bt[j][i] = bt[j][i]/(diag[i]+diag[j]);
    }
  }

  #pragma omp parallel for
  for (i=0; i < m; i++) { // MPI cut + OMP
    fst_(bt[i], &n, z, &nn);
  }

  transpose (b,bt,m);

  #pragma omp parallel for
  for (j=0; j < m; j++) { // MPI cut + OMP
    fstinv_(b[j], &n, z, &nn);
  }

  local_max = 0.0;
  // MPI, work in range (and handle last node overflow)
  #pragma omp parallel for
  // OMP - done here to avoid thread startup/shutdowns and cache misses
  for (j=mpi_rank * mpi_work; j < (mpi_rank + 1) * mpi_work && j < m; j++) {
    for (i=0; i < m; i++) {
      if (b[j][i] > local_max) local_max = b[j][i];
    }
  }

  MPI_Reduce(&local_max, &global_max, 1,
             MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  if (mpi_rank == 0) {
    printf (" umax = %e \n", global_max);
  }
 end:
  MPI_Finalize();
}

void transpose (Real **bt, Real **b, int m)
{
  int i, j;
  for (j=0; j < m; j++) {
    for (i=0; i < m; i++) {
      bt[j][i] = b[i][j];
    }
  }
}

Real *createRealArray (int n)
{
  Real *a;
  int i;
  a = (Real *)malloc(n*sizeof(Real));
  for (i=0; i < n; i++) {
    a[i] = 0.0;
  }
  return (a);
}

Real **createReal2DArray (int n1, int n2)
{
  int i, n;
  Real **a;
  a    = (Real **)malloc(n1   *sizeof(Real *));
  a[0] = (Real  *)malloc(n1*n2*sizeof(Real));
  for (i=1; i < n1; i++) {
    a[i] = a[i-1] + n2;
  }
  n = n1*n2;
  memset(a[0],0,n*sizeof(Real));
  return (a);
}
