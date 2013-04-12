#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <mpi.h>


MPI_Datatype send_t, recv_t;

#define ROWS 4
#define COLS 4
#define NODES 2
#define RPN 2
// Rows per node

int **create2DArray(int m, int n) {
  int i, total;
  total = m * n;
  int **a;
  a =    (int **) malloc(m * sizeof(int *));
  a[0] = (int *) malloc(total *sizeof(int));
  for (i=1; i < n; i++) {
    a[i] = a[i-1] + n;
  }
  memset(a[0], 0, total*sizeof(int));
  return a;
}

void prn(int *a, int size) {
  printf("[");
  int i;
  for (i = 0; i < size; i++) {
    printf("%d, ", a[i]);
  }
  printf("]\n");
}

void transpose (int **b, int rpn, int off) {
  int i, j, temp;
  for (j=0; j < rpn; j++) {
    for (i=j+1; i < rpn; i++) {
      temp = b[j][i+off];
      b[j][i+off] = b[i][j+off];
      b[i][j+off] = temp;
    }
  }
}

main(int argc, char **argv )
{
  int mpi_size, mpi_rank, mpi_work;
  
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  if (mpi_size != NODES) {
    if (mpi_rank == 0) {
      printf("Need exactly %d nodes\n", NODES);
    }
    goto end;
  }
  
  
  int **send = create2DArray(ROWS, COLS);
  int **recv = create2DArray(ROWS, COLS);
  int i;

  for (i = 0; i < RPN*COLS; i++) {
    send[0][i] = i + RPN*COLS * mpi_rank;
  }

  if (mpi_rank == 0) {
    for (i = 0; i < RPN; i++) {
      prn(send[i], COLS);
    }
  }
  
  for (i = 0; i < RPN; i++) {
    MPI_Alltoall(send[i], RPN, MPI_INT, recv[i], RPN, MPI_INT, MPI_COMM_WORLD);
  }
  for (i = 0; i < RPN; i++) {
    transpose(recv, RPN, i*RPN);    
  }

  
  if (mpi_rank == 0) {
    for (i = 0; i < RPN; i++) {
      prn(recv[i], COLS);
    }
  }
 end:
  MPI_Finalize();
}
