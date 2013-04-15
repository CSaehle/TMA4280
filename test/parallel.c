#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <mpi.h>


MPI_Datatype send_t, recv_t;

#define ROWS 3
#define COLS 3
#define NODES 2
#define RPN 2
#define PRINT_NODE 1
// Rows per node

void prn(int *a, int size);
int **create2DArray(int rows, int cols);

void local_transpose (int **b, int rpn, int off) {
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


  int **send = create2DArray(RPN, RPN*RPN);
  int **recv = create2DArray(RPN, RPN*RPN);
  int i, j;

  for (i = 0; i < RPN; i++) {
    for (j = 0; j < COLS; j++) {
      int size = i*COLS + j  + RPN*COLS* mpi_rank;
      if (size < ROWS*COLS) {
        send[i][j] = size + 1;
      }
    }
  }

  if (mpi_rank == PRINT_NODE) {
    for (i = 0; i < RPN; i++) {
      prn(send[i], RPN*RPN);
    }
    printf("\n");
  }

  for (i = 0; i < RPN; i++) {
    MPI_Alltoall(send[i], RPN, MPI_INT, recv[i], RPN, MPI_INT, MPI_COMM_WORLD);
  }
  if (mpi_rank == PRINT_NODE) {
    for (i = 0; i < RPN; i++) {
      prn(recv[i], RPN*RPN);
    }
    printf("\n");
  }
  for (i = 0; i < RPN; i++) {
    local_transpose(recv, RPN, i*RPN);
  }

  if (mpi_rank == PRINT_NODE) {
    for (i = 0; i < RPN; i++) {
      prn(recv[i], RPN*RPN);
    }
    printf("\n");
  }
  for (i = 0; i < RPN; i++) {
    MPI_Alltoall(recv[i], RPN, MPI_INT, send[i], RPN, MPI_INT, MPI_COMM_WORLD);
  }
  for (i = 0; i < RPN; i++) {
    local_transpose(send, RPN, i*RPN);
  }

  if (mpi_rank == PRINT_NODE) {
    for (i = 0; i < RPN; i++) {
      prn(send[i], RPN*RPN);
    }
    printf("\n");
  }
 end:
  MPI_Finalize();
}

int **create2DArray(int rows, int cols) {
  int i, total;
  total = rows * cols;
  int **a;
  a =    (int **) malloc(rows * sizeof(int *));
  a[0] = (int *) malloc(total *sizeof(int));
  for (i=1; i < rows; i++) {
    a[i] = a[i-1] + cols;
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
