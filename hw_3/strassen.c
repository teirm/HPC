/**
 * University of Pittsburgh
 * Department of Computer Science
 * CS1645: Introduction to HPC Systems
 * Instructor Bryan Mills, PhD
 * Student: Cyrus Ramavarapu 
 * Implement Pthreads version of Strassen algorithm for matrix multiplication.
 */

#include "strassen.h"
#include "timer.h"
#include "io.h"

#include <pthread.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// Make these globals so threads can operate on them. You will need to
// add additional matrixes for all the M and C values in the Strassen
// algorithms.
int **A;
int **B;
int **C;
// Reference matrix, call simpleMM to populate.
int **R;

/* A Matrices needed for Strassen Algorithm */
int **A_11;
int **A_12;
int **A_21;
int **A_22;

/* B Matrices needed for Strassen Algorithm */
int **B_11;
int **B_12;
int **B_21;
int **B_22;

/* B Matrices needed for Strassen Algorithm */
int **C_11;
int **C_12;
int **C_21;
int **C_22;

/* M Matrices Needed for Strassen Algorithm */ 
int **M_1;
int **M_2;
int **M_3;
int **M_4;
int **M_5;
int **M_6;
int **M_7;

/* Checks if the matrix dimensions are a power of 2 */
int is_power_two(int val)
{
    return (val != 0) && ((val & (val-1)) == 0);
}

int compute_next_power_two(int val)
{
    int next_power_two;

    next_power_two = 1;

    while (next_power_two < val) {
        next_power_two = next_power_two << 1;
    }

    return next_power_two;
}

void padded_split(int **M, int N, int padded_size, int m_pos)
{

    int i;
    int j;

    if (m_pos == 0) {
        A_11 = allocMatrix(padded_size/2);
        A_12 = allocMatrix(padded_size/2);
        A_21 = allocMatrix(padded_size/2);
        A_22 = allocMatrix(padded_size/2); 
        
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                if (i < padded_size/2 && j < padded_size/2) {
//                    printf("%d %d\n", i, j);
                    A_11[i][j] = M[i][j];
                } else if (i >= padded_size/2 && j < padded_size/2) {
//                    printf("%d %d\n", i - padded_size/2, j);
                    A_21[i - padded_size/2][j] = M[i][j];
                } else if (i < padded_size/2 && j >= padded_size/2) {
//                    printf("%d %d\n", i, j - padded_size/2);
                    A_12[i][j - padded_size/2] = M[i][j];
                } else { 
//                    printf("%d %d\n", i - padded_size/2, j - padded_size/2);
                    A_22[i - padded_size/2][j - padded_size/2] = M[i][j];
                }
            }
        }
    } else {
        B_11 = allocMatrix(padded_size/2);
        B_12 = allocMatrix(padded_size/2);
        B_21 = allocMatrix(padded_size/2);
        B_22 = allocMatrix(padded_size/2); 
        
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                if (i < padded_size/2 && j < padded_size/2) {
//                    printf("%d %d\n", i, j);
                    B_11[i][j] = M[i][j];
                } else if (i >= padded_size/2 && j < padded_size/2) {
//                    printf("%d %d\n", i - padded_size/2, j);
                    B_21[i - padded_size/2][j] = M[i][j];
                } else if (i < padded_size/2 && j >= padded_size/2) {
//                    printf("%d %d\n", i, j - padded_size/2);
                    B_12[i][j - padded_size/2] = M[i][j];
                } else { 
//                    printf("%d %d\n", i - padded_size/2, j - padded_size/2);
                    B_22[i - padded_size/2][j - padded_size/2] = M[i][j];
                }
            }
        }
    }
}

// Stupid simple Matrix Multiplication, meant as example.
void simpleMM(int N) {
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
      for (int k=0; k<N; k++) {
	R[i][j] += A[i][k] * B[k][j];
      }
    }
  }
}

void matrix_add(int ***T_1, int **T_2, int dim)
{
    int i;
    int j;

    for(i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            (*T_1)[i][j] = (*T_1)[i][j] + T_2[i][j];
        }
    }
}

void matrix_sub(int ***T_1, int **T_2, int dim)
{
    int i;
    int j;

    for(i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            (*T_1)[i][j] = (*T_1)[i][j] - T_2[i][j];
        }
    }
}

// WRITE YOUR CODE HERE, you will need to also add functions for each
// of the sub-matrixes you will need to calculate but you can create your
// threads in this fucntion.
void strassenMM(int N) {
    
    int new_size; 
    new_size = N; 
   
    if (is_power_two(N) == 0) {
        new_size = compute_next_power_two(N);
    } 
    
    padded_split(A,N, new_size, 0);    

    printf("%d\n", A_11[1][1]);

    printf("A_11\n"); 
    printMatrix(A_11,new_size/2);
    
    printf("A_12\n"); 
    printMatrix(A_12,new_size/2);

    matrix_add(&A_11, A_12, new_size/2);

    printf("A_11\n"); 
    printMatrix(A_11,new_size/2);

/*   printMatrix(A,N);
    
    
    printf("A_12\n"); 
    printMatrix(A_12,new_size/2);
    
    printf("A_21\n"); 
    printMatrix(A_21,new_size/2);
    
    printf("A_22\n"); 
    printMatrix(A_22,new_size/2);
*/

    padded_split(B,N,new_size,1);
/*    printMatrix(B,N);
    
    printf("B_11\n"); 
    printMatrix(B_11,new_size/2);
    
    printf("B_12\n"); 
    printMatrix(B_12,new_size/2);
    
    printf("B_21\n"); 
    printMatrix(B_21,new_size/2);
    
    printf("B_22\n"); 
    printMatrix(B_22,new_size/2);
*/

}

// Allocate square matrix.
int **allocMatrix(int size) {
  int **matrix;
  matrix = (int **)malloc(size * sizeof(int *));
  for (int row = 0; row < size; row++) {
    matrix[row] = (int *)malloc(size * sizeof(int));
  }
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      matrix[i][j] = 0;
    }
  }
  return matrix;
}

// Allocate memory for all the matrixes, you will need to add code
// here to initialize any matrixes that you need.
void initMatrixes(int N) {
  A = allocMatrix(N);
  B = allocMatrix(N);
  C = allocMatrix(N);
  R = allocMatrix(N);
}

// Free up matrixes.
void cleanup() {
  free(A);
  free(B);
  free(C);
  free(R);
}

// Main method
int main(int argc, char* argv[]) {
  int N;
  double elapsedTime;

  // checking parameters
  if (argc != 2 && argc != 4) {
    printf("Parameters: <N> [<fileA> <fileB>]\n");
    return 1;
  }
  N = atoi(argv[1]);
  initMatrixes(N);

  // reading files (optional)
  if(argc == 4){
    readMatrixFile(A,N,argv[2]);
    readMatrixFile(B,N,argv[3]);
  } else {
    // Otherwise, generate two random matrix.
    for (int i=0; i<N; i++) {
      for (int j=0; j<N; j++) {
	A[i][j] = rand() % 5;
	B[i][j] = rand() % 5;
      }
    }
  }

  // Do simple multiplication and time it.
  timerStart();
  simpleMM(N);
  printf("Simple MM took %ld ms\n", timerStop());

  // Do strassen multiplication and time it.
  timerStart();
  strassenMM(N);
  printf("Strassen MM took %ld ms\n", timerStop());

/* NOTE: Commented out for debugging purposes
  if (compareMatrix(C, R, N) != 0) {
    if (N < 20) {
      printf("\n\n------- MATRIX C\n");
      printMatrix(C,N);
      printf("\n------- MATRIX R\n");
      printMatrix(R,N);
    }
    printf("Matrix C doesn't match Matrix R, if N < 20 they will be printed above.\n");
  }
*/
  // stopping timer
  
  cleanup();
  return 0;
}
