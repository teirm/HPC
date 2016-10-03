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

/* C Matrices needed for Strassen Algorithm */
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

void matrix_free(int ***M, int dim)
{
    int i;

    for (i = 0; i < dim; i++) {
        free((*M)[i]);
    }

    free(*M);
}


void matrix_mult(int **M_1, int **M_2, int ***D, int dim)
{
    int i;
    int j;
    int k;

    for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            for (k = 0; k < dim; k++) {
                (*D)[i][j] += M_1[i][k] * M_2[k][j];
            }
        }
    }
}

void matrix_add(int **T_1, int **T_2, int ***D, int dim)
{
    int i;
    int j;

    for(i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            (*D)[i][j] = T_1[i][j] + T_2[i][j];
        }
    }
}

void matrix_sub(int **T_1, int **T_2, int ***D, int dim)
{
    int i;
    int j;

    for(i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            (*D)[i][j] = T_1[i][j] - T_2[i][j];
        }
    }
}

void strassen_allocate(int dim)
{
    M_1 = allocMatrix(dim);
    M_2 = allocMatrix(dim);
    M_3 = allocMatrix(dim);
    M_4 = allocMatrix(dim);
    M_5 = allocMatrix(dim);
    M_6 = allocMatrix(dim);
    M_7 = allocMatrix(dim);

    C_11 = allocMatrix(dim);
    C_12 = allocMatrix(dim);
    C_21 = allocMatrix(dim);
    C_22 = allocMatrix(dim);
}

void strassen_deallocate(int dim)
{
    matrix_free(&M_1, dim);
    matrix_free(&M_2, dim);
    matrix_free(&M_3, dim);
    matrix_free(&M_4, dim);
    matrix_free(&M_5, dim);
    matrix_free(&M_6, dim);
    matrix_free(&M_7, dim);
    
    matrix_free(&C_11, dim);
    matrix_free(&C_12, dim);
    matrix_free(&C_21, dim);
    matrix_free(&C_22, dim);            
}


void recombine_matrices(int **C_11, int **C_12, int **C_21, int **C_22, int ***D, int old_dim, int new_dim)
{
    int i;
    int j;

    printf("Going to recombine matrices\n");

    for (i = 0; i < old_dim; i++) {
        for (j = 0; j < old_dim; j++) {
            if (i < new_dim/2 && j < new_dim/2) { 
//                printf("Adding C_11 at %d %d\n", i,j);
                (*D)[i][j] = C_11[i][j];
            } else if (i >= new_dim/2 && j < new_dim/2) {
//                printf("Adding C_21 at %d %d\n", i-new_dim/2,j);
                (*D)[i][j] = C_21[i-new_dim/2][j];
            } else if (i < new_dim/2 && j >= new_dim/2) {
//                printf("Adding C_12 at %d %d\n", i,j-new_dim/2);
                (*D)[i][j] = C_12[i][j-new_dim/2];
            } else {
//                printf("Adding C_22 at %d %d\n", i-new_dim/2,j-new_dim/2);
                (*D)[i][j] = C_22[i-new_dim/2][j-new_dim/2];
            }
        }
    }
}



void calc_M1(int **A_11, int **A_22, int **B_11, int **B_22, int ***D, int dim)
{

    int **T_1;
    int **T_2;

    T_1 = allocMatrix(dim);
    T_2 = allocMatrix(dim);
    
    matrix_add(A_11, A_22, &T_1, dim);
    matrix_add(B_11, B_22, &T_2, dim);

    matrix_mult(T_1, T_2, D, dim);

    matrix_free(&T_1, dim);
    matrix_free(&T_2, dim);
}

void calc_M2(int **A_21, int **A_22, int **B_11, int ***D, int dim)
{
    
    int **T_1;
    
    T_1 = allocMatrix(dim);
    
    matrix_add(A_21, A_22, &T_1, dim);
    matrix_mult(T_1, B_11, D, dim); 
    
    matrix_free(&T_1, dim);
}

void calc_M3(int **A_11, int **B_12, int **B_22, int ***D, int dim)
{

    int **T_1;

    T_1 = allocMatrix(dim);

    matrix_sub(B_12, B_22, &T_1, dim);
    matrix_mult(A_11, T_1, D, dim);
    
    matrix_free(&T_1, dim);
}

void calc_M4(int **A_22, int **B_21, int **B_11, int ***D, int dim)
{
    
    int **T_1;
    
    T_1 = allocMatrix(dim);
    
    matrix_sub(B_21, B_11, &T_1, dim);
    matrix_mult(A_22, T_1, D, dim);    
    
    matrix_free(&T_1, dim);
}


void calc_M5(int **A_11, int **A_12, int **B_22, int ***D, int dim)
{
    int **T_1;
    
    T_1 = allocMatrix(dim);
    
    matrix_add(A_11, A_12, &T_1, dim);
    matrix_mult(T_1, B_22, D, dim);
    
    matrix_free(&T_1, dim);
}

void calc_M6(int **A_21, int **A_11, int **B_11, int **B_12, int ***D, int dim)
{
    int **T_1;
    int **T_2;
    
    T_1 = allocMatrix(dim);
    T_2 = allocMatrix(dim); 
    
    matrix_sub(A_21, A_11, &T_1, dim);
    matrix_add(B_11, B_12, &T_2, dim);
    matrix_mult(T_1, T_2, D, dim); 
    
    matrix_free(&T_1, dim);
    matrix_free(&T_2, dim);
}

void calc_M7(int **A_12, int **A_22, int **B_21, int **B_22, int ***D, int dim)
{
    int **T_1;
    int **T_2;

    T_1 = allocMatrix(dim);
    T_2 = allocMatrix(dim);

    matrix_sub(A_12, A_22, &T_1, dim);
    matrix_add(B_21, B_22, &T_2, dim);
    matrix_mult(T_1, T_2, D, dim);
    
    matrix_free(&T_1, dim);
    matrix_free(&T_2, dim);
}

void calc_C11(int **M_1, int **M_4, int **M_5, int **M_7, int ***D, int dim)
{
    matrix_add(M_1, M_4, D, dim);
    matrix_sub(*D, M_5, D, dim);
    matrix_add(*D, M_7, D, dim);    
}

void calc_C12(int **M_3, int **M_5, int ***D, int dim)
{
    matrix_add(M_3, M_5, D, dim);    
}
void calc_C21(int **M_2, int **M_4, int ***D, int dim)
{
    matrix_add(M_2, M_4, D, dim);    
}

void calc_C22(int **M_1, int **M_2, int **M_3, int **M_6, int ***D, int dim)
{
    matrix_sub(M_1, M_2, D, dim);
    matrix_add(*D, M_3, D, dim);
    matrix_add(*D, M_6, D, dim);    
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
    
    strassen_allocate(new_size/2); 
  
    printf("Matrix A:\n");
    printMatrix(A,N);
     
    printf("Matrix B:\n");
    printMatrix(B,N);
   
    padded_split(A,N,new_size,0);
    padded_split(B,N,new_size,1);
   
    calc_M1(A_11, A_22, B_11, B_22, &M_1, new_size/2); 
    calc_M2(A_21, A_22, B_11, &M_2, new_size/2);
    calc_M3(A_11, B_12, B_22, &M_3, new_size/2);
    calc_M4(A_22, B_21, B_11, &M_4, new_size/2);
    calc_M5(A_11, A_12, B_22, &M_5, new_size/2);
    calc_M6(A_21, A_11, B_11, B_12, &M_6, new_size/2);
    calc_M7(A_12, A_22, B_21, B_22, &M_7, new_size/2);

    calc_C11(M_1, M_4, M_5, M_7, &C_11, new_size/2);
    calc_C12(M_3, M_5, &C_12, new_size/2);
    calc_C21(M_2, M_4, &C_21, new_size/2);
    calc_C22(M_1, M_2, M_3, M_6, &C_22, new_size/2);
    
    recombine_matrices(C_11, C_12, C_21, C_22, &C, N, new_size);

    printMatrix(C, N);
    
    strassen_deallocate(new_size/2);
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
void cleanup(int dim) {
  matrix_free(&A, dim);
  matrix_free(&B, dim);
  matrix_free(&C, dim);
  matrix_free(&R, dim);
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

  if (compareMatrix(C, R, N) != 0) {
    if (N < 20) {
      printf("\n\n------- MATRIX C\n");
      printMatrix(C,N);
      printf("\n------- MATRIX R\n");
      printMatrix(R,N);
    }
    printf("Matrix C doesn't match Matrix R, if N < 20 they will be printed above.\n");
  }
  // stopping timer
  
  cleanup(N);
  return 0;
}
