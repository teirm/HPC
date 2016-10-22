/**
 * University of Pittsburgh
 * Department of Computer Science
 * CS1645/CS2045: Introduction to HPC Systems
 * Instructor: Xiaolong Cui 
 * STUDENTS: Implement OpenMP parallel shear sort.
 */

#include <omp.h>
#include <math.h>
#include "timer.h"
#include "io.h"
#include "shear.h"

#include <algorithm>
#include <vector>
#include <iostream>

#define MAX_VALUE 10000

using namespace std;

void shear_sort(int **A, int M) {

    int shear_iterations;

    shear_iterations = log2(M*M);

    for (int i = 0; i < shear_iterations; i++) {
        
#       pragma omp parallel for num_threads(M) 
        for (int j = 0; j < M; j++) {
            sort_row(&A, j, M);
        }

#       pragma omp parallel for num_threads(M) 
        for (int k = 0; k < M; k++) {
            sort_col(&A, k, M);
        }
    }
}

void sort_row(int ***M, int row_num, int size)
{
    vector<int> row_vector;
    
    for (int i = 0; i < size; i++) {
        row_vector.push_back((*M)[row_num][i]);
    }

    if (row_num % 2 == 0) {
        sort(row_vector.begin(), row_vector.end());
    } else {
        sort(row_vector.begin(), row_vector.end(), greater<int>());
    }

#       pragma omp parallel for num_threads(size) 
    for (int i = 0; i < size; i++) {
        (*M)[row_num][i] = row_vector[i];
    }
}

void sort_col(int ***M, int col_num, int size)
{
    vector<int> col_vector;

    for (int i = 0; i < size; i++) {
        col_vector.push_back((*M)[i][col_num]);
    }
    sort(col_vector.begin(), col_vector.end());

#       pragma omp parallel for num_threads(size) 
    for (int i = 0; i < size; i++) {
        (*M)[i][col_num] = col_vector[i];
    }
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

void matrix_free(int ***M, int dim)
{
    int i;

    for (i = 0; i < dim; i++) {
        free((*M)[i]);
    }

    free(*M);
}

// Main method      
int main(int argc, char* argv[]) {
  int N, M;
  int **A;
  double elapsedTime;

  // checking parameters
  if (argc != 2 && argc != 3) {
    printf("Parameters: <N> [<file>]\n");
    return 1;
  }
  N = atoi(argv[1]);
  M = (int) sqrt(N); 
  if(N != M*M){
    printf("N has to be a perfect square!\n");
    exit(1);
  }

  // allocating matrix A
  A = allocMatrix(M);

  // reading files (optional)
  if(argc == 3){
    readMatrixFile(A,M,argv[2]);
  } else {
    srand (time(NULL));
    // Otherwise, generate random matrix.
    for (int i=0; i<M; i++) {
      for (int j=0; j<M; j++) {
	A[i][j] = rand() % MAX_VALUE;
      }
    }
  }
  
  // starting timer
  timerStart();

  // calling shear sort function
  shear_sort(A,M);
  // stopping timer
  elapsedTime = timerStop();

  // print if reasonably small

  if (M <= 10) {
    printMatrix(A,M);
  }

  printf("Took %ld ms\n", timerStop());

  // releasing memory
  for (int i=0; i<M; i++) {
    delete [] A[i];
  }
  delete [] A;

  return 0;
}
