/**
 * University of Pittsburgh
 * Department of Computer Science
 * CS1645/CS2045: Introduction to HPC Systems
 * Instructor: Xiaolong Cui
 * Students:
 * Implement openmp verions of conway's game of life.
 */

#include <iostream>
#include "timer.h"
#include "io.h"
#include "life.h"

#define TOTAL_NEIGHBORS 8

using namespace std;

// Function implementing Conway's Game of Life
void conway(int **World, int N, int M){
  // STUDENT: IMPLEMENT THE GAME HERE, make it parallel!
    int **old_world;

    old_world = allocMatrix(N);            
    copy_matrix(&World, &old_world, N);        
    
    printf("World:\n");
    printMatrix(World, N);
   
    printf("\n\nOld_World:\n");
    printMatrix(old_world, N);

    cout << endl << "STARTING EVOLUTION" << endl;

    /* Outerloop to evolve system */
    for (int gen = 0; gen < M; gen++) {
        /* Innerloop to scan and check system */
#       pragma omp parallel for num_threads(N)         
        for (int i = 0; i < N; i++) {
#       pragma omp parallel for num_threads(N)
            for (int j = 0; j < N; j++) {
                World[i][j] = process_cell(i,j,N,old_world);  
            }
        }
       
        copy_matrix(&World, &old_world, N);  
        cout << "EVOVLED WORLD GENERATION: "  << gen << endl;
        printMatrix(World, N); 
    }
}

int process_cell(int x_pos, int y_pos, int size, int **world)
{
    int current_cell; 
     
    int left_neighbor;
    int right_neighbor;
    int top_neighbor;
    int bottom_neighbor;
   
    int top_left_neighbor;
    int top_right_neighbor;    
    int bottom_left_neighbor;
    int bottom_right_neighbor; 

    int live_neighbors;
    int dead_neighbors;
 
    left_neighbor = 0;
    right_neighbor = 0;
    top_neighbor = 0;
    bottom_neighbor = 0; 
  
    top_left_neighbor = 0;
    top_right_neighbor = 0;
    bottom_right_neighbor = 0;
    bottom_left_neighbor = 0;  
 
    current_cell = world[x_pos][y_pos]; 
    
    live_neighbors = 0;
    dead_neighbors = 0; 
/*
    printf("Current Position: [%d,%d]\n", x_pos, y_pos); 
*/
    if (x_pos - 1 >= 0) {
        right_neighbor = world[x_pos - 1][y_pos];
         
        if (y_pos - 1 >= 0) {
            top_right_neighbor = world[x_pos-1][y_pos-1];
        }
    
        if (y_pos + 1 < size) {
            bottom_right_neighbor = world[x_pos-1][y_pos+1];
        }

    }

    if (y_pos - 1 >= 0) {
        top_neighbor = world[x_pos][y_pos-1]; 
    }
    
    if (y_pos + 1 < size) {
        bottom_neighbor = world[x_pos][y_pos+1];
    }
    
    if (x_pos + 1 < size) {
        left_neighbor = world[x_pos+1][y_pos];
        
        if (y_pos - 1 >= 0) {
            top_left_neighbor = world[x_pos+1][y_pos-1];
        }
    
        if (y_pos + 1 < size) {
            bottom_left_neighbor = world[x_pos+1][y_pos+1];
        }
    }
/*   
    printf("Top: %d\nBottom: %d\nRight: %d\nLeft:%d\n",\
            top_neighbor, bottom_neighbor, right_neighbor, left_neighbor);
    printf("Top Left: %d\nBottom Left: %d\nTop Right: %d\nBottom Right:%d\n",\
            top_left_neighbor, bottom_left_neighbor, top_right_neighbor, bottom_right_neighbor);
*/   
    live_neighbors = right_neighbor + left_neighbor +\
                     top_neighbor + bottom_neighbor +\
                     top_right_neighbor + top_left_neighbor +\
                     bottom_right_neighbor + bottom_left_neighbor;
   
    dead_neighbors = TOTAL_NEIGHBORS - live_neighbors; 

    if (current_cell == 1) {
        if (live_neighbors < 2 || live_neighbors > 3) {
            current_cell = 0;
        } else {
            current_cell = 1;
        }
    } else {
        if (live_neighbors == 3) {
            current_cell = 1;
        }
    }
/*
    printf("Dead Neighbors: %d\nLive Neighbors:%d\nResult: %d\n\n",\
            dead_neighbors, live_neighbors, current_cell); 
*/
    return current_cell;
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

void copy_matrix(int ***S, int ***D, int size)
{
#   pragma omp parallel for num_threads(size)  
    for (int i = 0; i < size; i++) {
#   pragma omp parallel for num_threads(size)
        for (int j = 0; j < size; j++) {
            (*D)[i][j] = (*S)[i][j];
        }
    }
}

// Main method      
int main(int argc, char* argv[]) {
  int N,M;
  int **World;
  double elapsedTime;

  // checking parameters
  if (argc != 3 && argc != 4) {
    printf("Parameters: <N> <M> [<file>]\n");
    return 1;
  }
  N = atoi(argv[1]);
  M = atoi(argv[2]);

  // allocating matrices
  World = allocMatrix(N);

  // reading files (optional)
  if(argc == 4){
    readMatrixFile(World,N,argv[3]);
  } else {
    // Otherwise, generate two random matrix.
    srand (time(NULL));
    for (int i=0; i<N; i++) {
      for (int j=0; j<N; j++) {
	World[i][j] = rand() % 2;
      }
    }
  }

  // starting timer
  timerStart();

  // calling conway's game of life 
  conway(World,N,M);

  // stopping timer
  elapsedTime = timerStop();

//  printMatrix(World,N);

  printf("Took %ld ms\n", timerStop());

  // releasing memory
  for (int i=0; i<N; i++) {
    delete [] World[i];
  }
  delete [] World;

  return 0;
}
