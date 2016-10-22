/**
 * Header file for life.cc
 * Author: Cyrus Ramavarapu
 * Date 11 October 2016
 */

void conway(int **World, int N, int M);

/* Matrix Functions */
int **allocMatrix(int size);

void matrix_free(int ***M, int dim);
void copy_matrix(int ***S, int ***D, int size);
int process_cell(int x_pos, int y_pos, int size, int **world);
