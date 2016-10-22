/**
 * Header file for shear.cc
 * Author: Cyrus Ramavarapu
 * Date: 11 October 16
 */

void shear_sort(int **A, int M);

int **allocMatrix(int Size);

void matrix_free(int ***M, int dim);

void sort_row(int ***M, int row_num, int size);
void sort_col(int ***M, int col_num, int size);
