/**
 * Author:  Cyrus Ramavarapu
 * Purpose: Function declarations for additions
 *          to strassen.c
 * Date:    30 September 2016
 */

/* Algebraic Functions */ 

int is_power_two(int val);

int compute_next_power_two(int val);

/* Matrix Functions */

void padded_split(int **M, int N, int padded_size, int m_pos);

int **allocMatrix(int *size);

void matrix_add(int **M_1, int **M_2, int ***D, int *dim);
void matrix_sub(int **M_1, int **M_2, int ***D, int *dim);
void matrix_mult(int **M_1, int **M_2, int ***D, int *dim);

void recombine_matrices(int **C_11, int **C_12, int **C_21, int **C_22, int ***D, int old_dim, int new_dim);

void matrix_free(int ***M, int *dim);

void cleanup(int *dim);

/* Functions for Strassens Algorithm */
void strassen_allocate(int *dim);
void strassen_deallocate(int *dim);

void *calc_M1(void *dim);
void *calc_M2(void *dim);
void *calc_M3(void *dim);
void *calc_M4(void *dim);
void *calc_M5(void *dim);
void *calc_M6(void *dim);
void *calc_M7(void *dim);

void *calc_C11(void *dim);
void *calc_C12(void *dim);
void *calc_C21(void *dim);
void *calc_C22(void *dim);


