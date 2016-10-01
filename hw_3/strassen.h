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

int **allocMatrix(int size);

void matrix_add(int **M_1, int **M_2, int ***D, int dim);
void matrix_sub(int **M_1, int **M_2, int ***D, int dim);
void matrix_mult(int **M_1, int **M_2, int ***D, int dim);

/* Functions for Strassens Algorithm */
void strassen_allocate(int dim);

void calc_M1(int **A_11, int **A_22, int **B_11, int **B_22, int ***D, int dim);
void calc_M2(int **A_21, int **A_22, int **B_11, int ***D, int dim);
void calc_M3(int **A_11, int **B_12, int **B_22, int ***D, int dim);
void calc_M4(int **A_22, int **B_21, int **B_11, int ***D, int dim);
void calc_M5(int **A_11, int **A_12, int **B_22, int ***D, int dim);
void calc_M6(int **A_21, int **A_11, int **B_11, int **B_12, int ***D, int dim);
void calc_M7(int **A_12, int **A_22, int **B_21, int **B_22, int ***D, int dim);

void calc_C11(int **M_1, int **M_4, int **M_5, int **M_7);
void calc_C12(int **M_3, int **M_5);
void calc_C21(int **M_2, int **M_4);
void calc_C22(int **M_1, int **M_2, int **M_3, int **M_6);
