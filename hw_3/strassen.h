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

void matrix_add(int ***M_1, int **M_2, int dim);
void matrix_sub(int ***M_1, int **M_2, int dim);
