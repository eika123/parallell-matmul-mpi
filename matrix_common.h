
#include <stdlib.h>

#ifndef MAT_COMMON_H
#define MAT_COMMON_H

#define SELECT_DATA_TYPE 2

// enums only work compile time
#define SHORT 0
#define INT 1
#define FLOAT 2
#define DOUBLE 3

#if SELECT_DATA_TYPE == SHORT
typedef short DATA_TYPE;
#elif SELECT_DATA_TYPE == INT
typedef int DATA_TYPE;
#elif SELECT_DATA_TYPE == FLOAT
typedef float DATA_TYPE;
#elif SELECT_DATA_TYPE == DOUBLE
typedef double DATA_TYPE;
#endif


// trickery: store matrix as a 1D array
struct matrix
{
    int numRows; // height of column = the number of rows
    int numCols; // lenght of rows = the number of columsn
    DATA_TYPE *data;
} typedef matrix_t;

#define get_val(matrix, row, col) ((matrix->data)[row * (matrix->numCols) + col])
#define set_val(matrix, row, col, value) ((matrix->data)[row * (matrix->numCols) + col] = (DATA_TYPE)value)


matrix_t *matrix_factory(int m, int n);

void populate_vector_constants(matrix_t *vec, DATA_TYPE val);
void populate_vector_random(matrix_t *vec, int range);
void populate_matrix_random(matrix_t *mat, int range);


#endif