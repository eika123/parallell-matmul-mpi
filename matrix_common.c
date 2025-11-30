#include <stdlib.h>
#include <stdio.h>
#include "matrix_common.h"

matrix_t *matrix_factory(int m, int n)
{
    matrix_t *matrix = malloc(sizeof(matrix_t));
    if (!matrix)
    {
        fprintf(stderr, "malloc failed\n");
        exit(1);
    }
    matrix->numRows = m;
    matrix->numCols = n;
    matrix->data = calloc((size_t)m * (size_t)n, sizeof(DATA_TYPE));
    if (!matrix->data)
    {
        fprintf(stderr, "calloc failed\n");
        free(matrix);
    }
    return matrix;
}


void populate_vector_constants(matrix_t *vec, DATA_TYPE val) {
    for (int row = 0; row < vec->numRows; row++) {
        set_val(vec, row, 0, val);
    }
}

void populate_vector_random(matrix_t *vec, int range) {
    for (int row = 0; row < (vec->numRows); row++) {
        set_val(vec, row, 0, (rand() % range));
    }
}

void populate_matrix_random(matrix_t *mat, int range) {
    for (int row = 0; row < mat->numRows; row++) {
        for (int col = 0; col < mat->numCols; col++) {
            set_val(mat, row, col, (DATA_TYPE)(rand() % range));
        }
    }
}