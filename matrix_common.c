#include "matrix_common.h"

#include <stdio.h>
#include <stdlib.h>

matrix_t *matrix_factory(int numRows, int numCols)
{
    matrix_t *matrix = malloc(sizeof(matrix_t));
    if (!matrix) {
        printf("malloc failed\n");
        exit(1);
    }
    matrix->numRows = numRows;
    matrix->numCols = numCols;
    matrix->data = calloc((size_t) numRows * (size_t) numCols, sizeof(DATA_TYPE));
    if (!matrix->data) {
        printf("calloc failed to allocate matrix->data\n");
        free(matrix);
        exit(1);
    }
    return matrix;
}

void matrix_destroy(matrix_t *matrix)
{
    if (matrix == NULL) {
        printf(
            "matrix_common::matrix_destroy: expected argument matrix to be not NULL\n");
        exit(1);
    }
    if (matrix->data == NULL) {
        printf(
            "matrix_common::matrix_destroy: expected argument matrix to have not NULL "
            "data field\n");
        exit(1);
    }
    free(matrix->data);
    free(matrix);
}

void populate_vector_constants(matrix_t *vec, DATA_TYPE val)
{
    for (int row = 0; row < vec->numRows; row++) {
        set_val(vec, row, 0, val);
    }
}

void populate_vector_random(matrix_t *vec, int range)
{
    for (int row = 0; row < (vec->numRows); row++) {
        set_val(vec, row, 0, (rand() % range));
    }
}

void populate_matrix_random(matrix_t *mat, int range)
{
    for (int row = 0; row < mat->numRows; row++) {
        for (int col = 0; col < mat->numCols; col++) {
            set_val(mat, row, col, (DATA_TYPE) (rand() % range));
        }
    }
}

/**
 * @brief Multiply vector by a scalar value
 *
 * @param vector
 * @param scalar_value
 * @return int
 */
void vecmul_scalar(matrix_t *vector, DATA_TYPE scalar_value)
{
    for (int row = 0; row < vector->numRows; row++) {
        (vector->data[row]) *= scalar_value;
    }
}

/**
 * @brief Multiply vector by inverse of scalar value. Useful when scalar value is not a
 * floating point number or for generic code
 *
 * @param vector
 * @param scalar_value
 * @return int
 */
void vecdivide_scalar(matrix_t *vector, DATA_TYPE scalar_value)
{
    for (int row = 0; row < vector->numRows; row++) {
        (vector->data[row]) /= scalar_value;
    }
}