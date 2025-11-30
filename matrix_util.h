#include <stddef.h>

#include "matrix_common.h"

#ifndef MATRIX_UTIL
#define MATRIX_UTIL

void raw_pprint_matrix(matrix_t *A, char *row_container_left, char *row_container_right, char *elem_delim, char *outer_container_left, char *outer_container_right);

void pprint_matrix(matrix_t *A)
{
    raw_pprint_matrix(A, "|", "|", NULL, NULL, NULL);
}

void python_pprint_matrix(matrix_t *A)
{
    raw_pprint_matrix(A, "[", "],", ",", "[", "]");
}

#endif