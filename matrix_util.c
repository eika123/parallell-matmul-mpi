#include <stdio.h>

#ifndef MAT_COMMON_H
#include "matrix_common.h"
#endif


/**
 * @brief raw pretty printing
 * 
 * @param A matrix_t: matrix to be printed
 * @param row_container_left left delimiter of rows, typically [ or |
 * @param row_container_right right delimiter of rows, typically ] or |
 * @param elem_delim element delimiter, typically comma, whitespace or a tab
 * @param outer_container_left
 * @param outer_container_right 
 */
void raw_pprint_matrix(matrix_t *A, char *row_container_left, char *row_container_right, char *elem_delim, char *outer_container_left, char *outer_container_right)
{
    int M = A->numRows, N = A->numCols, i = 0, j = 0;
    if (outer_container_left) {
        printf("%s\n", outer_container_left);
    }
    for (int i = 0; i < M; i++)
    {
        printf("%s", row_container_left);
        for (j = 0; j < N; j++)
        {
            switch (SELECT_DATA_TYPE)
            {
            case SHORT:
                printf(" %8hd", (short)get_val(A, i, j));
                break;
            case INT:
                printf(" %8d", (int)get_val(A, i, j));
                break;
            case FLOAT:
                printf(" %8.4f", (float)get_val(A, i, j));
                break;
            case DOUBLE:
                printf(" %8.4f", (double)get_val(A, i, j));
                break;
            default:
                break;
            }
            if (elem_delim && j < M - 1)
            {
                printf("%s", elem_delim);
            }
        }
        printf(" %s\n", row_container_right);
    }
    if (outer_container_right) {
        printf("%s", outer_container_right);
    }
}
