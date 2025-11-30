#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>

typedef enum
{
    SHORT,
    INT,
    FLOAT,
    DOUBLE,
} NUMBER_TYPE;

size_t numberSize(NUMBER_TYPE num_type)
{
    switch (num_type)
    {
    case SHORT:
        return sizeof(short);
    case INT:
        return sizeof(int);
    case FLOAT:
        return sizeof(float);
    case DOUBLE:
        return sizeof(double);
    default:
        fprintf(stderr, "unknown number type\n");
        exit(1);
    }
}

// trickery: store matrix as a long array
struct matrix
{
    int numRows; // height of column = the number of rows
    int numCols; // lenght of rows = the number of columsn
    NUMBER_TYPE dtype;
    void *data; //
} typedef matrix_t;

matrix_t *matrix_factory(int m, int n, NUMBER_TYPE num_type)
{
    matrix_t *matrix = malloc(sizeof(matrix_t));
    if (!matrix)
    {
        fprintf(stderr, "malloc failed\n");
        exit(1);
    }
    matrix->numRows = m;
    matrix->numCols = n;
    matrix->dtype = num_type;
    matrix->data = calloc((size_t)m * (size_t)n, numberSize(num_type));
    if (!matrix->data)
    {
        fprintf(stderr, "calloc failed\n");
        free(matrix);
    }
    return matrix;
}

#define get_val(matrix, row, col, dtype) (((dtype *)matrix->data)[row * (matrix->numCols) + col])
#define set_val(matrix, row, col, value, dtype) (((dtype *)matrix->data)[row * (matrix->numCols) + col] = (dtype)value)

void pprint_matrix(matrix_t *A)
{
    int N = A->numRows;
    int M = A->numCols;
    int i = 0, j = 0;
    for (int i = 0; i < N; i++)
    {
        printf("|");
        for (j = 0; j < M; j++)
        {
            switch (A->dtype)
            {
            case SHORT:
                printf(" %8hd", get_val(A, i, j, short));
                break;
            case INT:
                printf(" %8hd", get_val(A, i, j, int));
                break;
            case FLOAT:
                printf(" %8.4f", get_val(A, i, j, float));
                break;
            case DOUBLE:
                printf(" %8.4f", get_val(A, i, j, double));
                break;
            default:
                break;
            }
        }
        printf(" |\n");
    }
}


int main(int argc, char *argv[])
{
    matrix_t *mat_floats = matrix_factory(5, 10, FLOAT);
    matrix_t *mat_ints = matrix_factory(5, 10, INT);
    matrix_t *mat_shorts = matrix_factory(5, 10, SHORT);

    set_val(mat_floats, 1, 1, 1.0, float);
    set_val(mat_ints, 1, 1, 1, int);
    set_val(mat_shorts, 1, 1, 25, short);

    pprint_matrix(mat_floats);
    printf("\n");
    pprint_matrix(mat_ints);
    printf("\n");
    pprint_matrix(mat_shorts);
}