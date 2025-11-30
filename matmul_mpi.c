#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>

#define SELECT_DATA_TYPE 2
typedef enum number_types
{
    SHORT = 0,
    INT,
    FLOAT,
    DOUBLE
} NUMBER_TYPE;

#if SELECT_DATA_TYPE == 0
typedef short DATA_TYPE;
#elif SELECT_DATA_TYPE == 1
typedef int DATA_TYPE;
#elif SELECT_DATA_TYPE == 2
typedef float DATA_TYPE;
#elif SELECT_DATA_TYPE == 3
typedef double DATA_TYPE;
#endif

// trickery: store matrix as a 1D array
struct matrix
{
    int numRows; // height of column = the number of rows
    int numCols; // lenght of rows = the number of columsn
    DATA_TYPE *data;
} typedef matrix_t;

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

#define get_val(matrix, row, col) ((matrix->data)[row * (matrix->numCols) + col])
#define set_val(matrix, row, col, value) ((matrix->data)[row * (matrix->numCols) + col] = (DATA_TYPE)value)

void raw_pprint_matrix(matrix_t *A, char *row_container_left, char *row_container_right, char *elem_delim, char *outer_container_left, char *outer_container_right)
{
    int N = A->numRows, M = A->numCols, i = 0, j = 0;
    if (outer_container_left) {
        printf("%s\n", outer_container_left);
    }
    for (int i = 0; i < N; i++)
    {
        printf("%s", row_container_left);
        for (j = 0; j < M; j++)
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

void pprint_matrix(matrix_t *A)
{
    raw_pprint_matrix(A, "|", "|", NULL, NULL, NULL);
}

void python_pprint_matrix(matrix_t *A)
{
    raw_pprint_matrix(A, "[", "],", ",", "[", "]");
}

inline void sanitize_serial_matvecmul_input(matrix_t *mat, matrix_t *vec, matrix_t *result)
{
    if (!mat || !vec || !result)
    {
        fprintf(stderr, "serial_matvecmul expects args to be not NULL, received "
                        "mat=%p, vec=%p, result=%p\n",
                mat, vec, result);
        exit(1);
    }
    if (vec->numCols != 1)
    {
        fprintf(stderr, "serial_matvecmul expects vector to have one column only");
        exit(1);
    }
    if (vec->numRows != mat->numCols)
    {
        fprintf(stderr, "matrix and vector dimensions mismatch for multiplication, "
                        "matrix has %d columns, vector has %d rows\n",
                vec->numRows,
                mat->numCols);
        exit(1);
    }
}

int serial_matvecmul(matrix_t *mat, matrix_t *vec, matrix_t *result)
{
    sanitize_serial_matvecmul_input(mat, vec, result);
    int M = mat->numRows, N = mat->numCols;
    DATA_TYPE mres;
    for (int row = 0; row < M; row++)
    {
        mres = (DATA_TYPE)0;
        for (int col = 0; col < N; col++)
        {
            mres += get_val(mat, row, col) * get_val(vec, col, 0);
        }
        set_val(result, row, 0, mres);
    }
}

void populate_vector_constants(matrix_t *vec, DATA_TYPE val)
{
    int M = vec->numRows;

    for (int row = 0; row < M; row++)
    {
        set_val(vec, row, 0, val);
    }
}

void populate_vector_random(matrix_t *vec, int range)
{
    int M = vec->numRows;

    for (int row = 0; row < M; row++)
    {
        set_val(vec, row, 0, (DATA_TYPE)(rand() % range));
    }
}

void populate_matrix_random(matrix_t *mat, int range)
{
    int M = mat->numRows, N = mat->numCols;
    for (int row = 0; row < M; row++)
    {
        for (int col = 0; col < N; col++)
        {
            set_val(mat, row, col, (DATA_TYPE)(rand() % range));
        }
    }
}

int main(int argc, char *argv[])
{
    matrix_t *mat = matrix_factory(20, 10);
    matrix_t *vec = matrix_factory(10, 1);
    matrix_t *res = matrix_factory(20, 1);

    populate_matrix_random(mat, 100);
    populate_vector_random(vec, 20);

    serial_matvecmul(mat, vec, res);

    python_pprint_matrix(mat);
    printf("\n");
    printf("\n");
    python_pprint_matrix(vec);
    printf("\n");
    printf("\n");
    python_pprint_matrix(res);
}