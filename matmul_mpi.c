#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "matrix_common.h"
#include "matrix_util.h"

#if SELECT_DATA_TYPE == SHORT
#    define SELECT_DATA_TYPE_MPI MPI_SHORT
#elif SELECT_DATA_TYPE == INT
#    define SELECT_DATA_TYPE_MPI MPI_INT
#elif SELECT_DATA_TYPE == FLOAT
#    define SELECT_DATA_TYPE_MPI MPI_FLOAT
#elif SELECT_DATA_TYPE == DOUBLE
#    define SELECT_DATA_TYPE_MPI MPI_DOUBLE
#endif

inline void sanitize_serial_matvecmul_input(matrix_t *mat,
                                            matrix_t *vec,
                                            matrix_t *result)
{
    if (!mat || !vec || !result) {
        fprintf(stderr,
                "serial_matvecmul expects args to be not NULL, received "
                "mat=%p, vec=%p, result=%p\n",
                mat,
                vec,
                result);
        exit(1);
    }
    if (vec->numCols != 1) {
        fprintf(stderr, "serial_matvecmul expects vector to have one column only");
        exit(1);
    }
    if (vec->numRows != mat->numCols) {
        fprintf(stderr,
                "matrix and vector dimensions mismatch for multiplication, "
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
    for (int row = 0; row < M; row++) {
        mres = (DATA_TYPE) 0;
        for (int col = 0; col < N; col++) {
            mres += get_val(mat, row, col) * get_val(vec, col, 0);
        }
        set_val(result, row, 0, mres);
    }
}

/**
 * @brief
 *
 * @param send_matrix: matrix_t pointer to matrix that is to be scattered
 * @param sendcounts: (Root only) Array where sendcounts[i] is the number of elements to
 *  send to process i.
 * @param dspls: (Root only) Array where displs[i] is the offset from the beginning of
 *  sendbuf for data intended for process i
 * @param recv_matrix
 * @param root
 * @param comm MPI_Comm communicator in which to scatter the send_matrix
 * @return int
 */
int scatterv_matrix_mpi(matrix_t *send_matrix,
                        const int sendcounts[],
                        const int dspls[],
                        matrix_t *recv_matrix,
                        int root,
                        MPI_Comm comm)
{
    // rank is specific to communicator
    int comm_size, rank;
    MPI_Comm_size(comm, &comm_size);
    MPI_Comm_rank(comm, &rank);

    if (rank == root && send_matrix == NULL) {
        fprintf(stderr, "Rank[%d]: root sender received NULL as send_matrix\n", rank);
        exit(1);
    }

    DATA_TYPE *sendbuf = rank == root ? (DATA_TYPE *) send_matrix->data : NULL;

    int recv_count = (recv_matrix->numRows) * (recv_matrix->numCols);
    MPI_Scatterv(send_matrix->data,
                 sendcounts,
                 dspls,
                 SELECT_DATA_TYPE_MPI,
                 recv_matrix->data,
                 recv_count,
                 SELECT_DATA_TYPE_MPI,
                 root,
                 comm);
}

/**
 * Wrapper for gathering a matrix with MPI using MPI_Gatherv
 */
int gatherv_matrix_mpi(const matrix_t *send_matrix,
                       int sendcount,
                       matrix_t *recv_matrix,
                       const int recvcount[],
                       const int displs[],
                       int root,
                       MPI_Comm comm,
                       MPI_Request *request)
{
    return 0;
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

    return 0;
}