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
        result->data[row] = mres;
    }
}

DATA_TYPE vector_sum_component_squares(matrix_t *vector)
{
    int M = vector->numRows;
    DATA_TYPE sum_squares = (DATA_TYPE) 0;
    DATA_TYPE component;
    for (int i = 0; i < M; i++) {
        component = vector->data[i];
        sum_squares += component * component;
    }
    return sum_squares;
}

DATA_TYPE total_sum_squares(int sum_squares_components[], int length) {
    DATA_TYPE sum_squares = 0;
    for (int i=0; i < length; i++) {
        sum_squares += (sum_squares_components[i]) * (sum_squares_components[i]);
    }
    return sum_squares;
}

/**
 * @brief Helper for scatter functions. Calculates how much to send to each process in a
 * communicator of comm_size, and stores the information in sendcounts and displacements.
 * See @see MPI_Scatterv
 *
 * @param send_matrix
 * @param sendcounts
 * @param displacements
 * @param comm_size
 */
inline void plan_even_distribution(matrix_t *send_matrix,
                                   matrix_t *recv_matrix,
                                   int sendcounts[],
                                   int displacements[],
                                   int comm_size)
{
    int M = send_matrix->numRows, N = send_matrix->numCols;
    int snd_cnt = (M / comm_size) * N;

    for (int rnk = 0; rnk < comm_size; rnk++) {
        sendcounts[rnk] = (M / comm_size) * N;
        displacements[rnk] = rnk * (M / comm_size) * N;
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
                        int sendcounts[],
                        int displacements[],
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
    int recv_count = (recv_matrix->numRows) * (recv_matrix->numCols);

    DATA_TYPE *sendbuf = rank == root ? (DATA_TYPE *) send_matrix->data : NULL;

    MPI_Scatterv(send_matrix->data,
                 sendcounts,
                 displacements,
                 SELECT_DATA_TYPE_MPI,
                 recv_matrix->data,
                 recv_count,
                 SELECT_DATA_TYPE_MPI,
                 root,
                 comm);
}

/**
 * @brief Responsibility of caller to provide and free sendcounts, dispelcount + matrices.
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
int scatterv_matrix_mpi_world(matrix_t *send_matrix,
                              int sendcounts[],
                              int displacements[],
                              matrix_t *recv_matrix,
                              const int rank,
                              const int root,
                              const int world_size)
{
    if (rank == root && send_matrix == NULL) {
        fprintf(stderr, "Rank[%d]: root sender received NULL as send_matrix\n", rank);
        exit(1);
    }
    if (recv_matrix == NULL) {
        printf("Rank[%d]: expected recv_matrix to be not NULL\n", rank);
        exit(1);
    }

    void *send_data;
    if (rank != root && send_matrix == NULL) {
        send_data = NULL;
    } else {
        send_data = send_matrix->data;
    }
    int recv_count = (recv_matrix->numRows) * (recv_matrix->numCols);

    DATA_TYPE *sendbuf = rank == root ? (DATA_TYPE *) send_matrix->data : NULL;
    MPI_Scatterv(send_data,
                 sendcounts,
                 displacements,
                 SELECT_DATA_TYPE_MPI,
                 recv_matrix->data,
                 recv_count,
                 SELECT_DATA_TYPE_MPI,
                 root,
                 MPI_COMM_WORLD);
}

inline void gatherv_matrix_mpi_world_sanitize_input(matrix_t *recv_matrix,
                                                    int sendcount,
                                                    const matrix_t *send_matrix,
                                                    const int recvcounts[],
                                                    const int recv_displacements[],
                                                    int rank,
                                                    int root,
                                                    MPI_Comm comm)
{
    if (rank == root && recvcounts == NULL) {
        fprintf(stderr,
                "Rank[%d]: matmul_mpi::gatherv_matrix_mpi_world: root process expected "
                "recvcounts to be not NULL\n",
                rank);
        exit(1);
    }
    if (rank == root && recv_displacements == NULL) {
        fprintf(stderr,
                "Rank[%d]: matmul_mpi::gatherv_matrix_mpi_world: root process expected "
                "recv_displacements to be not NULL\n",
                rank);
        exit(1);
    }
    if (rank == root && recv_matrix == NULL) {
        fprintf(stderr,
                "Rank[%d]: matmul_mpi::gatherv_matrix_mpi_world: root process expected "
                "recv_matrix to be not NULL\n",
                rank);
        exit(1);
    }
    if (send_matrix == NULL) {
        fprintf(stderr,
                "Rank[%d]: matmul_mpi::gatherv_matrix_mpi_world: expected send_matrix to "
                "be not NULL\n",
                rank);
        exit(1);
    }
}

/**
 * @brief
 *
 * @param send_matrix The matrix to send
 * @param sendcount The number of elements in the send buffer send_matrix->data
 * @param recv_matrix
 * @param recvcount
 * @param displs
 * @param root
 * @param comm
 * @param request
 * @return int
 */
int gatherv_matrix_mpi_world(matrix_t *recv_matrix,
                             int sendcount,
                             const matrix_t *send_matrix,
                             const int recvcounts[],
                             const int recv_displacements[],
                             int rank,
                             int root,
                             MPI_Comm comm)
{
    gatherv_matrix_mpi_world_sanitize_input(recv_matrix,
                                            sendcount,
                                            send_matrix,
                                            recvcounts,
                                            recv_displacements,
                                            rank,
                                            root,
                                            comm);
    MPI_Gatherv(send_matrix->data,
                sendcount,
                SELECT_DATA_TYPE_MPI,
                recv_matrix->data,
                recvcounts,
                recv_displacements,
                SELECT_DATA_TYPE_MPI,
                root,
                MPI_COMM_WORLD);
    return 0;
}

void simple_check()
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

void allocate_receive_matrix(const int rank,
                             const int comm_size,
                             matrix_t **recv_matrix,
                             int send_matrix_num_rows,
                             int send_matrix_num_cols)
{
    int M = send_matrix_num_rows, N = send_matrix_num_cols;
    if (rank == comm_size - 1) {
        *recv_matrix = matrix_factory(M / comm_size + (M % comm_size), N);
    } else {
        *recv_matrix = matrix_factory(M / comm_size, N);
    }
}

int main(int argc, char *argv[])
{
    int world_size, my_rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    int M = 18, N = 7;

    matrix_t *send_matrix = NULL, *vector = NULL, *result_vector;
    int *sendcounts = NULL, *displacements = NULL;
    DATA_TYPE *vector_norm_contributions = calloc(world_size, sizeof(DATA_TYPE));

    matrix_t *recv_matrix, *recv_vector;
    vector = matrix_factory(N, 1);
    allocate_receive_matrix(my_rank, world_size, &recv_matrix, M, N);
    allocate_receive_matrix(my_rank, world_size, &recv_vector, M, 1);
    allocate_receive_matrix(my_rank, world_size, &result_vector, M, 1);

    if (my_rank == 0) {
        send_matrix = matrix_factory(M, N);

        populate_matrix_random(send_matrix, 5);
        populate_vector_random(vector, 4);

        sendcounts = calloc(world_size, sizeof(int));
        displacements = calloc(world_size, sizeof(int));

        plan_even_distribution(
            send_matrix, recv_matrix, sendcounts, displacements, world_size);
        sendcounts[world_size - 1] += M % world_size;
        scatterv_matrix_mpi_world(
            send_matrix, sendcounts, displacements, recv_matrix, my_rank, 0, world_size);

    } else {
        scatterv_matrix_mpi_world(
            send_matrix, sendcounts, displacements, recv_matrix, my_rank, 0, world_size);
    }

    MPI_Bcast(vector->data, vector->numRows, SELECT_DATA_TYPE_MPI, 0, MPI_COMM_WORLD);
    serial_matvecmul(recv_matrix, vector, result_vector);
    vector_norm_contributions[my_rank] = vector_sum_component_squares(result_vector);

    MPI_Allgather((vector_norm_contributions + my_rank),
                  1,
                  SELECT_DATA_TYPE_MPI,
                  vector_norm_contributions,
                  world_size,
                  SELECT_DATA_TYPE_MPI,
                  MPI_COMM_WORLD);

    DATA_TYPE total_squares = total_sum_squares(vector_norm_contributions, world_size);
    printf("Rank[%d]: total_squares=%d\n", my_rank, total_squares);
    

    // MPI_Allgather(result_vector->data,
    //               result_vector->numRows,
    //               SELECT_DATA_TYPE_MPI,
    //               vector->data,
    //               vector->numRows,
    //               SELECT_DATA_TYPE_MPI,
    //               MPI_COMM_WORLD);

    // printf("Rank[%d]: result_vector_length=%d ... result sum squares=%d\n",
    //        my_rank,
    //        result_vector->numRows,
    //        vector_norm_contributions[my_rank]);

    // MPI_Barrier(MPI_COMM_WORLD);
    //  if (my_rank == 0) {
    //      printf("Rank[%d]: M / world_size=%d\n", my_rank, world_size);
    //      pprint_matrix(vector);
    //      printf("\n\n");
    //  }
    //  MPI_Barrier(MPI_COMM_WORLD);
    //  if (my_rank == 1) {
    //      printf("Rank[%d]: M / world_size=%d\n", my_rank, world_size);
    //      pprint_matrix(vector);
    //      printf("\n\n");
    //  }
    //  MPI_Barrier(MPI_COMM_WORLD);
    //  if (my_rank == 2) {
    //      printf("Rank[%d]: M / world_size=%d\n", my_rank, world_size);
    //      pprint_matrix(vector);
    //      printf("\n\n");
    //  }
    //  MPI_Barrier(MPI_COMM_WORLD);
    //  if (my_rank == 3) {
    //      printf("Rank[%d]: M / world_size=%d\n", my_rank, world_size);
    //      pprint_matrix(recv_matrix);
    //      printf("\n");
    //  }

    // MPI_Barrier(MPI_COMM_WORLD);
    // if (my_rank == 0) {
    //     printf("\n\n\n ============= send_matrix ============ :\n");
    //     pprint_matrix(send_matrix);
    // }

done:
    MPI_Finalize();

    // matrix_destroy(vector);
    // matrix_destroy(recv_matrix);
    if (my_rank == 0) {
        matrix_destroy(send_matrix);
    }

    return 0;
}