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

inline void plan_even_distribution(matrix_t *send_matrix, matrix_t *recv_matrix, int sendcounts[], int displacements[], int comm_size) {

    int M = send_matrix -> numRows, N = send_matrix -> numCols;
    int snd_cnt = (M / comm_size) * N;

    for (int rnk=0; rnk<comm_size; rnk++) {
        sendcounts[rnk] = (M/comm_size) * N;
        displacements[rnk] = rnk*(M/comm_size)*N;
    }
}

/**
 * @brief Helper for scatter functions. Calculates how much to send to each process in a communicator of comm_size,
 *        and stores the information in sendcounts and displacements. See @see MPI_Scatterv
 * 
 * @param send_matrix 
 * @param sendcounts 
 * @param displacements 
 * @param comm_size 
 */
inline void plan_distribution(matrix_t *send_matrix, matrix_t *recv_matrix, int sendcounts[], int displacements[], int comm_size) {
    int M = send_matrix -> numRows, N = send_matrix -> numCols;
    int snd_cnt = (M / comm_size) * N;
    int rest_snd_cnt = (M % comm_size) * N;

    for (int rnk=0; rnk<comm_size; rnk++) {
        sendcounts[rnk] = (M/comm_size) * N;
        displacements[rnk] = rnk*(M/comm_size)*N;
    }
    sendcounts[comm_size - 1] += rest_snd_cnt;
    printf("sendcounts and displacements:\n");
    for (int i = 0; i < comm_size; i++) {
        printf("sendcounts[%d] = %2d    ..... displacements[%d] = %3d\n", i, sendcounts[i], i, displacements[i]);
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
                 MPI_COMM_WORLD);
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



void simple_check() {
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

int main(int argc, char *argv[])
{
    int world_size, my_rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    int M = 21, N = 12;

    matrix_t *send_matrix = matrix_factory(M, N);
    matrix_t *vec = matrix_factory(N, 1);
    matrix_t *res = matrix_factory(M, 1);

    populate_matrix_random(send_matrix, 100);
    populate_vector_random(vec, 20);

    int *sendcounts = calloc(world_size, sizeof(int));
    int *displacements = calloc(world_size, sizeof(int));

    matrix_t *recv_matrix;
    if (my_rank == world_size - 1) {
        //printf("Rank[%d]: making extra large receive matrix of size %d x %d\n", my_rank, M / world_size + (M % world_size), N);
        recv_matrix = matrix_factory(M / world_size + (M % world_size), N);
    } else {
        //printf("Rank[%d]: making receive matrix of size %d x %d\n", my_rank,  M / world_size , N);
        recv_matrix = matrix_factory(M / world_size, N);
    }
    plan_even_distribution(send_matrix, recv_matrix, sendcounts, displacements, world_size);
    sendcounts[world_size - 1] += M % world_size;
    scatterv_matrix_mpi_world(send_matrix, sendcounts, displacements, recv_matrix, my_rank, 0, world_size);

    MPI_Barrier(MPI_COMM_WORLD);
    if (my_rank==0) {
        printf("Rank[%d]: M / world_size=%d\n", my_rank, world_size);
        pprint_matrix(recv_matrix);
        printf("\n\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (my_rank==1) {
        printf("Rank[%d]: M / world_size=%d\n", my_rank, world_size);
        pprint_matrix(recv_matrix);
        printf("\n\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (my_rank==2) {
        printf("Rank[%d]: M / world_size=%d\n", my_rank, world_size);
        pprint_matrix(recv_matrix);
        printf("\n\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (my_rank==3) {
        printf("Rank[%d]: M / world_size=%d\n", my_rank, world_size);
        pprint_matrix(recv_matrix);
        printf("\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (my_rank==0) {
        printf("\n\n\n ============= send_matrix ============ :\n");
        pprint_matrix(send_matrix);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // char lc[12];
    // sprintf(lc, "rank[%d][", my_rank);
    // printf("\n\n");
    // printf("\n\n");

    done:
    MPI_Finalize();

    return 0;
}