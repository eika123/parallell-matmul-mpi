

NUM_PROCESSES=4
BINDIR=bin

mpirun -np $NUM_PROCESSES $BINDIR/matmul_mpi
