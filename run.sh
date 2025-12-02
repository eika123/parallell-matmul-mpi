#/bin/bash

if [[ -n "$1" ]]; then
    NUM_PROCESSES=$1
else
    NUM_PROCESSES=4
fi

BINDIR=bin

mpirun -np $NUM_PROCESSES $BINDIR/matmul_mpi
