CC=mpicc
CFLAGS=-O2


BINDIR=bin


.PHONY: clean

all: matmul_mpi

matmul_mpi: matmul_mpi.c matrix_common.c matrix_util.c
	mkdir -p $(BINDIR)
	$(CC) $^ $(CFLAGS) -o $(BINDIR)/$@  -lm

clean:
	rm -f $(BINDIR)/*