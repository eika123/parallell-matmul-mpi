CC=mpicc
CFLAGS=-O2


BINDIR=bin

.PHONY: clean

all: matmul_mpi

matmul_mpi: matmul_mpi.c
	mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) -o $(BINDIR)/$@ $^ -lm

clean:
	rm -f $(BINDIR)/*