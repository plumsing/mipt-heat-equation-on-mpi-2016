CFLAGS = -Wall -g -std=c99

all: 
	mpicc $(CFLAGS) -o heat_mpi heat_mpi.c