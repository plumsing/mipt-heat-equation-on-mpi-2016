CFLAGS=-Wall -g -std=c99

all: heat_mpi heat_no_mpi heat_simple_mpi 

heat_mpi: 
	mpicc $(CFLAGS) -o heat_mpi heat_mpi.c

heat_simple_mpi:
	mpicc $(CFLAGS) -o heat_simple_mpi heat_simple_mpi.c

heat_no_mpi:
	gcc $(CFLAGS) -o heat_no_mpi heat_no_mpi.c

clean:
	rm heat_mpi heat_simple_mpi heat_no_mpi	*.txt