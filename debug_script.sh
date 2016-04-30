# Script shows how to use gdb with MPI.
make
mpirun -np 2 xterm -e gdb --args ./heat_mpi 100 0.02 1 1 2 0 1 1 1000000 100 0.3 
