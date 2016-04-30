# Script that demonstrate that program efficienty resizes tasks. 
make
for numproc in {1..2}
do
	echo $numproc
	mpirun -np $numproc ./heat_mpi 5000 0.05 1 1 2 0 1 0 1000 0 0.2
	mpirun -np $numproc ./heat_mpi 5000 0.05 1 1 2 0 1 0 1000 100 0.2
done



