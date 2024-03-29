#!./roundup
# Checking of heat_no_mpi and heat_simple_mpi. Doesn't cover main program. 
make
heat='heat_no_mpi'
heat_mpi='heat_simple_mpi'

it_check_N10P1() {
	./$heat 10 0.1 > s.txt
	mpirun -np 1 ./$heat_mpi 10 0.1 > p.txt
	diff -q s.txt p.txt
}

it_check_N10P2() {
	./$heat 10 0.1 > s.txt
	mpirun -np 2 ./$heat_mpi 10 0.1 > p.txt
	diff -q s.txt p.txt
}

it_check_N10P3() {
	./$heat 10 0.1 > s.txt
	mpirun -np 3 ./$heat_mpi 10 0.1 > p.txt
	diff -q s.txt p.txt
}

it_check_N10P7() {
	./$heat 10 0.1 > s.txt
	mpirun -np 7 ./$heat_mpi 10 0.1 > p.txt
	diff -q s.txt p.txt
}

it_check_N10P9() {
	./$heat 10 0.1 > s.txt
	mpirun -np 9 ./$heat_mpi 10 0.1 > p.txt
	diff -q s.txt p.txt
}

it_check_N100P2() {
	./$heat 100 0.1 > s.txt
	mpirun -np 2 ./$heat_mpi 100 0.1 > p.txt
	diff -q s.txt p.txt
}

it_check_N100P7() {
	./$heat 100 0.1 > s.txt
	mpirun -np 7 ./$heat_mpi 100 0.1 > p.txt
	diff -q s.txt p.txt
}
