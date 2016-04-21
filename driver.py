import subprocess
import sys

if __name__ == '__main__':
	procnum = 4
	if len(sys.argv) > 1:
		procnum = int(sys.argv[1])
	

	# int N = atoi(argv[1]);
	# double T = atof(argv[2]);
	# double Length = atof(argv[3]);
	# double T1 = atof(argv[4]);
	# double T2 = atof(argv[5]);
	# double TS = atof(argv[6]);
	# double C = atof(argv[7]);
	# double out = atoi(argv[8]);
	# int test_diffic = atoi(argv[9]);
	# int iters_test = atoi(argv[10]);
	# double regularization = atof(argv[11]);

	N = 100
	Time = 1
	Length = 1
	T1 = 1 # Temp on border 1v
	T2 = 2 # Temp on border 2
	TS = 0 # Starting temperature
	C = 1 # Constant to set  dt = 0.3 * h * h / C, where h = L/(N-1)
	out = 1 # Don't need to print results
	test_difficulty = 1000000 # Approximately 0.01 sec for this test.
	iters_between_tests = 0
	regularization = 0.1
	print("N = {}, Time = {}".format(N, Time))

	cmd = ["make"]
	subprocess.call(cmd) 

	cmd = ["mpirun", "-np", str(procnum), "./heat_mpi", str(N), str(Time), str(Length), \
			str(T1), str(T2), str(TS), str(C), str(out), str(test_difficulty),\
			 str(iters_between_tests), str(regularization)]

	print(' '.join(map(str,cmd)))

	print(cmd)
	print(subprocess.check_output(cmd).decode("utf-8")) # Change to string