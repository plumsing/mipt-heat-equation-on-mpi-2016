#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <assert.h>

void test_task(int param)
{
	int i = 0;
	int a[2] = {};
	a[0] = 1;
	a[1] = 1;
	for (i = 0; i < param; i++) {
		a[0] = a[0] + a[1];
		a[1] = a[0] * a[1];
	}
}

void divide_proportion(int size, double *rates, int *sts, int *fns, int N, int out) 
{
	double total_rate = 0;
	for (int i = 0; i < size; i++)
		total_rate += rates[i];

	if (out) {
		printf("Rates:\n");
		for (int i = 0; i < size; i++)
			printf("%f ", rates[i]);
		printf("\ntotal_rate = %f\n", total_rate);
	}
	sts[0] = 1;
	for (int i = 0; i < size - 1; i++) {
		fns[i] = sts[i] + (int) (rates[i] * N / total_rate);
		sts[i + 1] = fns[i] + 1;		
	}
	
	fns[size - 1] = N;
}

int get_borders(int rank, int size, int *st, int *fn, double rate, int N, int out) {
	if (!rank) {
		double *rates = (double *)calloc(size, sizeof(double));
		assert(rates);
		rates[0] = rate;
		int *sts = (int *)calloc(size, sizeof(int));
		int *fns = (int *)calloc(size, sizeof(int));
		int i = 0;
		for (i = 1; i < size; i++)
			MPI_Recv(rates + i, 1, MPI_DOUBLE, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		divide_proportion(size, rates, sts, fns, N, out);	

		for (i = 1; i < size; i++) {
			int params[2] = {};
			params[0] = sts[i];
			params[1] = fns[i];
			MPI_Send(params, 2, MPI_INT, i, 0, MPI_COMM_WORLD); // Try MPI_Ibsend 
		}
		*st = sts[0];
		*fn = fns[0];
		free(rates);
		free(sts);
		free(fns);
	} else {
		MPI_Send(&rate, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		int params[2] = {};
		MPI_Recv(params, 2, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		*st = params[0];
		*fn = params[1];
	}

	return 0;
}

double get_rate(int param)
{
	double time_start = MPI_Wtime();
	test_task(param);
	double time_duration = MPI_Wtime() - time_start;
	return 1/time_duration;
}

int manage_memory(double *u[2], int mem, int req, double mem_multrate) 
{
	if (mem >= req) {
		if (mem > req * mem_multrate)
			mem = realloc(u[0], (int)req * mem_multrate);
		return mem;
	} else 
		mem = realloc(u[0], (int)req * mem_multrate);
	return 0;
}

int main(int argc, char **argv)
{
	if (argc < 11) {
		printf("Not enough arguments.\n");
	}
	int N = atoi(argv[1]);
	double T = atof(argv[2]);
	double Length = atof(argv[3]);
	double T1 = atof(argv[4]);
	double T2 = atof(argv[5]);
	double TS = atof(argv[6]);
	double C = atof(argv[7]);
	double out = atoi(argv[8]);
	int test_difficulty = atoi(argv[9]);
	int iters_test = atoi(argv[10]);

	if (N <= 0 || T <= 0 || Length <= 0 || T1 < 0 || T2 < 0 || TS < 0 || C <= 0) {
		printf("Bad arguments\n");
		return 0;
	}


	double *u[2] = {};
	double h =  Length / N; // I don't no why division by N - 1 in previous program 
	double dt = 0.3 * h * h / C;
	int steps = T / dt;
	int un = 1;
	int mem = 0;
	int rank, size;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);



	for (int i = 0; i < steps; i++) {
		if (i % iters_test == 0) {

			double rate = get_rate(test_difficulty);

			int st_new = 0;
			int fn_new = 0;
			get_borders(rank, size, &st_new, &fn_new, rate, N, out);
			if (out)
				printf("Step = %d. Recalibration = %d. Rank = %d. Rate = %f. Borders: st = %d, fn = %d.\n",\
				i, i / iters_test, rank, rate, st_new, fn_new);
			int new_mem = manage_memory(u, mem, fn_new - st_new + 1);
		}
	}
	
	// /* Выделение памяти. */
	// u[0] = (double*)calloc(N, sizeof(double));
	// u[1] = (double*)calloc(N, sizeof(double));

	// /* Граничные условия. */
	// u[0][0] = u[1][0] = T1;
	// u[0][N - 1] = u[1][N - 1] = T2;

	// int perproc = N / size;
	// int st = perproc * rank;
	// int fn = st + perproc;
	// if (rank == size - 1) fn = N - 1;
	// if (rank == 0) st = 1;
	// /* Цикл интегрирования. */
	// for (i = 0; i < steps; i++) {
	// 	/* Начало обмена. */
	// 	if (rank%2) {
	// 		if (rank != size - 1) {
	// 			MPI_Ssend(&u[un][fn-1], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
	// 			MPI_Recv(&u[un][fn], 1, MPI_DOUBLE, rank + 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	// 		}
	// 		MPI_Recv(&u[un][st-1], 1, MPI_DOUBLE, rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	// 		MPI_Ssend(&u[un][st], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
	// 	} else {
	// 		if (rank) {
	// 			MPI_Recv(&u[un][st-1], 1, MPI_DOUBLE, rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	// 	 		MPI_Ssend(&u[un][st], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
	// 		} 
	// 		if (rank != size - 1) {
	// 			MPI_Ssend(&u[un][fn-1], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
	// 			MPI_Recv(&u[un][fn], 1, MPI_DOUBLE, rank + 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	// 		}
	// 	}
			
	// 	/* Конец обмена. */
	// 	for (j = st; j < fn; j++) {
	// 		u[!un][j] = u[un][j] + 0.3 * (u[un][j-1] - 2.0 * u[un][j] + u[un][j+1]);
	// 	}
	// 	un = !un;
	// }
	// if (rank == 0) {
	// 	for (i = 1; i < size; i++) {
	// 		MPI_Recv(u[un] + perproc * i, perproc  + (i == size - 1 ? N % size : 0), MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	// 	}
	// 	/* Вывод на экран. */
	// 	if (out) {
	// 		for (i = 0; i < N; i++) {
	// 			printf("%f %f\n", h * i, u[un][i]);
	// 		}
	// 	}
	// } else {
	// 	MPI_Send(u[un] + perproc * rank, perproc + (rank == size - 1 ? N % size : 0), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	// }
	// free(u[0]);
	// free(u[1]);
	MPI_Finalize();
	return 0;
}
