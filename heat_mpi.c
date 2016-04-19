#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <assert.h>

int reallocate(int rank, int to_fill, double *u[2], int st, int fn, int st_new, int fn_new);
int reallocate_rcv(int rank, double *u_dst, double *u_src, int st, int fn, int st_new, int fn_new);
int reallocate_snd(int rank, double *u_dst, double *u_src, int st, int fn, int st_new, int fn_new);
int get_borders_initial(int rank, int size, int *st, int *fn, double rate, int N, int out);
void get_borders_initial_divide_proportion(int size, double *rates, int *sts, int *fns, int N, int out);


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

double get_rate(int param)
{
	double time_start = MPI_Wtime();
	test_task(param);
	double time_duration = MPI_Wtime() - time_start;
	return 1/time_duration;
}

int max(int a, int b)
{
	if (a > b)
		return a;
	return b;
}

int manage_memory(double *u[2], int mem, int req) 
{	
	if (mem >= req) {
		if (mem > req * 2) {
			u[0] = realloc(u[0], (int)mem / 2);
			u[1] = realloc(u[1], (int)mem / 2);
		return mem / 2;
	} else {
		u[0] = realloc(u[0], max(req, mem * 2));
		u[1] = realloc(u[1], max(req, mem * 2);
		return max(req, mem * 2);
	}
	return 0; // Shouldn't get here.
}



int	init_memory(double *u, st, fn, N, TS, T1, T2)
{
	int mem = fn - st + 1 + 2;
	for (int i = 0; i < mem; i++)
		u[i] = TS;
	if (st == 1)
		u[0] = T1;
	if (fn == N)
		u[mem - 1] = T2;
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

	// Number of iterations between tests. If zero then no testing is made.
	int iters_test = atoi(argv[10]);

	if (N <= 0 || T <= 0 || Length <= 0 || T1 < 0 || T2 < 0 || TS < 0 || C <= 0 || iters_test < 0) {
		printf("Bad arguments\n");
		return 0;
	}


	double *u[2] = {};
	double h =  Length / N; // I don't no why division by N - 1 in previous program 
	double dt = 0.3 * h * h / C;
	int steps = T / dt;
	assert(steps > 0);

	int un = 1;
	int mem = 0;
	int rank, size;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	double rate = get_rate(test_difficulty);

	int st = 0; 
	int fn = 0;	

	get_borders(rank, size, st, fn, rate, N, out);
	if (out)
		printf("Step = %d. Recalibration = %d. Rank = %d. Rate = %f. Borders: st = %d, fn = %d.\n",\
		i, i / iters_test, rank, rate, st_new, fn_new);
	
	
	int req_mem = fn - st + 1 + 2;
	int mem = manage_memory(u, 0, req_mem); // This variable will held the amount of memory process holds.
	int test_kind = 1; // 1 - means that even sends to process with bigger number. 0 - vice verse. 

	init_memory(u[0], mem, st, fn, N, TS, T1, T2);
	int to_fill = 1; // 1 - means, that on the next step we fill array with number 1. 0 -//- with 0.
	for (int i = 0; i < steps; i++) {
		if (iters_test && i && !(i % iters_test)) {
			resize_tasks(test_kind, to_fill, rank, size, &st, &fn, test_difficulty, out);
			int get_borders_initial(int rank, int size, int *st, int *fn, double rate, int N, int out) {

			test_kind = !test_kind;
		}	
		make_step(u, mem, rank, size)

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


int get_new_borders_even(int test_kind, int rank, int size, int st, int fn, int *st_new, int *fn_new, int test_difficulty)
{
	double rate = get_rate(test_difficulty);
	double rate_neigh = 0;
	if (test_kind && rank != size - 1) {
		MPI_Recv(&rate_neigh, 1, MPI_DOUBLE, rank + 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&fn, 1, MPI_INT, rank + 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		
		int middle = (int) rate * (fn - st) / (rate + rate_neigh);
		
		*st_new = st;
		*fn_new = st + middle;

		int st_neigh = st + middle + 1;
		MPI_Send(&st_neigh, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);	
	} else if (!test_kind && rank) { 
		MPI_Recv(&rate_neigh, 1, MPI_DOUBLE, rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&st, 1, MPI_INT, rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		int middle = (int) rate * (fn - st) / (rate + rate_neigh);
		*st_new = st + middle + 1;
		*fn_new = fn;
		int fn_neigh = st + middle;
		MPI_Send(&fn_neigh, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD);
	} else {
		*st_new = st;
		*fn_new = fn;
	}
	return 0;
}


int get_new_borders_odd(int test_kind, int rank, int size, int st, int fn, int *st_new, int *fn_new, int test_difficulty)
{
	double rate = get_rate(test_difficulty);
	if (test_kind) {
		MPI_Send(&rate, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
		MPI_Send(&fn, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD);

		*fn_new = fn;
		MPI_Recv(st_new, 1, MPI_INT, rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	} else if (rank != size - 1) {
		MPI_Send(&rate, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
		MPI_Send(&st, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);

		*st_new = st;
		MPI_Recv(fn_new, 1, MPI_INT, rank + 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	} else {
		*st_new = st;
		*fn_new = fn;
	}
	return 0;
}

int should_reallocate(int st, int fn, int st_new, int fn_new, double regularization) {
	if (fn - st <= 0)
		if (fn_new - st_new <= 0)
			return 0;
		else 
			return 1;
	if (1 - (fn_new - st_new) / (fn - st) > regularization)
		return 1;
	return 0; 
}


int reallocate_status_exchange(int test_kind, int rank, int size, int should)
{
	int should_neigh = 0;
	if (!(rank%2)) {
		if (test_kind && rank != size - 1) {
			MPI_Recv(&should_neigh, 1, MPI_INT, rank + 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Send(&should, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
		} else if (rank) {
			MPI_Recv(&should_neigh, 1, MPI_INT, rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Send(&should, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD);
		}
	} else {
		if (test_kind) {
			MPI_Send(&should, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD);
			MPI_Recv(&should_neigh, 1, MPI_INT, rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		} else if (rank != size - 1) { 
			MPI_Send(&should, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
			MPI_Recv(&should_neigh, 1, MPI_INT, rank + 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	}
	return max(should, should_neigh);
}

int resize_tasks(int test_kind, int to_fill, int rank, int size, int *st, int *fn, int test_difficulty, double regularization, int out) 
{
	int st_new = 0;
	int fn_new = 0;

	if (!(rank%2))
		get_new_borders_even(test_kind, rank, size, *st, *fn, &st_new, &fn_new, test_difficulty);
	else 
		get_new_borders_odd(test_kind, rank, size, *st, *fn, &st_new, &fn_new, test_difficulty);

	if (out)
		printf("resize_tasks: rank = %d, size = %d, test_kind = %d, st = %d, fn = %d, st_new = %d, fn_new = %d\n",\
			rank, size, test_kind, *st, *fn, st_new, fn_new);


	int should = should_reallocate(*st, *fn, st_new, fn_new, regularization);
	if (out)
		printf("resize_tasks: rank = %d, size = %d, test_kind = %d, should (before exchange) = %d\n", rank, size, test_kind, should);
	should = reallocate_status_exchange(test_kind, rank, size, should);
	if (out)
		printf("resize_tasks: rank = %d, size = %d, test_kind = %d, should (after exchange) = %d\n", rank, size, test_kind, should);
	
	if (!should)	
		return 0;


	// Exchange should happen here.

	return 1;

}

int reallocate(int rank, int to_fill, double *u[2], int st, int fn, int st_new, int fn_new) 
{
	double *u_new[2] = {}; // To initiate new mem.

	u_new[0] = (double *)calloc(fn_new - st_new + 3, sizeof(double));
	assert(u_new[0]); // Can handle this case. Just continue with the same arrays.
	u_new[1] = (double *)calloc(fn_new - st_new + 3, sizeof(double));
	assert(u_new[1]);

	u_src = u[!to_fill]; // Array with actual data.
	u_dst = u_new[!to_fill]; // Array, where we want actual data after reallocation.
	if (fn - st < fn_new - st_new) 
		reallocate_rcv(rank, u_dst, u_src, st, fn, st_new, fn_new);
	else
		reallocate_snd(rank, u_dst, u_src, st, fn, st_new, fn_new);

	free(u[0]);
	free(u[1]);

	u[0] = u_new[0];
	u[1] = u_new[1];

	return 0;
}

int reallocate_rcv(int rank, double *u_dst, double *u_src, int st, int fn, int st_new, int fn_new)
{
	int to_rec = (fn_new - st_new) - (fn - st); // How many doubles should receive
	double *dst_rcmem = u_dst + 1; // address, to where doubles will be received.
	double *dst_cpmem = u_dst + 1; // address, to where doubles will be copied.
	
	int rank_neigh = rank - 1; // From whom I wan't to receive. 
	if (fn != fn_new) { // st == st_new; fn < fn_new;
		dst_rcmem += fn - st + 1;
		rank_neigh += 2;
	} else // st_new < st; fn == fn_new;
		dst_cpmem += st - st_new;


	MPI_Recv(dst_rcmem, to_rec, MPI_DOUBLE, rank_neigh, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	// fn - st + 1 -> number of bytes that come from old array.
	memcpy(dst_cpmem, u_src + 1, sizeof(double) * (fn - st + 1));

	return 0;
}

int reallocate_snd(int rank, double *u_dst, double *u_src, int st, int fn, int st_new, int fn_new)
{
	int to_send = (fn - st) - (fn_new - st_new); // How many doubles should send.
	double *src_sdmem = u_src + 1; // address, from which doubles will be sent to neigh.
	double *src_cpmem = u_src + 1; // address, from which doubles will be copied to u_dst.

	int rank_neigh = rank - 1;
	if (fn != fn_new) { // Here fn > fn_new; st == st_new
		src_sdmem += fn_new - st_new + 1;
		rank_neigh += 2;
	} else  // Here fn == fn_new; st < st_new
		src_cpmem += st_new - st;
	
	MPI_Send(src_sdmem, to_send, MPI_DOUBLE, rank_neigh, 0, MPI_COMM_WORLD);
	// fn_new - st_new + 1 -> number of bytes that come from old array.
	memcpy(u_dst + 1, src_cpmem, sizeof(double) * (fn_new - st_new + 1)); 

	return 0;
}

int get_borders_initial(int rank, int size, int *st, int *fn, double rate, int N, int out)
{
	if (!rank) {
		double *rates = (double *)calloc(size, sizeof(double));
		assert(rates);
		rates[0] = rate;
		int *sts = (int *)calloc(size, sizeof(int));
		int *fns = (int *)calloc(size, sizeof(int));
		int i = 0;
		for (i = 1; i < size; i++)
			MPI_Recv(rates + i, 1, MPI_DOUBLE, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		get_borders_initial_divide_proportion(size, rates, sts, fns, N, out);	

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


void get_borders_initial_divide_proportion(int size, double *rates, int *sts, int *fns, int N, int out) 
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

