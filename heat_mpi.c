#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <assert.h>
#include <string.h>
#include "heat_mpi.h"
	
#define DETAILED_PRINT 2


int main(int argc, char **argv)
{
	if (argc != 12) {
		printf("Bad number of arguments: %d", argc);
		return 0;
	}

	int N = atoi(argv[1]);
	double T = atof(argv[2]);
	double Length = atof(argv[3]);
	double T1 = atof(argv[4]);
	double T2 = atof(argv[5]);
	double TS = atof(argv[6]);
	double C = atof(argv[7]);
	int out = atoi(argv[8]);
	int test_diffic = atoi(argv[9]);
	int iters_test = atoi(argv[10]);
	double regularization = atof(argv[11]);

	if (!check_args(N, T, Length, T1, T2, TS, C, iters_test, regularization)) {
		printf("Bad arguments\n");
		return 0;
	}

	double h =  Length / N;  
	// double dt = 0.3 * h * h / C; 
	// int steps = T / dt;

	// To see if program scales. In current difference scheme it will not work, because t ~ O(h^2)
	int steps = 10000;

	assert(steps > 0);
	
	int rank, size;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (out)
		printf("MPI_Init: Rank = %d, Size = %d\n", rank, size);

	// Starting and finishing indexes of the grid. Can be from 1 to N. 
	// Algorithms evaluates cells from starting to finishing inclusively.
	int st = 0;
	int fn = 0; 

	get_borders_initial(rank, size, &st, &fn, N, out);
	
	double *u[2] = {};
	assert(alloc_mem(u, fn - st + 3));
	init_temperature(u[0], st, fn, N, TS, T1, T2);

	double time_start = MPI_Wtime();
	evaluate(rank, size, u, &st, &fn, steps, iters_test, test_diffic, regularization, out);	

	if (!rank) {
		printf("%f", MPI_Wtime() - time_start);
		printf("\n");
	}

	double *res = 0;

	if (!rank) 
		res = collect_results(size, st, fn, N, u[0], out);
	else 
		send_results(rank, st, fn, u[0], out);
	dealloc_mem(u);

	// if (!rank)
	// 	print_results(res, N, out);
	free(res);

	MPI_Finalize();
	return 0;
}

int check_args(int N, double T, double Length, double T1, double T2, double TS, double C, int iters_test, double regularization)
{
	if (N > 0 && T > 0 && Length > 0 && T1 >= 0 && T2 >= 0 && TS >= 0 && C > 0 &&\
	 iters_test >= 0 && regularization <= 1 && regularization >= 0)
		return 1;
	return 0;
}

void evaluate(int rank, int size, double *u[2], int *st, int *fn, int steps, int iters_test, int test_diffic, double regularization, int out)
/* Main function that perfoms all computation. Assumes that u[0] contains valid temperature data. 
Makes number of steps = steps. After each number of iter_test steps performs recalculation of cells number
for current processor. Uses по-очереди? on the right on the left diffusion to do this. 
test_diffic computes the nodes rate. regularization - percentage of rate difference to be tolerated.*/
{
	int s = 0; // denotes the index of array with actual data.
	int exchange_kind = 1; // 1 - means that even sends to process with bigger number. 0 - vice verse. 
	
	if (out)
		printf("evaluate: rank = %d, total steps = %d\n", rank, steps);
	
	for (int i = 0; i < steps; i++) {
		if (out == DETAILED_PRINT)
			printf("evaluate: rank = %d, step = %d, borders: st = %d, fn = %d.\n", rank, i, *st, *fn);
		if (iters_test && !(i % iters_test)) {
			double rate = get_rate(test_diffic);

			//Here is intentional reducing of rate to test if algorithm holds this situation.
			if (rank == 0)
				rate /= 8;	

			resize_tasks(exchange_kind, s, rank, size, u, st, fn, rate, regularization, out);
			if (out)
				printf("evaluate: rank = %d, step = %d, recalib = %d, borders: st_new = %d, fn_new = %d.\n",\
				rank, i, i / iters_test, *st, *fn);
			exchange_kind = !exchange_kind; // exchange direction reversed.
		}	
		// print_results(u[s], *fn - *st + 3);
		make_step(u, s, *fn - *st + 3);

		// Here is intentional reducting the performance
		if (rank == 0) {
			for (int i = 0; i < 7; i++)
				make_step(u, s, *fn - *st + 3);
		}

		borders_temperature_exchange(rank, size, u, s, *fn - *st + 3, out);
		s = !s;
	}
}

void make_step(double *u[2], int s, int mem)
/* Makes one step of evaluation.
s - source array index, !s - destination array index. 
mem - total size of current array. */
{
	for (int i = 1; i < mem - 1; i++)
		u[!s][i] =  u[s][i] + 0.3 * (u[s][i-1] - 2.0 * u[s][i] + u[s][i+1]);
}

int get_borders_initial(int rank, int size, int *st, int *fn, int N, int out)
/* Assigns borders to processes. Uniform distribution 
(each process has approximatelly equal number of cells to evaluate). */
{
	int perproc = N / size;
	*st = 1 + rank * perproc;
	*fn = (rank + 1) * perproc;
	if (rank == size - 1)
		*fn = N;

	if (out)
		printf("get_borders_initial: rank = %d, borders: st = %d, fn = %d\n",\
		rank, *st, *fn);
	return 0;
}

void borders_temperature_exchange(int rank, int size, double *u[2], int s, int mem, int out)
/* Sends st and fn elements and receives st-1, fn+1 (if not nodes on edge).
Even rank: 1) snd fn to bigger 2) rcv fn+1 from bigger 3) rcv  st-1 from smaller 4) snd st to smaller.
Odd rank: 1) rcv st-1 from smaller 2) snd st to smaller 3) snd fn to bigger 4) rcv fn + 1 from smaller. */
{
	double *st_snd = &u[!s][1]; // pointer to the element to pass to smaller rank.
	double *fn_snd = &u[!s][mem - 2]; // pointer to the element to pass to bigger rank.

	double *st_rec = &u[!s][0]; // pointer to the element to rec from smaller rank.
	double *fn_rec = &u[!s][mem - 1]; // pointer to the element to rec from bigger rank.

	if (!(rank % 2)) {
		if (rank != size - 1) { // send to bigger rank. rec from bigger rank.
			MPI_Send(fn_snd, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
			MPI_Recv(fn_rec, 1, MPI_DOUBLE, rank + 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		} else 
			u[!s][mem - 1] = u[s][mem - 1];

		if (rank) { // rec from smaller. send to smaller.
			MPI_Recv(st_rec, 1, MPI_DOUBLE, rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Send(st_snd, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
		} else // Node on left border.
			u[!s][0] = u[s][0];
	} else {
		// rec from smaller. send to smaller.
		MPI_Recv(st_rec, 1, MPI_DOUBLE, rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Send(st_snd, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);

		if (rank != size - 1) { // send to bigger. rec from bigger.
			MPI_Send(fn_snd, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
			MPI_Recv(fn_rec, 1, MPI_DOUBLE, rank + 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		} else // Node on right border. 
			u[!s][mem - 1] = u[s][mem - 1];
	}
}


int max(int a, int b)
{
	if (a > b)
		return a;
	return b;
}

int	init_temperature(double *u, int st, int fn, int N, double TS, double T1, double T2)
/* Initiates temperature in u array with TS - starting temperature. 
If it is a border node, makes an appropriate temperature on border: T1 - left, T2 - right.*/
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


int alloc_mem(double *u[2], int req_mem) 
/* Allocates 2 arrays. Returns 0 if allocation fails, 1 if success.*/ 
{
	u[0] = (double *)calloc(req_mem, sizeof(double));
	if (!u[0])
		return 0;
	u[1] = (double *)calloc(req_mem, sizeof(double));
	if (!u[1]) {
		free(u[0]);
		return 0;
	}
	return 1;
}

void dealloc_mem(double *u[2])
{
	free(u[0]);
	free(u[1]);
}


int resize_tasks(int exchange_kind, int s, int rank, int size, double *u[2], int *st, int *fn, double rate, double regularization, int out)
/* Exchange rates and makes number of cells resizing. */ 
{
	int st_new = 0;
	int fn_new = 0;

	if (!(rank%2))
		get_borders_even(exchange_kind, rank, size, *st, *fn, &st_new, &fn_new, rate);
	else 
		get_borders_odd(exchange_kind, rank, size, *st, *fn, &st_new, &fn_new, rate);

	if (out)
		printf("resize_tasks: rank = %d, rate = %f, exchange_kind = %d, st = %d, fn = %d, st_new = %d, fn_new = %d\n",\
			rank, rate, exchange_kind, *st, *fn, st_new, fn_new);


	int should = should_reallocate(exchange_kind, rank, size, *st, *fn, st_new, fn_new, regularization, out);

	if (!should)
		return 0;

	reallocate(rank, s, u, *st, *fn, st_new, fn_new);
	*st = st_new;
	*fn = fn_new;

	return 1;
}

int get_borders_even(int exchange_kind, int rank, int size, int st, int fn, int *st_new, int *fn_new, double rate)
{
	double rate_neigh = 0;
	if (exchange_kind && rank != size - 1) {
		MPI_Recv(&rate_neigh, 1, MPI_DOUBLE, rank + 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&fn, 1, MPI_INT, rank + 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		
		int middle = (int) (rate * (double)(fn - st) / (double)(rate + rate_neigh));
		
		*st_new = st;
		*fn_new = st + middle;

		int st_neigh = st + middle + 1;
		MPI_Send(&st_neigh, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);	
	} else if (!exchange_kind && rank) { 
		MPI_Recv(&rate_neigh, 1, MPI_DOUBLE, rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&st, 1, MPI_INT, rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		int middle = (int) (rate * (double)(fn - st) / (double)(rate + rate_neigh));
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


int get_borders_odd(int exchange_kind, int rank, int size, int st, int fn, int *st_new, int *fn_new, double rate)
{
	if (exchange_kind) {
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

int should_reallocate(int exchange_kind, int rank, int size, int st, int fn, int st_new, int fn_new, double regularization, int out) 
/* Checks own reallocate status and exchanges it with partner. If one of nodes wants to reallocate - reallocation happens.*/
{
	int should = 0;
	if (fn - st <= 0) {
		if (fn_new - st_new <= 0)
			should = 0;
		else 
			should = 1;
	}
	if (1 - (fn_new - st_new) / (fn - st) > regularization)
		should = 1;

	if (out)
		printf("resize_tasks: rank = %d, should (before exchange) = %d\n", rank, should);
	should = reallocate_status_exchange(exchange_kind, rank, size, should);
	if (out)
		printf("resize_tasks: rank = %d, should (after exchange) = %d\n", rank, should);
	
	return should; 
}


int reallocate_status_exchange(int exchange_kind, int rank, int size, int should)
{
	int should_neigh = 0;

	if (!(rank%2)) {
		if (exchange_kind && rank != size - 1) {
			MPI_Recv(&should_neigh, 1, MPI_INT, rank + 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Send(&should, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
		} else if (!exchange_kind && rank) {
			MPI_Recv(&should_neigh, 1, MPI_INT, rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Send(&should, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD);
		}
	} else {
		if (exchange_kind) {
			MPI_Send(&should, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD);
			MPI_Recv(&should_neigh, 1, MPI_INT, rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		} else if (!exchange_kind && rank != size - 1) { 
			MPI_Send(&should, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
			MPI_Recv(&should_neigh, 1, MPI_INT, rank + 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	}
	return max(should, should_neigh);
}

int reallocate(int rank, int s, double *u[2], int st, int fn, int st_new, int fn_new) 
/* Reallocation procedure. */
{
	double *u_new[2] = {}; // To initiate new mem.

	u_new[0] = (double *)calloc(fn_new - st_new + 3, sizeof(double));
	assert(u_new[0]); // Can handle this case. Just continue with the same arrays.
	u_new[1] = (double *)calloc(fn_new - st_new + 3, sizeof(double));
	assert(u_new[1]);

	double *u_src = u[s]; // Array with actual data.
	double *u_dst = u_new[s]; // Array, where we want actual data after reallocation.
	if (fn - st < fn_new - st_new) // Node needs to receive data.
		reallocate_rcv(rank, u_dst, u_src, st, fn, st_new, fn_new);
	else // Node needs to send data.
		reallocate_snd(rank, u_dst, u_src, st, fn, st_new, fn_new);

	free(u[0]);
	free(u[1]);

	u[0] = u_new[0];
	u[1] = u_new[1];

	return 0;
}


int reallocate_rcv(int rank, double *u_dst, double *u_src, int st, int fn, int st_new, int fn_new)
/* This function is called if node needs to receive data after resizing. */
{

	int to_rec = (fn_new - st_new) - (fn - st); // How many doubles should receive
	double *dst_rcmem = u_dst; // address, to where doubles will be received.
	double *dst_cpmem = u_dst; // address, to where doubles will be copied.
	
	int rank_neigh = rank - 1; // From whom I wan't to receive. 
	if (fn != fn_new) { // st == st_new; fn < fn_new;
		dst_rcmem += fn - st + 3;
		rank_neigh += 2;
	} else // st_new < st; fn == fn_new;
		dst_cpmem += st - st_new;


	MPI_Recv(dst_rcmem, to_rec, MPI_DOUBLE, rank_neigh, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	memcpy(dst_cpmem, u_src, sizeof(double) * (fn - st + 3));

	return 0;
}

int reallocate_snd(int rank, double *u_dst, double *u_src, int st, int fn, int st_new, int fn_new)
/* This function is called if node needs to send data after resizing. */
{
	int to_send = (fn - st) - (fn_new - st_new); // How many doubles should send.
	double *src_sdmem = u_src + 2; // address, from which doubles will be sent to neigh.
	double *src_cpmem = u_src; // address, from which doubles will be copied to u_dst.

	int rank_neigh = rank - 1;
	if (fn != fn_new) { // Here fn_new < fn; st == st_new
		src_sdmem += fn_new - st_new - 1;
		rank_neigh += 2;
	} else  // Here fn == fn_new; st < st_new
		src_cpmem += st_new - st;
	
	MPI_Send(src_sdmem, to_send, MPI_DOUBLE, rank_neigh, 0, MPI_COMM_WORLD);
	// fn_new - st_new + 1 -> number of bytes that come from old array.
	memcpy(u_dst, src_cpmem, sizeof(double) * (fn_new - st_new + 3)); 

	return 0;
}

double get_rate(int param)
/* Function that returns rate ("speed") of current node. */
{
	double time_start = MPI_Wtime();
	test_task(param);
	double time_duration = MPI_Wtime() - time_start;
	return 1/time_duration;
}

void test_task(int param)
/* Testing task that is performed to measure speed. */
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

void send_results(int rank, int st, int fn, double *u, int out)
/* Sending results to zero node. */
{
	int data_size = max(fn - st + 1, 0);

	if (out == DETAILED_PRINT)
		printf("send_results: rank = %d, data_size = %d (st = %d, fn = %d) (sending)\n",\
		 rank, data_size, st, fn);
	
	MPI_Send(&data_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	MPI_Send(u + 1, data_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	if (out)
		printf("send_results: rank = %d, data_size = %d\n (sent)", rank, data_size);
}

double* collect_results(int size, int st, int fn, int N, double *u, int out) 
/* When zero node calls this function it collects results from all other nodes. */
{
	double *res = (double *)calloc(N, sizeof(double));
	assert(res);

	double *pointer = res;

	int data_size = max(fn - st + 1, 0);
	memcpy(pointer, u + 1, data_size * sizeof(double));
	pointer += data_size;
	
	if (out)
		printf("collect_results: rank = 0, data_size = %d (copied)\n", data_size);	
	for (int i = 1; i < size; i++) {
		if (out == DETAILED_PRINT)
			printf("collect_results: rank = %d (receiving)\n", i);
			
		MPI_Recv(&data_size, 1, MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(pointer, data_size, MPI_DOUBLE, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		pointer += data_size;

		if (out)
			printf("collect_results: rank = %d, data_size = %d (received)\n", i, data_size);
	}
	return res;
}

void print_results(double *res, int N, int out)
{
	if (out)
		printf("print_res:\n");
	for (int i = 0; i < N; i++)
		printf("%f ", res[i]);
	printf("\n");
}
