void evaluate(int rank, int size, double *u[2], int *st, int *fn, int steps,\
 int iters_test, int test_diffic, double regularization, int out);

int get_borders_even(int test_kind, int rank, int size, int st, int fn,\
 	int *st_new, int *fn_new, double rate);
int get_borders_odd(int test_kind, int rank, int size, int st, int fn,\
 	int *st_new, int *fn_new, double rate);

void make_step(double *u[2], int s, int mem);
void borders_temperature_exchange(int rank, int size, double *u[2], int s, int mem, int out);
// void border_filling_neigh(int rank, int rank_neigh, double *u[2], int s, int index);

int should_reallocate(int test_kind, int rank, int size, int st, int fn, int st_new, int fn_new, double regularization, int out);
int resize_tasks(int test_kind, int s, int rank, int size,\
 	double *u[2], int *st, int *fn, double rate, double regularization, int out);
int reallocate_status_exchange(int test_kind, int rank, int size, int should);
int reallocate(int rank, int s, double *u[2], int st, int fn, int st_new, int fn_new);
int reallocate_rcv(int rank, double *u_dst, double *u_src, int st, int fn, int st_new, int fn_new);
int reallocate_snd(int rank, double *u_dst, double *u_src, int st, int fn, int st_new, int fn_new);

int get_borders_initial(int rank, int size, int *st, int *fn, int N, int out);
// int get_borders_initial(int rank, int size, int *st, int *fn, double rate, int N, int out);
// void get_borders_initial_divide_proportion(int size, double *rates,
//  	int *sts, int *fns, int N, int out);

void test_task(int param);
double get_rate(int param);
int max(int a, int b);

int	init_temperature(double *u, int st, int fn, int N, double TS, double T1, double T2);
int alloc_mem(double *u[2], int req_mem);
void dealloc_mem(double *u[2]);


double* collect_results(int size, int st, int fn, int N, double *u, int out);
void send_results(int rank, int st, int fn, double *u, int out);
void print_results(double *res, int N, int out);

int check_args(int N, double T, double Length, double T1, double T2, double TS, double C, int iters_test, double regularization);

