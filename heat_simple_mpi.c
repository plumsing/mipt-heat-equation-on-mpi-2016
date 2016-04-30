#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define L 1.0
#define C 1.0
#define T1 1.0
#define T2 2.0

/* Computes heat equation using MPI. No task resizing. */
int main(int argc, char **argv)
{
	if (argc < 3) {
		printf("Usage: %s N T\n", argv[0]);
	}
	int N = atoi(argv[1]);
	double T = atof(argv[2]);
	double *u[2];
	double h = L / (N - 1);
	double dt = 0.3 * h * h / C;
	int steps = T / dt;
	int i, j, un = 1;
	int rank, size;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	/* Выделение памяти. */
	u[0] = (double*)calloc(N, sizeof(double));
	u[1] = (double*)calloc(N, sizeof(double));

	/* Граничные условия. */
	u[0][0] = u[1][0] = T1;
	u[0][N - 1] = u[1][N - 1] = T2;

	int perproc = N / size;
	int st = perproc * rank;
	int fn = st + perproc;
	if (rank == size - 1) fn = N - 1;
	if (rank == 0) st = 1;
	/* Цикл интегрирования. */
	for (i = 0; i < steps; i++) {
		/* Начало обмена. */
		if (rank%2) {
			if (rank != size - 1) {
				MPI_Ssend(&u[un][fn-1], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
				MPI_Recv(&u[un][fn], 1, MPI_DOUBLE, rank + 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
			MPI_Recv(&u[un][st-1], 1, MPI_DOUBLE, rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Ssend(&u[un][st], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
		} else {
			if (rank) {
				MPI_Recv(&u[un][st-1], 1, MPI_DOUBLE, rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		 		MPI_Ssend(&u[un][st], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
			} 
			if (rank != size - 1) {
				MPI_Ssend(&u[un][fn-1], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
				MPI_Recv(&u[un][fn], 1, MPI_DOUBLE, rank + 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
		}
			
		/* Конец обмена. */
		for (j = st; j < fn; j++) {
			u[!un][j] = u[un][j] + 0.3 * (u[un][j-1] - 2.0 * u[un][j] + u[un][j+1]);
		}
		un = !un;
	}
	if (rank == 0) {
		for (i = 1; i < size; i++) {
			MPI_Recv(u[un] + perproc * i, perproc  + (i == size - 1 ? N % size : 0), MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		/* Вывод на экран. */
		for (i = 0; i < N; i++) {
			printf("%f %f\n", h * i, u[un][i]);
		}
	} else {
		MPI_Send(u[un] + perproc * rank, perproc + (rank == size - 1 ? N % size : 0), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}
	free(u[0]);
	free(u[1]);
	MPI_Finalize();
	return 0;
}
