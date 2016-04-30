# MPI_heat_equation

Program that allows to compute heat equation using MPI.

Main program is heat_mpi.c. It uses explicit finite difference scheme for evaluation. 
The main feature of the program - diffusion rebalancing of loads on node based on the nodes performance. 

Diffusion rebalancing is made by comparing performance of two neighboring nodes and redistribution of the
tasks according to there speed. This technique allows to continue effective evaluation even if some of the 
nodes has problems with performance. 

To launch program you can use run_script.sh or test_script.sh.

Programs require gets parameters:

1. N - number of cells to use in evaluation.
2. T - time of evolution.
3. Length - length of the string.
4. T1 - temperature on the left border.
5. T2 - temperature on the right border.
6. TS - starting temperature of the string.
7. C - constant in formula that define dt (time step) to h (spacial step) ratio: dt = 0.3 * h * h / C.
8. out - variable that stands for logging.
9. test_diffic - difficulty of the test that is used to determine nodes rate.
10. iters_test - number of time steps between task resizing.
11. regularization - percentage of task difference to be tolerated.
