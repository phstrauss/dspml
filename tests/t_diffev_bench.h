/*  Differential evolution solver heurisitc
    Lester E. Godwin
    Partially rewrote, enhanced (EitherOr and Probabilist Parent Centric evolution strategies)
    plus maintained by Philippe Strauss, Spring 2012.  */


/* see p. 55 Table 2.1 for others tests ideas */


#ifndef _OPTIM_BENCH_H
#define _OPTIM_BENCH_H

#ifdef __cplusplus
extern "C" {
#endif


#define MAX_DIM 200
#define PCX_STRATEGY Rand1Exp
#define PCX_PROB 0.07
#define PCX_SIGMA2_A 0.7
#define PCX_SIGMA2_B 0.3
#define STD_F 0.7
#define STD_XOVER 0.6
#define STD_STRATEGY Best1Exp
#define DIM30_POP 350
#define DIM100_POP 350
#define DIM30_GEN 100000
#define DIM100_GEN 100000
#define F_RND_FACT 0.001
#define F_POP_FACT 0.0


double Rosenbrock2D(int n, double bestEnergy, double *trial, bool *bAtSolution);
double Ackley(int n, double bestEnergy, double *trial, bool *bAtSolution);
double Griewangk(int n, double bestEnergy, double *trial, bool *bAtSolution);
double Rastrigin(int n, double bestEnergy, double *trial, bool *bAtSolution);
double Schwefel(int n, double bestEnergy, double *trial, bool *bAtSolution);
double Michalewicz(int n, double bestEnergy, double *trial, bool *bAtSolution);


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif


