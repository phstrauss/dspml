/*  Differential evolution solver heurisitc
    Lester E. Godwin
    Partially rewrote, plain ANSII C port,
    enhanced (EitherOr and Probabilist Parent Centric evolution strategies)
    plus maintained by Philippe Strauss, Spring 2012.  */


#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include "de.h"
#include "perr.h"
#include "monotime.h"
#include "optim_bench.h"


typedef enum { STD, PCX, EITHER_OR } de_algo_t;


void print_coeffs(de_t *st){
    int i;

    printf("Best Coefficients:\n");
    for (i = 0; i < st->nDim; i++)
        printf("[%d]: %E\n", i, st->bestSolution[i]);
    printf("\n");
}

/* start gently : rosenbrock's 2D */
double t_rosenbrock(de_algo_t algo) {
    double min[MAX_DIM];
    double max[MAX_DIM];
    int i;
    de_t prob_data;
    double begin, end;

    int ndim = 2;
    int npop = 30;
    int ngen = 100;

    for (i=0; i<ndim; i++) {
        max[i] =  5.12;
        min[i] = -5.12;
    }

    de_Init_Uniform(&prob_data, ndim, npop, min, max, DE_MINIMIZE);

    switch (algo) {
    case PCX:
        de_SetupPCX(&prob_data, PCX_PROB, STD_F, STD_XOVER, PCX_STRATEGY, Rosenbrock2D, PCX_SIGMA2_A, PCX_SIGMA2_B);
        de_SetupRandomFK(&prob_data, F_RND_FACT, F_POP_FACT, F_RND_FACT, F_POP_FACT);
    case STD:
        de_SetupStd(&prob_data, STD_F, STD_XOVER, STD_STRATEGY, Rosenbrock2D);
        de_SetupRandomFK(&prob_data, F_RND_FACT, F_POP_FACT, F_RND_FACT, F_POP_FACT);
        break;
    case EITHER_OR:
        de_SetupEitherOr(&prob_data, 0.6, 0.8, 0.01, Rosenbrock2D);
        de_SetupRandomFK(&prob_data, 0.05, F_POP_FACT, 0.05, F_POP_FACT);
    }

    printf("\nRosenbrock's 2D : Calculating...\n");
    begin = monotime();
    de_Solve(&prob_data, ngen);
    end = monotime();

    print_coeffs(&prob_data);

    de_Finalize(&prob_data);

    return (end - begin);  
}

double t_ackley(de_algo_t algo) {
    double min[MAX_DIM];
    double max[MAX_DIM];
    int i;
    de_t prob_data;
    double begin, end;

    int ndim = 100;
    int npop = DIM100_POP;
    int ngen = 10000;

    for (i=0; i<ndim; i++) {
        max[i] =  32.768;
        min[i] = -32.768;
    }

    de_Init_Uniform(&prob_data, ndim, npop, min, max, DE_MINIMIZE);

    switch (algo) {
    case PCX:
        de_SetupPCX(&prob_data, 0.02, STD_F, STD_XOVER, PCX_STRATEGY, Ackley, PCX_SIGMA2_A, PCX_SIGMA2_B);
        de_SetupRandomFK(&prob_data, F_RND_FACT, F_POP_FACT, F_RND_FACT, F_POP_FACT);
        break;
    case STD:
        de_SetupStd(&prob_data, STD_F, STD_XOVER, STD_STRATEGY, Ackley);
        de_SetupRandomFK(&prob_data, F_RND_FACT, F_POP_FACT, F_RND_FACT, F_POP_FACT);
        break;
    case EITHER_OR:
        eprintf("Not implemented\n");
        return 0.0;
    }
    printf("\nAckley's : Calculating...\n");
    begin = monotime();
    de_Solve(&prob_data, ngen);
    end = monotime();

    print_coeffs(&prob_data);

    de_Finalize(&prob_data);

    return (end - begin);  
}

double t_griewangk(de_algo_t algo) {
    double min[MAX_DIM];
    double max[MAX_DIM];
    int i;
    de_t prob_data;
    double begin, end;

    int ndim = 100;
    int npop = DIM100_POP;
    int ngen = 10000;

    for (i=0; i<ndim; i++) {
        max[i] =  600.0;
        min[i] = -600.0;
    }

    de_Init_Uniform(&prob_data, ndim, npop, min, max, DE_MINIMIZE);

    switch (algo) {
    case PCX:
        de_SetupPCX(&prob_data, PCX_PROB, STD_F, STD_XOVER, PCX_STRATEGY, Griewangk, PCX_SIGMA2_A, PCX_SIGMA2_B);
        de_SetupRandomFK(&prob_data, F_RND_FACT, F_POP_FACT, F_RND_FACT, F_POP_FACT);
        break;
    case STD:
        de_SetupStd(&prob_data, STD_F, STD_XOVER, STD_STRATEGY, Griewangk);
        de_SetupRandomFK(&prob_data, F_RND_FACT, F_POP_FACT, F_RND_FACT, F_POP_FACT);
        break;
    case EITHER_OR:
        de_SetupEitherOr(&prob_data, 0.6, 0.8, 0.01, Griewangk);
        de_SetupRandomFK(&prob_data, 0.05, F_POP_FACT, 0.05, F_POP_FACT);
    }
    printf("\nGriewangk's : Calculating...\n");
    begin = monotime();
    de_Solve(&prob_data, ngen);
    end = monotime();

    print_coeffs(&prob_data);

    de_Finalize(&prob_data);

    return (end - begin);  
}

double t_rastrigin(de_algo_t algo) {
    double min[MAX_DIM];
    double max[MAX_DIM];
    int i;
    de_t prob_data;
    double begin, end;

    int ndim = 100;
    int npop = DIM100_POP;
    int ngen = DIM100_GEN;

    for (i=0; i<ndim; i++) {
        max[i] =  5.12;
        min[i] = -5.12;
    }

    de_Init_Uniform(&prob_data, ndim, npop, min, max, DE_MINIMIZE);

    switch (algo) {
    case PCX:
        de_SetupPCX(&prob_data, PCX_PROB, STD_F, STD_XOVER, PCX_STRATEGY, Rastrigin, PCX_SIGMA2_A, PCX_SIGMA2_B);
        de_SetupRandomFK(&prob_data, F_RND_FACT, F_POP_FACT, F_RND_FACT, F_POP_FACT);
        break;
    case STD:
        de_SetupStd(&prob_data, STD_F, STD_XOVER, STD_STRATEGY, Rastrigin);
        de_SetupRandomFK(&prob_data, F_RND_FACT, F_POP_FACT, F_RND_FACT, F_POP_FACT);
        break;
    case EITHER_OR:
        eprintf("Not implemented\n");
        return 0.0;
    }
    printf("\nRastrigin's : Calculating...\n");
    begin = monotime();
    de_Solve(&prob_data, ngen);
    end = monotime();

    print_coeffs(&prob_data);

    de_Finalize(&prob_data);

    return (end - begin);  
}

double t_schwefel(de_algo_t algo) {
    double min[MAX_DIM];
    double max[MAX_DIM];
    int i;
    de_t prob_data;
    double begin, end;

    int ndim = 100;
    int npop = DIM100_POP;
    int ngen = 10000;

    for (i=0; i<ndim; i++) {
        max[i] =  500.0;
        min[i] = -500.0;
    }

    de_Init_Uniform(&prob_data, ndim, npop, min, max, DE_MINIMIZE);

    switch (algo) {
    case PCX:
        de_SetupPCX(&prob_data, 0.02, STD_F, STD_XOVER, PCX_STRATEGY, Schwefel, PCX_SIGMA2_A, PCX_SIGMA2_B);
        de_SetupRandomFK(&prob_data, F_RND_FACT, F_POP_FACT, F_RND_FACT, F_POP_FACT);
        break;
    case STD:
        de_SetupStd(&prob_data, STD_F, STD_XOVER, Rand1Exp, Schwefel);
        de_SetupRandomFK(&prob_data, F_RND_FACT, F_POP_FACT, F_RND_FACT, F_POP_FACT);
        break;
    case EITHER_OR:
        eprintf("Not implemented\n");
        return 0.0;    
    }
    printf("\nSchwefel's : Calculating...\n");
    begin = monotime();
    de_Solve(&prob_data, ngen);
    end = monotime();

    print_coeffs(&prob_data);

    de_Finalize(&prob_data);

    return (end - begin);  
}

int main(int argc, char *argv[]) {
    double el[20]; /* elapsed */


    el[0] = t_rosenbrock(STD);
    dprintf("elapsed: %f s\n", el[0]);

    el[1] = t_rosenbrock(PCX);
    dprintf("elapsed: %f s\n", el[1]);

    el[2] = t_rosenbrock(EITHER_OR);
    dprintf("elapsed: %f s\n", el[2]);

    /***/
    el[3] = t_ackley(STD);
    dprintf("elapsed: %f s\n", el[3]);

    el[4] = t_ackley(PCX);
    dprintf("elapsed: %f s\n", el[4]);

    /***/
    // el[5] = t_griewangk(STD);
    // dprintf("elapsed: %f s\n", el[5]);

    // el[6] = t_griewangk(PCX);
    // dprintf("elapsed: %f s\n", el[6]);

    // el[7] = t_griewangk(EITHER_OR);
    // dprintf("elapsed: %f\n", el[7]);

    /***/
    // el[8] = t_rastrigin(STD);
    // dprintf("elapsed: %f s\n", el[8]);

    // el[9] = t_rastrigin(PCX);
    // dprintf("elapsed: %f s\n", el[9]);

    /***/
    el[10] = t_schwefel(STD);
    dprintf("elapsed: %f s\n", el[10]);

    el[11] = t_schwefel(PCX);
    dprintf("elapsed: %f s\n", el[11]);


    /* summary */
    dprintf("\nRuntime summary\n");
    dprintf("Rosenbrock's 2D:\n");
    dprintf("\tStd: %.3E [s]\n", el[0]);
    dprintf("\tPCX: %.3E [s]\n", el[1]);
    dprintf("\tE/O: %.3E [s]\n\n", el[2]);

    dprintf("Ackley's 100D:\n");
    dprintf("\tStd: %.3E [s]\n", el[3]);
    dprintf("\tPCX: %.3E [s]\n\n", el[4]);

    /*dprintf("Griewangk's 100D:\n");
    dprintf("\tStd: %.3E [s]\n", el[5]);
    dprintf("\tPCX: %.3E [s]\n", el[6]);
    //dprintf("\tE/O: %.3E [s]\n", el[7]);
    dprintf("\n");

    dprintf("Rastrigin's 100D:\n");
    dprintf("\tStd: %.3E [s]\n", el[8]);
    dprintf("\tPCX: %.3E [s]\n\n", el[9]);*/

    dprintf("Schwefel's 100D:\n");
    dprintf("\tStd: %.3E [s]\n", el[10]);
    dprintf("\tPCX: %.3E [s]\n", el[11]);


    return 0;
}