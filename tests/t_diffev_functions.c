/*  Differential evolution solver heurisitc
    Lester E. Godwin
    Partially rewrote, enhanced (EitherOr and Probabilist Parent Centric evolution strategies)
    plus maintained by Philippe Strauss, Spring 2012.  */


#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include "de.h"
#include "perr.h"
#include "monotime.h"
#include "optim_bench.h"


const double epsilon = 1e-8;


double Rosenbrock2D(int n, double bestEnergy, double *trial, bool *bAtSolution) {
    return (pow(1.0 - trial[0], 2.0) + 100.0 * pow(trial[1] - pow(trial[0], 2.0), 2.0)); 
}

/*  Ackley's
    narrow, needle-like global minima
    y = 0.0 for xn = 0.0 for all n's */
double Ackley(int n, double bestEnergy, double *trial, bool *bAtSolution) {
    int i;
    double sumx2 = 0.0;
    double sumcos = 0.0;
    double c = 2 * M_PI;
    double a = 20.0;
    double b = 0.2;
    double ret;

    for (i = 0; i < n; ++i) {
        sumx2 += pow(trial[i], 2.0);
        sumcos += cos(c * trial[i]);
    }
    sumx2 /= (double) n;
    sumcos /= (double) n;
    ret = -a * exp(-b * sqrt(sumx2)) - exp(sumcos) + a + exp(1.0);

    if (fabs(ret) < epsilon)
        *bAtSolution = true;

    return ret;
}

/*  Griewangk's
    y = 0.0 for xn = 0.0 for all n's */
double Griewangk(int n, double bestEnergy, double *trial, bool *bAtSolution) {
    int i;
    double sumx2 = 0.0;
    double prodcos = 0.0;
    double ret;

    for (i = 0; i < n; ++i) {
        prodcos *= cos(trial[i] / sqrt((double) (i+1)));
        sumx2 += pow(trial[i], 2.0);
    }
    ret = 1.0 + (sumx2 / 4000.0) - prodcos;

    if (fabs(ret) < (1.0 + epsilon))
        *bAtSolution = true;

    return ret; 
}

double Rastrigin(int n, double bestEnergy, double *trial, bool *bAtSolution) {
    int i;
    double sum = 0.0;
    double ret;

    for (i=0; i < n; ++i)
        sum += pow(trial[i], 2.0) - 10.0 * cos(2.0 * M_PI * trial[i]);

    ret = sum + n * 10.0;

    if (fabs(ret) < epsilon)
        *bAtSolution = true;

    return ret;
}

/*  Schwefel
    deceptive local minimas far from global */
double Schwefel(int n, double bestEnergy, double *trial, bool *bAtSolution) {
    int i;
    double sum = 0.0;
    double ret;

    for (i = 0; i < n; ++i)
        sum += -1.0 * trial[i] * sin(sqrt(fabs(trial[i])));

    /* see http://www.it.lut.fi/ip/evo/functions/node10.html */
    ret = 418.9828872724 * n + sum;

    if (fabs(ret) < epsilon)
        *bAtSolution = true;

    return ret;
}

/*  Michalewicz's
    c.f.: GEATbx_ObjFunExpl_v37.pdf
    solutions bound's : 0.0 -> Pi; global minima of ~ -4.687 for 5 dims, ~ -9.66 for 10 dims */
double Michalewicz(int n, double bestEnergy, double *trial, bool *bAtSolution) {
    int i;
    double sum = 0.0;
    double m = 10.0;

    for(i = 0; i < n; ++i) {
        sum += sin(trial[i]) * pow(sin(((double) i+1) * pow(trial[i], 2.0)/ M_PI), 2.0 * m);
    }

    return (-sum);
}

/*  Langermann's :

    http://extreme.adorio-research.org/download/mvf/html/node29.html :

    The matrix A and column vector c are identical for Shekel's foxhole mvfFoxhole function except
    that c2 = -1.5 for Langerman's function. Langerman's function has a global minimum value of -1.4

    http://www.cs.cmu.edu/afs/cs/project/jair/pub/volume24/ortizboyer05a-html/node23.html :

    BDL+96
    H. Bersini, M. Dorigo, S. Langerman, G. Seront, and L. M. Gambardella. 
    Results of the first international contest on evolutionary optimisation (1st iceo). 
    In Proceedings of IEEE International Conference on Evolutionary Computation,
    IEEE-EC 96, pages 611-615, Nagoya, Japan, May 20-22 1996. IEEE Press.  */
