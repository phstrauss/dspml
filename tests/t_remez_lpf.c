#include "remez.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

main() {
   double *weights, *desired, *transitions;
   double *h;
   int i, ret;
   int n = 101;

   transitions   = (double *) malloc(4 * sizeof(double));
   weights = (double *) malloc(2 * sizeof(double));
   desired = (double *) malloc(2 * sizeof(double));
   h       = (double *) malloc(n * sizeof(double));

   desired[0] = 1;
   desired[1] = 0;

   weights[0] = 1;
   weights[1] = 1;

   transitions[0] = 0;
   transitions[1] = 0.1;
   transitions[2] = 0.2;
   transitions[3] = 0.5;

   ret = remez_from_bands(h, n, 2, transitions, desired, weights, STANDARD, 0.0001);

   printf("return from remez: %d (0 = failed, no convergence, otherwise number of iterations)\n", ret);
   for (i=0; i<n; i++)
       printf("%23.20f\n", h[i]);

   free(transitions);
   free(weights);
   free(desired);
   free(h);
}
