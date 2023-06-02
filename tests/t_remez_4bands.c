/**************************************************************************
 * Parks-McClellan algorithm for FIR filter design (C version)
 *-------------------------------------------------
 *  Copyright (C) 1995  Jake Janovetz (janovetz@coewl.cen.uiuc.edu)
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 *************************************************************************/

/***************************************************************************
 * Test program for the remez() function.  Sends appropriate arguments to
 * remez() to generate a filter.  Then, prints the resulting coefficients
 * to stdout.
 **************************************************************************/

#include "remez.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

main()
{
   double *weights, *desired, *transitions;
   double *h;
   int i, ret;

   transitions = (double *)malloc(10 * sizeof(double));
   weights = (double *)malloc(5 * sizeof(double));
   desired = (double *)malloc(5 * sizeof(double));
   h = (double *)malloc(300 * sizeof(double));

   desired[0] = 0;
   desired[1] = 1; 
   desired[2] = 0;
   desired[3] = 1;
   desired[4] = 0;

   weights[0] = 10;
   weights[1] = 1;
   weights[2] = 3;
   weights[3] = 1;
   weights[4] = 20;

   transitions[0] = 0;
   transitions[1] = 0.05;
   transitions[2] = 0.1;
   transitions[3] = 0.15;
   transitions[4] = 0.18;
   transitions[5] = 0.25;
   transitions[6] = 0.3;
   transitions[7] = 0.36;
   transitions[8] = 0.41;
   transitions[9] = 0.5;

   ret = remez_from_bands(h, 104, 5, transitions, desired, weights, STANDARD, 0.0001);
   printf("return from remez: %d (0 = failed, no convergence, otherwise number of iterations)\n", ret);
   for (i=0; i<104; i++)
       printf("%23.20f\n", h[i]);

   free(transitions);
   free(weights);
   free(desired);
   free(h);
}

