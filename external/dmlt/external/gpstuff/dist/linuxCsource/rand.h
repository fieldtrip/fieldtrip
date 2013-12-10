/* RAND.H - Interface to random number generation procedures. */


/* Copyright (c) 1995-2003 by Radford M. Neal 
 *
 * Permission is granted for anyone to copy, use, modify, or distribute this
 * program and accompanying programs and documents for any purpose, provided 
 * this copyright notice is retained and prominently displayed, along with
 * a note saying that the original programs are available from Radford Neal's
 * web page, and note is made of any changes made to the programs.  The
 * programs and documents are distributed without any warranty, express or
 * implied.  As the programs were written for research purposes only, they have
 * not been tested to the degree that would be advisable in any important
 * application.  All use of these programs is entirely at the user's own risk.
 */

/* GENERATORS FOR VARIOUS DISTRIBUTIONS. */

double rand_uniform (void);	/* Uniform from [0,1) */
double rand_uniopen (void);	/* Uniform from (0,1) */

int rand_int (int);		/* Uniform from 0, 1, ... (n-1) */
int rand_pickd (double *, int);	/* From 0 ... (n-1), with given distribution */
int rand_pickf (float *, int);	/* Same as above, but with floats */

double rand_gaussian (void);	/* Gaussian with mean zero and unit variance */
double rand_exp (void);		/* Exponential with mean one */
double rand_cauchy (void);	/* Cauchy centred at zero with unit width */
double rand_gamma (double);	/* Gamma with given shape parameter */
double rand_beta (double, double); /* Beta with given parameters */
