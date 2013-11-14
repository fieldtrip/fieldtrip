/* RAND.C - Random number generation module. */

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

#include <math.h>
#include "mex.h"
#include "rand.h"

/*
 * 2000-05-31  Aki Vehtari  <Aki.Vehtari@hut.fi>
 *  * Use Matlab's rand and randn functions
 *  * rand_uniform and rand_uniopen both generate uniformly form (0,1)
 */

/*    
   Many of the methods used in this module may be found in the following
   reference:

      Devroye, L. (1986) Non-Uniform Random Variate Generation, 
        New York: Springer-Verlag.

   The methods used here are not necessarily the fastest available.  They're
   selected to be reasonably fast while also being easy to write.
*/


/* CONSTANT PI.  Defined here if not in <math.h>. */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* GENERATE UNIFORMLY FROM [0,1). */
/* NOT... This too currently generates uniformly from (0,1) !!!
   Maybe someday I'll fix this..? */
double rand_uniform (void)
{
  mxArray *R;
  mexCallMATLAB(1,&R,0,NULL,"rand");
  return(mxGetScalar(R));
}


/* GENERATE UNIFORMLY FORM (0,1). */
double rand_uniopen (void)
{
  mxArray *R;
  mexCallMATLAB(1,&R,0,NULL,"rand");
  return(mxGetScalar(R));
}


/* GENERATE RANDOM INTEGER FROM 0, 1, ..., (n-1). */
int rand_int (int n)
{ 
  return (int) (n * rand_uniform());
}


/* GENERATE INTEGER FROM 0, 1, ..., (n-1), WITH GIVEN DISTRIBUTION.  The
   distribution is given by an array of probabilities, which are not 
   necessarily normalized, though they must be non-negative, and not all 
   zero. */
int rand_pickd (double *p, int n)
{ 
  double t, r;
  int i;

  t = 0;
  for (i = 0; i<n; i++)
  { if (p[i]<0) abort();
    t += p[i];
  }

  if (t<=0) abort();

  r = t * rand_uniform();

  for (i = 0; i<n; i++)
  { r -= p[i];
    if (r<0) return i;
  }

  /* Return value with non-zero probability if we get here due to roundoff. */

  for (i = 0; i<n; i++) 
  { if (p[i]>0) return i;
  }

  abort(); 
}


/* SAME PROCEDURE AS ABOVE, BUT WITH FLOAT ARGUMENT. */

int rand_pickf (float *p, int n)
{ 
  double t, r;
  int i;

  t = 0;
  for (i = 0; i<n; i++)
  { if (p[i]<=0) abort();
    t += p[i];
  }

  if (t<=0) abort();

  r = t * rand_uniform();

  for (i = 0; i<n; i++)
  { r -= p[i];
    if (r<0) return i;
  }

  /* Return value with non-zero probability if we get here due to roundoff. */

  for (i = 0; i<n; i++) 
  { if (p[i]>0) return i;
  }

  abort(); 
}


/* GAUSSIAN GENERATOR.
   Done by using the Matlab's randn function. */
double rand_gaussian (void)
{
  mxArray *R;
  mexCallMATLAB(1,&R,0,NULL,"randn");
  return(mxGetScalar(R));
}


/* EXPONENTIAL GENERATOR.  See Devroye, p. 29.  Written so as to never
   return exactly zero. */
double rand_exp (void)
{
  return -log(rand_uniopen());
}


/* CAUCHY GENERATOR.  See Devroye, p. 29. */
double rand_cauchy (void)
{
  return tan (M_PI * (rand_uniopen()-0.5));
}


/* GAMMA GENERATOR.  Generates a positive real number, r, with density
   proportional to r^(a-1) * exp(-r).  See Devroye, p. 410 and p. 420. 
   Things are fiddled to avoid ever returning a value that is very near 
   zero. */
double rand_gamma (double a)
{
  double b, c, X, Y, Z, U, V, W;

  if (a<0.00001)
  { X = a;
  }

  else if (a<=1) 
  { 
    U = rand_uniopen();
    X = rand_gamma(1+a) * pow(U,1/a);
  }

  else if (a<1.00001)
  { X = rand_exp();
  }

  else
  {
    b = a-1;
    c = 3*a - 0.75;
  
    for (;;)
    {
      U = rand_uniopen();
      V = rand_uniopen();
    
      W = U*(1-U);
      Y = sqrt(c/W) * (U-0.5);
      X = b+Y;
  
      if (X>=0)
      { 
        Z = 64*W*W*W*V*V;
  
        if (Z <= 1 - 2*Y*Y/X || log(Z) <= 2 * (b*log(X/b) - Y)) break;
      }
    }
  }

  return X<1e-15 && X<a ? (a<1e-15 ? a : 1e-15) : X;
}


/* BETA GENERATOR. Generates a real number, r, in (0,1), with density
   proportional to r^(a-1) * (1-r)^(b-1).  Things are fiddled to avoid
   the end-points, and to make the procedure symmetric between a and b. */
double rand_beta 
( double a, 
  double b
)
{
  double x, y, r;

  do
  { x = rand_gamma(a);
    y = rand_gamma(b);
    r = 1.0 + x/(x+y);
    r = r - 1.0;
  } while (r<=0.0 || r>=1.0);

  return r;
}
