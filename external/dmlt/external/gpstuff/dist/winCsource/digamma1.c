#include <math.h>
#include "mex.h"

/* DIGAMMA1.C - Compute the digamma function. */

double digamma (double);

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{

  if (nrhs != 1 )
    mexErrMsgTxt( "Wrong number of input arguments." );
  
  if (nlhs > 1 )
    mexErrMsgTxt( "Too many output arguments." );
  
  {
    double *x, *y; 
    const int *dims;
    int len, i;

    dims = mxGetDimensions(prhs[0]);
    len = dims[0]*dims[1];
    plhs[0] = mxCreateDoubleMatrix(dims[0], dims[1], mxREAL); 
    x = mxGetPr(prhs[0]);
    y = mxGetPr(plhs[0]);

    for (i = 0; i < len; i++) {
      y[i] = digamma(x[i]);
    }
  }
  
  return;
  
}

/* DIGAMMA.C - Compute the digamma function. */

/* Copyright (c) 1995-2004 by Radford M. Neal 
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

/* digamma(x) is defined as (d/dx) log Gamma(x).  It is computed here
   using an asymptotic expansion when x>5.  For x<=5, the recurrence
   relation digamma(x) = digamma(x+1) - 1/x is used repeatedly.  See
   Venables & Ripley, Modern Applied Statistics with S-Plus, pp. 151-152. */


/* COMPUTE THE DIGAMMA FUNCTION.  Returns -inf if the argument is an integer 
   less than or equal to zero. */

double digamma (double x)
{
  double r, f, t;

  r = 0;

  while (x<=5) 
  { r -= 1/x;
    x += 1;
  }

  f = 1/(x*x);

  t = f*(-1/12.0 + f*(1/120.0 + f*(-1/252.0 + f*(1/240.0 + f*(-1/132.0
       + f*(691/32760.0 + f*(-1/12.0 + f*3617/8160.0)))))));

  return r + log(x) - 0.5/x + t;
}
