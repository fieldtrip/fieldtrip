/*
 * TANH   75 % faster hyperbolic tangent.
 * TANH(X) is the hyperbolic tangent of the elements of X.
 *
 */

/* (c) 1998-2001 Aki Vehtari <Aki.Vehtari@hut.fi>
 * 
 *This software is distributed under the GNU General Public 
 *License (version 3 or later); please refer to the file 
 *License.txt, included with the software, for details.
 *
 */

#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  double *x, *y, xa; 
  const int *dims;
  int len, i;
  
  dims = mxGetDimensions(prhs[0]);
  len = dims[0]*dims[1];
  plhs[0] = mxCreateDoubleMatrix(dims[0], dims[1], mxREAL); 
  x = mxGetPr(prhs[0]);
  y = mxGetPr(plhs[0]);
  
  for (i = len; i > 0; i--, x++, y++) {
    xa=fabs(*x);
    *y = (xa > 18) ? (*x>0 ? 1.0 : -1.0)
      : ((xa > 1.5e-8) ? 1.0-2.0/(exp(2.0**x)+1.0) : *x);
  }
  
  return;
}
