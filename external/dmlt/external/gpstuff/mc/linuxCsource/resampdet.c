/* 
 * RESAMPDET Deterministic resampling
 *
 *  Description
 *  S = RESAMPDET(P) returns a new set of indices according to the
 *  probabilities P. P is array of probabilities, which are not
 *  necessarily normalized, though they must be non-negative, and
 *  not all zero. The size of S is the size of P. 
 *
 *  S = RESAMPDET(P,M,N) returns M by N matrix.
 *
 *  S = RESAMPDET(P,M) returns M by M matrix.
 *
 *  Default is to use no-sort resampling. For sorted resampling use
 *   [PS,PI]=SORT(P);
 *   S=PI(RESAMPDET(PS));
 *  Sorted re-sampling is slower but has smaller variance. Note
 *  that deterministic resampling is not unbiased. Stratified
 *  resampling (RESAMPSTR) is unbiased, almost as fast as
 *  deterministic resampling, and has only slightly larger
 *  variance.
 *
 *  In deterministic resampling indices are sampled using
 *  deterministic numbers u_j~(j-a)/n, for fixed a in [0,1) and
 *  n is length of P. Compare this to simple random resampling
 *  where u_j~U[0,n]. See, Kitagawa, G., Monte Carlo Filter and
 *  Smoother for Non-Gaussian Nonlinear State Space Models,
 *  Journal of Computational and Graphical Statistics, 5(1):1-25,
 *  1996. 
 *
 *  See also RESAMPSIM, RESAMPRES, RESAMPSTR
 *
 * Last modified: 2003-03-20 12:53:57 EET
 *
 */

/* Copyright (C) 2003 Aki Vehtari
 * 
 *This software is distributed under the GNU General Public 
 *License (version 3 or later); please refer to the file 
 *License.txt, included with the software, for details.
 *
 */

#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{

  if (nlhs > 1 ) 
    mexErrMsgTxt( "Too many output arguments." );
  
  if (nrhs < 1 || nrhs > 3)
    mexErrMsgTxt( "Wrong number of input arguments." );

  {
    double *p, *q, *s, sum=0, b, c=0.5; 
    const int *dims;
    int i, len, a, n, m, mn;

    dims = mxGetDimensions(prhs[0]);
    len = dims[0]*dims[1];
    p = mxGetPr(prhs[0]);
    q = mxMalloc(len*sizeof(double));

    if (nrhs==1) {
      m = dims[0];
      n = dims[1];
    }
    else { /* nrhs > 1 */
      dims = mxGetDimensions(prhs[1]);
      if (dims[0]>1 || dims[1]>1)
	    mexErrMsgTxt("M must be a scalar.");
        m = mxGetScalar(prhs[1]);
        if (nrhs==2) {
	      n = m;
        }	
        else {
	      dims = mxGetDimensions(prhs[2]);
	      if (dims[0]>1 || dims[1]>1)
	        mexErrMsgTxt("M must be a scalar.");
	        n = mxGetScalar(prhs[2]);
        }
    }
    mn = m*n;
    
    plhs[0]=mxCreateDoubleMatrix(m,n,mxREAL);
    s = mxGetPr(plhs[0]);

    /* integer part */
    for (i = 0; i < len; i++)
      sum+=p[i];
    for (i = 0; i < len; i++) {
      q[i]=p[i]/sum*mn;
      if (q[i]>=1.0) {
	a=(int)(q[i]);
	q[i]-=a;
	for (b=i+1.0;a>0;a--)
	  *s++=b;
      }
      c+=q[i];
      if (c>=1.0) {
	*s++=i+1.0;
	c-=1.0;
      }
    }
  }

  return;
}
