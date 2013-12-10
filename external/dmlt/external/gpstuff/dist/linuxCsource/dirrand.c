/* * DIRRAND - Uniform Dirichlet random vectors
 * Generates random vectors from uniform Dirichlet distribution (1,1,...,1)
 *
 *
 *   R = DIRRAND(M) returns vector of length M chosen from
 *   Dirichlet(1,...,1) distribution.
 *   R = DIRRAND(M,N) returns N such vectors in MxN matrix.
 *
 *   Sampling from Dirichlet is done using exponentials. The
 *   exponentials are generated using a simple inversion method.
 *   See Gelman et al., pages 480,482
 */

/*
 * Copyright (C) 1998-2003 Aki Vehtari
 *
 * Last modified: 2003-03-21 13:45:01 EET
 *
 */

/* 
 *This software is distributed under the GNU General Public 
 *License (version 3 or later); please refer to the file 
 *License.txt, included with the software, for details.
 *
 */

#include <math.h>
#include "mex.h"

void mexFunction(const int nlhs, mxArray *plhs[],
		 const int nrhs, const mxArray *prhs[])
{
  if (nlhs > 1 ) 
    mexErrMsgTxt( "Too many output arguments.");
  
  if (nrhs < 1 || nrhs > 2)
    mexErrMsgTxt( "Wrong number of input arguments." );
  
  {
    mxArray *MN[2];
    double *p=NULL, *r, rsum;
    const int *dims;
    int i, j, m, n, mn;
    
    dims = mxGetDimensions(prhs[0]);
    if (dims[0] > 1 || dims[1]>1)
      mexErrMsgTxt( "First argument must be a scalar." );
    m = mxGetScalar(prhs[0]);
    
    if (nrhs>1) {
      /* second argument tells how many vectors we want */
      dims = mxGetDimensions(prhs[1]);
      if (dims[0] > 1 || dims[1]>1)
	mexErrMsgTxt( "Second argument must be a scalar." );
      n = mxGetScalar(prhs[1]);
    } else {
      n=1;
    }
    
    mn = m*n;
    MN[0]=mxCreateDoubleScalar((double)m);
    MN[1]=mxCreateDoubleScalar((double)n);
    mexCallMATLAB(1,&plhs[0],2,MN,"rand");
    mxDestroyArray(MN[0]);
    mxDestroyArray(MN[1]);
    r = mxGetPr(plhs[0]);
    
    for (j = 0; j < mn; j+=m) {
      for (i = 0, rsum=0; i < m; i++) {
	r[j+i] = -log(r[j+i]);
	rsum += r[j+i];
      }
      for (i = 0; i < m; i++) 
	r[j+i] /= rsum;
    }
    
  }
  
  return;
}     
