/* BBMEAN     Bayesian bootstrap mean
 *
 *   Description
 *   bbm = bbmean(x,B,w)
 *   x MxN matrix of data
 *   B number of bootstrap replicates
 *   w Mx1 vector of weights 
 *
 * Last modified: 2011-09-07 10:47:40 EEST
 *
 */

/* Copyright (c) 1998-2004 Aki Vehtari
 *
 *This software is distributed under the GNU General Public 
 *License (version 3 or later); please refer to the file 
 *License.txt, included with the software, for details.
 *
 */

#include "mex.h"
#include <math.h>

void mexFunction(const int nlhs, mxArray *plhs[],
		 const int nrhs, const mxArray *prhs[])
{
  if (nlhs > 1 ) 
    mexErrMsgTxt( "Too many output arguments.");
  
  if (nrhs < 1 || nrhs > 3)
    mexErrMsgTxt( "Wrong number of input arguments." );
  
  {
    double *x, *w=NULL, *r, rsum, *bbm;
    mxArray *R, *MN[2];
    const int *dims;
    int B, b, i, j, m, n;
    
    dims = mxGetDimensions(prhs[0]);
    x = mxGetPr(prhs[0]);
    m = dims[0]; 
    n = dims[1];
    
    if (nrhs>1) {
      dims = mxGetDimensions(prhs[1]);
      if (dims[0]>1 || dims[1]>1)
	mexErrMsgTxt("B must be a scalar.");
      B = mxGetScalar(prhs[1]);
    }
    else {
      B = 1;
    }
    
    if (nrhs>2) {
      dims = mxGetDimensions(prhs[2]);
      if (dims[1]>1)
	mexErrMsgTxt("W must be a column vector.");
      if (dims[0]!=m)
	mexErrMsgTxt("W must have same length as X.");
      w = mxGetPr(prhs[2]);
    }

    plhs[0]=mxCreateDoubleMatrix(B, n, mxREAL);
    bbm = mxGetPr(plhs[0]);

    MN[0]=mxCreateDoubleScalar((double)m);
    MN[1]=mxCreateDoubleScalar(1.0);

    if (w!=NULL) {
      for (b = 0; b < B; b++) {
	mexCallMATLAB(1,&R,2,MN,"rand");
	r = mxGetPr(R);
	for (i = 0, rsum=0; i < m; i++) {
	  r[i] = -log(r[i])*w[i];
	  rsum += r[i];
	}
	for (i = 0; i < m; i++) {
	  r[i] /= rsum;
	  for (j = 0; j < n; j++) {
	    bbm[b+j*B] += x[i+j*m]*r[i];
	  }
	}
      }
    }
    else {
      for (b = 0; b < B; b++) {
	mexCallMATLAB(1,&R,2,MN,"rand");
	r = mxGetPr(R);
	for (i = 0, rsum=0; i < m; i++) {
	  r[i] = -log(r[i]);
	  rsum += r[i];
	}
	for (i = 0; i < m; i++) {
	  r[i] /= rsum;
	  for (j = 0; j < n; j++) {
	    bbm[b+j*B] += x[i+j*m]*r[i];
	  }
	}
      }
    }

    mxDestroyArray(MN[0]);
    mxDestroyArray(MN[1]);
    mxDestroyArray(R);

  }

  return;
}     
