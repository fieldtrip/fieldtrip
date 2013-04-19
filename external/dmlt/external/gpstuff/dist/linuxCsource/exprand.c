/* 
 *    EXPRAND Random matrices from exponential distribution.
 *
 *    R = EXPRAND(MU) returns a matrix of random numbers chosen   
 *    from the exponential distribution with parameter MU.
 *    The size of R is the size of MU.
 *    Alternatively, R = EXPRAND(MU,M,N) returns an M by N matrix. 
 *
 * Last modified: 2011-09-07 10:43:19 EEST
 *
 */

#include <math.h>
#include "mex.h"

/* Copyright (C) 1998-2000 Aki Vehtari
 * 
 *This software is distributed under the GNU General Public 
 *License (version 3 or later); please refer to the file 
 *License.txt, included with the software, for details.
 *
 */

void mexFunction(const int nlhs, mxArray *plhs[],
		 const int nrhs, const mxArray *prhs[])
{
  if (nlhs > 1 ) 
    mexErrMsgTxt( "Too many output arguments." );
  
  if ((nrhs < 1) | (nrhs == 2) | (nrhs > 3))
    mexErrMsgTxt( "Wrong number of input arguments." );
  
  {
    mxArray *MN[2];
    double mu, *mup, *r;
    const int *dims;
    int i, m, n, mn;
    
    dims = mxGetDimensions(prhs[0]);
    if (nrhs == 1) {
      m = dims[0];
      n = dims[1];
    }
    else {
      m=(int)mxGetScalar(prhs[1]);
      n=(int)mxGetScalar(prhs[2]);
      if ((dims[0] > 1 || dims[1] > 1) && dims[0]!=m && dims[1]!=n)
	mexErrMsgTxt( "Size information is inconsistent." );
    }
    MN[0]=mxCreateDoubleScalar((double)m);
    MN[1]=mxCreateDoubleScalar((double)n);
    mexCallMATLAB(1,&plhs[0],2,MN,"rand");
    mxDestroyArray(MN[0]);
    mxDestroyArray(MN[1]);

    r = mxGetPr(plhs[0]);
    mn=m*n;
    
    if (dims[0] > 1 || dims[1] > 1) {
      mup=mxGetPr(prhs[0]);
      for (i = 0; i < mn; i++) 
	r[i] = -mup[i]*log(r[i]);
    }
    else {
      mu=mxGetScalar(prhs[0]);
      if (mu==1) {
	for (i = 0; i < mn; i++) 
	  r[i] = -log(r[i]);
      }
      else {
	for (i = 0; i < mn; i++) 
	  r[i] = -mu*log(r[i]);
      }
    }
  }
  return;
}     

