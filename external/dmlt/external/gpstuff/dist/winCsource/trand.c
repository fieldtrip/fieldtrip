/*TRAND  Generates random numbers from Students t-distribution
 *
 *   R = TRAND(V) returns a matrix of random numbers chosen   
 *   from the T distribution with V degrees of freedom.
 *   The size of R is the size of V.
 *   Alternatively, R = TRAND(V,M,N) returns an M by N matrix. 
 */


/* Copyright (C) 1998-2000 Aki Vehtari
 * 
 *This software is distributed under the GNU General Public 
 *License (version 3 or later); please refer to the file 
 *License.txt, included with the software, for details.
 *
 */

#include <math.h>
#include "mex.h"
#include "rand.h"

void mexFunction(const int nlhs, mxArray *plhs[],
		 const int nrhs, const mxArray *prhs[])
{
  if (nlhs > 1 ) 
    mexErrMsgTxt( "Too many output arguments." );
  
  if ((nrhs < 1) | (nrhs == 2) | (nrhs > 3))
    mexErrMsgTxt( "Wrong number of input arguments." );
  
  {
    mxArray *MN[2];
    double nu, *nup, a, *r;
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
    mexCallMATLAB(1,&plhs[0],2,MN,"randn");
    mxDestroyArray(MN[0]);
    mxDestroyArray(MN[1]);

    r = mxGetPr(plhs[0]);
    mn=m*n;
    
    if (dims[0] > 1 || dims[1] > 1) {
      nup=mxGetPr(prhs[0]);
      for (i = 0; i < mn; i++) {
	nu=nup[i];
	if (nu == 1) {
	  /* Cauchy distribution */
	  r[i] = r[i]/rand_gaussian();
	}
	else if (nu == 2) {
	  a=rand_uniopen();
	  r[i] = sqrt(2.0)*(a-0.5)/sqrt(a-a*a);
	}
	else {
	  a=sqrt(nu);
	  r[i] = r[i]*a/sqrt(rand_gamma(nu/2)*2);
	}
      }
    }
    else {
      nu=mxGetScalar(prhs[0]);
      if (nu == 1) {
	for (i = 0; i < mn; i++) {
	  /* Cauchy distribution */
	  r[i] = r[i]/rand_gaussian();
	}
      }
      else if (nu == 2) {
	for (i = 0; i < mn; i++) {
	  a=rand_uniopen();
	  r[i] = sqrt(2.0)*(a-0.5)/sqrt(a-a*a);
	}
      }
      else {
	a=sqrt(nu);
	for (i = 0; i < mn; i++) {
	  r[i] = r[i]*a/sqrt(rand_gamma(nu/2)*2);
	}
      }
    }
  }
  return;
}     

