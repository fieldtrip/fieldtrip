/* BBPRCTILE  Bayesian bootstrap percentile
 *
 *   Description:
 *   bbp = bbprctile(x,p,B,w)
 *   x   = Mx1 matrix of data
 *   p   = Px1 or 1xP vector of percentiles
 *   B   = number of bootstrap replicates
 *   w   = Mx1 vector of weights 

 *
 * Copyright (C) 1998-2003 Aki Vehtari
 *
 * Last modified: 2011-09-07 10:36:44 EEST
 *
 */


/* 
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
  
  if (nrhs < 2 || nrhs > 4)
    mexErrMsgTxt( "Wrong number of input arguments." );
  
  {
    double *x, *p, *w=NULL, *r, rsum, *bbp;
    mxArray *R, *P, *MN[2];
    const int *dims;
    int B, b, i, j, jj, m, n, plen;
    
    dims = mxGetDimensions(prhs[0]);
    x = mxGetPr(prhs[0]);
    m = dims[0];
    n = dims[1];
    
    dims = mxGetDimensions(prhs[1]);
    if (dims[0]>1 && dims[1]>1)
      mexErrMsgTxt("P must be a vector.");
    plen = dims[0]*dims[1];
    mexCallMATLAB(1, &P, 1, (mxArray **) &prhs[1], "sort");
    p = mxGetPr(P);
    for (j = 0; j < plen; j++) {
      p[j]=p[j]/100;
    }
    
    if (nrhs>2) {
      dims = mxGetDimensions(prhs[2]);
      if (dims[0]>1 || dims[1]>1)
	mexErrMsgTxt("B must be a scalar.");
      B = mxGetScalar(prhs[2]);
    }
    else
      B = 1;
    
    if (nrhs>3) {
      dims = mxGetDimensions(prhs[3]);
      if (dims[1]>1)
	mexErrMsgTxt("W must be a column vector.");
      if (dims[0]!=m)
	mexErrMsgTxt("W must be a same size as X.");
      w = mxGetPr(prhs[3]);
    }

    plhs[0]=mxCreateDoubleMatrix(B, plen, mxREAL);
    bbp = mxGetPr(plhs[0]);

    MN[0]=mxCreateDoubleScalar((double)m);
    MN[1]=mxCreateDoubleScalar(1.0);

    if (w!=NULL) {
      for (b = 0; b < B; b++) {
	mexCallMATLAB(1,&R,2,MN,"rand");
	r = mxGetPr(R);
	jj=0;
	for (i = 0, rsum=0; i < m; i++) {
	  r[i] = -log(r[i])*w[i];
	  rsum += r[i];
	}
	r[0]=r[0]/rsum;
	for (j = 0; j < plen; j++) {
	  if (r[0] >= p[j]) {
	    bbp[b+j*B]=x[0];
	    jj++;
	  }
	  else {
	    break;
	  }
	}
	for (i = 1; i < m-1; i++) {
	  r[i]=r[i-1]+r[i]/rsum;
	  for (j = jj; j < plen; j++) {
	    if (r[i] >= p[j]) {
	      bbp[b+j*B]=x[i-1]+(x[i]-x[i-1])*(p[j]-r[i-1])/(r[i]-r[i-1]);
	      jj++;
	    }
	    else {
	      break;
	    }
	  }
	}
	r[i]=1;
	for (j = jj; j < plen; j++) {
	  if (r[i] >= p[j])
	    bbp[b+j*B]=x[i-1]+(x[i]-x[i-1])*(p[j]-r[i-1])/(r[i]-r[i-1]);
	  else {
	    break;
	  }
	}
      }
    }
    else {
      for (b = 0; b < B; b++) {
	mexCallMATLAB(1,&R,2,MN,"rand");
	r = mxGetPr(R);
	jj=0;
	for (i = 0, rsum=0; i < m; i++) {
	  r[i] = -log(r[i]);
	  rsum += r[i];
	}
	r[0]=r[0]/rsum;
	for (j = 0; j < plen; j++) {
	  if (r[0] >= p[j]) {
	    bbp[b+j*B]=x[0];
	    jj++;
	  }
	  else
	    break;
	}
	for (i = 1; i < m-1; i++) {
	  r[i]=r[i-1]+r[i]/rsum;
	  for (j = jj; j < plen; j++) {
	    if (r[i] >= p[j]) {
	      bbp[b+j*B]=x[i-1]+(x[i]-x[i-1])*(p[j]-r[i-1])/(r[i]-r[i-1]);
	      jj++;
	    }
	    else
	      break;
	  }
	}
	r[i]=1;
	for (j = jj; j < plen; j++) {
	  if (r[i] >= p[j])
	    bbp[b+j*B]=x[i-1]+(x[i]-x[i-1])*(p[j]-r[i-1])/(r[i]-r[i-1]);
	  else
	    break;
	}
	
      }
    }

    mxDestroyArray(MN[0]);
    mxDestroyArray(MN[1]);
    mxDestroyArray(R);
  }
    
  return;
}     
