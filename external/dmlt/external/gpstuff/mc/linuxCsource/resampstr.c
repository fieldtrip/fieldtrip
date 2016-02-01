/* Stratified Resampling No-Sort
 *
 * (C) 2003 Aki Vehtari
 *
 * Last modified: 2011-09-07 10:48:47 EEST
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
    mxArray *MN[2], *R;
    double *p, *q, *r, *s, sum=0, b, c=0; 
    const int *dims;
    int i, m, n, mn, pn, a, k=0;

    dims = mxGetDimensions(prhs[0]);
    pn = dims[0]*dims[1];            /* number of weights */
    p = mxGetPr(prhs[0]);            /* pointer to unormalized weights */
    q = mxMalloc(pn*sizeof(double)); /* pointer to normalized and scaled weights */

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
      }	else {
	dims = mxGetDimensions(prhs[2]);
	if (dims[0]>1 || dims[1]>1)
	  mexErrMsgTxt("M must be a scalar.");
	n = mxGetScalar(prhs[2]);
      }
    }
    mn = m*n;

    MN[0]=mxCreateDoubleScalar(m);
    MN[1]=mxCreateDoubleScalar(n);
    mexCallMATLAB(1,&R,2,MN,"rand");
    mxDestroyArray(MN[0]);
    mxDestroyArray(MN[1]);
    r = mxGetPr(R);           /* pointer to uniform random numbers */

    /* allocate matrix for return value */
    plhs[0]=mxCreateDoubleMatrix(m,n,mxREAL);
    s = mxGetPr(plhs[0]);     /* pointer to samples */

    for (i = 0; i < pn; i++)
      sum+=p[i];              /* compute the sum for normalization */
    for (i = 0; i < pn; i++) {
      q[i]=p[i]/sum*mn;       /* normalize and scale the weight */
      c+=q[i];                /* the cumulative weight */
      if (c>=1.0) {           /* if the cumulative weight over 1 */
	a=(int)(c);           /* integer part of the cumulative weight */
	c-=a;                 /* subract integer part from the cumulative weight */
	k+=a;                 
	for (b=i+1.0;a>0;a--)	  
	  *s++=b;             /* fill the return matrix 'a' times with sample 'b' */
      }
      if (k<mn && c>=r[k]) {  /* if the cumulative larger than a random number...*/
	*s++=i+1.0;           /* .. add sample... */
	c-=1.0;               /* ...and subtract one from the cumulative */
	k+=1;                 
      }
    }
    mxDestroyArray(R);        /* no need for this anymore */
    
  }

  return;
}
