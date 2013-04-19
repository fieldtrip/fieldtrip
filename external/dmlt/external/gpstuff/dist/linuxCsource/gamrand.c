/* 
 * GAMRAND   Random matrices from gamma distribution.
 *
 *  R = GAMRAND(A,B) returns a matrix of random numbers chosen   
 *  from the gamma distribution with parameters A and B.
 *  The size of R is the common size of A and B if both are matrices.
 *  If either parameter is a scalar, the size of R is the size of the other
 *  parameter. Alternatively, R = GAMRAND(A,B,M,N) returns an M by N matrix. 
 *
 *   Note: Parameterization as in (Neal, 1996).
 *     A is mean of the distribution
 *     B is degrees of freedom
 *
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
  
  if ((nrhs < 2) | (nrhs == 3) | (nrhs > 4))
    mexErrMsgTxt( "Wrong number of input arguments." );

   {
     double a, b, *aa, *bb, *r;
     const int *dims1, *dims2;
     int i, m, n, mn;
    
     dims1 = mxGetDimensions(prhs[0]);
     dims2 = mxGetDimensions(prhs[1]);
     if ((dims1[0] > 1 || dims1[1] > 1) &&
	 (dims2[0] > 1 || dims2[1] > 1) &&
	 (dims1[0]!=dims2[0] || dims1[1]!=dims2[1]))
       mexErrMsgTxt( "Size information is inconsistent." );
     if (nrhs < 3) {
       if (dims1[0] > 1 || dims1[1] > 1) {
	 m = dims1[0];
	 n = dims1[1];
       }
       else {
	 m = dims2[0];
	 n = dims2[1];
       }
     }
     else {
       m=(int)mxGetScalar(prhs[2]);
       n=(int)mxGetScalar(prhs[3]);
       if (((dims1[0] > 1 || dims1[1] > 1) && dims1[0]!=m && dims1[1]!=n) ||
	   ((dims2[0] > 1 || dims2[1] > 1) && dims2[0]!=m && dims2[1]!=n))
	 mexErrMsgTxt( "Size information is inconsistent." );
     }
     plhs[0]=mxCreateDoubleMatrix(m,n,mxREAL);
     r = mxGetPr(plhs[0]);
     mn=m*n;
    
     if (dims1[0] > 1 || dims1[1] > 1) {
       if (dims2[0] > 1 || dims2[1] > 1) {
	 aa=mxGetPr(prhs[0]);
	 bb=mxGetPr(prhs[1]);
	   for (i = 0; i < mn; i++) 
	     r[i] = rand_gamma(bb[i]/2)*2*aa[i]/bb[i];
       }
       else {
	 aa=mxGetPr(prhs[0]);
	 b=mxGetScalar(prhs[1]);
	 for (i = 0; i < mn; i++) 
	   r[i] = rand_gamma(b/2)*2*aa[i]/b;
       }
     }
     else {
       if (dims2[0] > 1 || dims2[1] > 1) {
	 a=mxGetScalar(prhs[0]);
	 bb=mxGetPr(prhs[1]);
	 for (i = 0; i < mn; i++) 
	   r[i] = rand_gamma(bb[i]/2)*2*a/bb[i];
       }
       else {
	 a=mxGetScalar(prhs[0]);
	 b=mxGetScalar(prhs[1]);
	 for (i = 0; i < mn; i++) 
	   r[i] = rand_gamma(b/2)*2*a/b;
       }
     }
   }
  return;
}     

