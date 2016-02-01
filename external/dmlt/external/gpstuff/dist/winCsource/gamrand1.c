/* 
 * GAMRAND1 Random matrices from gamma distribution.
 *
 *  R = GAMRAND1(A,B) returns a matrix of random numbers chosen   
 *  from the gamma distribution with parameters A and B.
 *  The size of R is the common size of A and B if both are matrices.
 *  If either parameter is a scalar, the size of R is the size of the other
 *  parameter. Alternatively, R = GAMRAND(A,B,M,N) returns an M by N matrix. 
 *
 *   Note: Parameterization as in (Neal, 1996).
 *     A is mean of the distribution
 *     B is degrees of freedom
 *
 *
 * Last modified: 2000-05-31 12:16:57 EEST
 *
 */

/* Copyright (C) 1998-2000 Aki Vehtari
 * 
 *This software is distributed under the GNU General Public 
 *License (version 2 or later); please refer to the file 
 *License.txt, included with the software, for details.
 *
 */

#include <math.h>
#include "mex.h"
#include "rand.h"

void mexFunction(const int nlhs_, mxArray *plhs_[],
		 const int nrhs_, const mxArray *prhs_[])
{
   if (nrhs_ != 2 )
      mexErrMsgTxt( "Wrong number of input arguments." );

   if (nlhs_ > 1 )
      mexErrMsgTxt( "Too many output arguments." );

   {
     double a, b;
     a=mxGetScalar(prhs_[0]);
     b=mxGetScalar(prhs_[1]);
     
     plhs_[0]=mxCreateDoubleMatrix(1,1,mxREAL);
     *mxGetPr(plhs_[0]) =
       rand_gamma(b/2)*2*a/b;
   }
  return;
}     

