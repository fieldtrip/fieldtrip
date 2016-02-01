/* ADAPTIVE REJECTION SAMPLING FROM CONDITIONAL DISTRIBUTION FOR A SIGMA VALUE.
 * Draws a random value from the conditional distribution for a sigma that is 
 * defined by its top-down prior and by the mean of the lower-level precision 
 * values that it controls, using the Adaptive Rejection Sampling method. 
 *
 *          Description
 *          R = COND_INVGAM_INVGAM(A, A1, A2, X) generates one sample
 *          from the conditional distribution of A given
 *          parameter structure X of lower level, structure A1 of
 *          same level hyper-parameters and A2 of higher level, i.e 
 *          is r~P(A|A1,A2,X). Returns one new sample R from the 
 *          distribution above
 */

/* Copyright (c) 1998-2004 Aki Vehtari  
 *
 *This software is distributed under the GNU General Public 
 *License (version 3 or later); please refer to the file 
 *License.txt, included with the software, for details.
 *
 */

#include <math.h>
#include "mex.h"

#include "rand.h"
#include "ars.h"

typedef struct { double w, a, a0, a1, s; } logp_data;
static double logp (double l, double *d, void *vp);

void mexFunction(const int nlhs_, mxArray *plhs_[],
		 const int nrhs_, const mxArray *prhs_[])
{
   if (nrhs_ != 5 )
      mexErrMsgTxt( "Wrong number of input arguments." );

   if (nlhs_ > 1 )
      mexErrMsgTxt( "Too many output arguments." );

   {
     double n;
     logp_data data;
     
     data.w=1./mxGetScalar(prhs_[0]);
     data.a0=mxGetScalar(prhs_[1]);
     data.a1=mxGetScalar(prhs_[2]);
     n=mxGetScalar(prhs_[4]);
     data.s=1./mxGetScalar(prhs_[3])*n;
     data.a=data.a0-n*data.a1;

     plhs_[0]=mxCreateDoubleMatrix(1,1,mxREAL);
     *mxGetPr(plhs_[0])= 
        exp(-ars(log(data.w),log(1+1/sqrt(data.a0)),logp,&data));  

   }
  return;
}     

/* Copyright (c) 1996 by Carl Edward Rasmussen and Radford M. Neal
 *
 * Permission is granted for anyone to copy, use, or modify this program 
 * for purposes of research or education, provided this copyright notice 
 * is retained, and note is made of any changes that have been made. 
 *
 * This program is distributed without any warranty, express or implied.
 * As this program was written for research purposes only, it has not been
 * tested to the degree that would be advisable in any important application.
 * All use of this program is entirely at the user's own risk.
 */

static double logp (double l, double *d, void *vp)
{ logp_data *p = vp;
  double t = exp(l); 
  double v;
  *d = p->a/2 - t*p->a0/(2*p->w) + p->a1*p->s/(2*t);
  v = l*p->a/2 - t*p->a0/(2*p->w) - p->a1*p->s/(2*t);
  /* fprintf(stderr,"%.3f %g %g\n",t,v,*d); */
  return v;
}
/*
static double logp (double l, double *d, void *vp)
{
  logp_data *p = vp;
  double t = exp(l); 
  double v;
  *d = p->a/2 -  t*p->a0/(2*p->w) + p->a1*p->s/(2*t);
  v = l*p->a/2 - t*p->a0/(2*p->w) - p->a1*p->s/(2*t);
  return v;
}
*/
