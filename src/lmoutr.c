#include <math.h>
#include "mex.h"
#include "geometry.h"

void
mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
  mxArray *la, *mu, *ze;
  double *la_p, *mu_p, *ze_p;
  double *v1, *v2, *v3, *r;

  if (nrhs != 4)
    mexErrMsgTxt ("Invalid number of input arguments");

  if (mxGetM(prhs[0])!=1 || mxGetN(prhs[0])!=3)
    mexErrMsgTxt ("Invalid dimension for input argument 1");
  if (mxGetM(prhs[1])!=1 || mxGetN(prhs[1])!=3)
    mexErrMsgTxt ("Invalid dimension for input argument 2");
  if (mxGetM(prhs[2])!=1 || mxGetN(prhs[2])!=3)
    mexErrMsgTxt ("Invalid dimension for input argument 3");
  if (mxGetM(prhs[3])!=1 || mxGetN(prhs[3])!=3)
    mexErrMsgTxt ("Invalid dimension for input argument 4");

  v1 = mxGetData (prhs[0]);
  v2 = mxGetData (prhs[1]);
  v3 = mxGetData (prhs[2]);
  r  = mxGetData (prhs[3]);

  la = mxCreateDoubleMatrix (1, 1, mxREAL);
  mu = mxCreateDoubleMatrix (1, 1, mxREAL);
  ze = mxCreateDoubleMatrix (1, 1, mxREAL);
  la_p = mxGetData (la);
  mu_p = mxGetData (mu);
  ze_p = mxGetData (ze);

  lmoutr(v1, v2, v3, r, la_p, mu_p, ze_p);

  /* assign the output parameters */
  plhs[0] = la;
  if (nlhs>1)
    plhs[1] = mu;
  if (nlhs>2)
    plhs[2] = ze;

  return;
}

