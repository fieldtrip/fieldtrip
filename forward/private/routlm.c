#include <math.h>
#include "mex.h"
#include "matrix.h"
#include "geometry.h"

void
mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
  mxArray *proj;
  double la, mu;
  double *v1, *v2, *v3, *proj_p;

  if (nrhs != 5)
    mexErrMsgTxt ("Invalid number of input arguments");

  if (mxGetM(prhs[0])!=1 || mxGetN(prhs[0])!=3)
    mexErrMsgTxt ("Invalid dimension for input argument 1");
  if (mxGetM(prhs[1])!=1 || mxGetN(prhs[1])!=3)
    mexErrMsgTxt ("Invalid dimension for input argument 2");
  if (mxGetM(prhs[2])!=1 || mxGetN(prhs[2])!=3)
    mexErrMsgTxt ("Invalid dimension for input argument 3");

  v1 = mxGetData (prhs[0]);
  v2 = mxGetData (prhs[1]);
  v3 = mxGetData (prhs[2]);
  la = mxGetScalar (prhs[3]);
  mu = mxGetScalar (prhs[4]);

  proj = mxCreateDoubleMatrix (1, 3, mxREAL);
  proj_p = mxGetData (proj);

  routlm(v1, v2, v3, la, mu, proj_p);

  /* assign the output parameters */
  plhs[0] = proj;

  return;
}

