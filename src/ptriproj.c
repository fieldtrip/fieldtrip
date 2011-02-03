#include <math.h>
#include "mex.h"
#include "matrix.h"
#include "geometry.h"

void
mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
  char str[256];
  mxArray *proj, *dist;
  double *v1, *v2, *v3, *r, *proj_p, *dist_p;
  int flag;

  if (nrhs <4 || nrhs>5)
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

  if (nrhs==5)
    flag = (int)mxGetScalar (prhs[4]);
  else
    flag = 0;

  proj = mxCreateDoubleMatrix (1, 3, mxREAL);
  dist = mxCreateDoubleMatrix (1, 1, mxREAL);
  proj_p = mxGetData (proj);
  dist_p = mxGetData (dist);

  (*dist_p) = ptriproj(v1, v2, v3, r, proj_p, flag);

  /* assign the output parameters */
  plhs[0] = proj;
  if (nlhs>1)
    plhs[1] = dist;

  return;
}

