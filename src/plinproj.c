#include <math.h>
#include "mex.h"
#include "geometry.h"

void
mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
  char str[256];
  mxArray *proj, *dist;
  double *l1, *l2, *r, *proj_p, *dist_p;
  int flag;

  if (nrhs <3 || nrhs>4)
    mexErrMsgTxt ("Invalid number of input arguments");

  if (mxGetM(prhs[0])!=1 || mxGetN(prhs[0])!=3)
    mexErrMsgTxt ("Invalid dimension for input argument 1");
  if (mxGetM(prhs[1])!=1 || mxGetN(prhs[1])!=3)
    mexErrMsgTxt ("Invalid dimension for input argument 2");
  if (mxGetM(prhs[2])!=1 || mxGetN(prhs[2])!=3)
    mexErrMsgTxt ("Invalid dimension for input argument 3");

  l1 = mxGetData (prhs[0]);
  l2 = mxGetData (prhs[1]);
  r  = mxGetData (prhs[2]);

  if (nrhs==4)
    flag = (int)mxGetScalar (prhs[3]);
  else
    flag = 0;

  proj = mxCreateDoubleMatrix (1, 3, mxREAL);
  dist = mxCreateDoubleMatrix (1, 1, mxREAL);
  proj_p = mxGetData (proj);
  dist_p = mxGetData (dist);

  (*dist_p) = plinproj(l1, l2, r, proj_p, flag);

  /* assign the output parameters */
  plhs[0] = proj;
  if (nlhs>1)
    plhs[1] = dist;

  return;
}

