#include <math.h>
#include "mex.h"
#include "geometry.h"

void
mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
  char str[256];
  mxArray *sect;
  double *v1, *v2, *v3, *l1, *l2, *sect_p;

  if (nrhs!=5)
    mexErrMsgTxt ("Invalid number of input arguments");

  if (mxGetM(prhs[0])!=1 || mxGetN(prhs[0])!=3)
    mexErrMsgTxt ("Invalid dimension for input argument 1");
  if (mxGetM(prhs[1])!=1 || mxGetN(prhs[1])!=3)
    mexErrMsgTxt ("Invalid dimension for input argument 2");
  if (mxGetM(prhs[2])!=1 || mxGetN(prhs[2])!=3)
    mexErrMsgTxt ("Invalid dimension for input argument 3");
  if (mxGetM(prhs[3])!=1 || mxGetN(prhs[3])!=3)
    mexErrMsgTxt ("Invalid dimension for input argument 4");
  if (mxGetM(prhs[4])!=1 || mxGetN(prhs[4])!=3)
    mexErrMsgTxt ("Invalid dimension for input argument 5");

  v1 = mxGetData (prhs[0]);
  v2 = mxGetData (prhs[1]);
  v3 = mxGetData (prhs[2]);
  l1 = mxGetData (prhs[3]);
  l2 = mxGetData (prhs[4]);

  sect = mxCreateDoubleMatrix (1, 3, mxREAL);
  sect_p = mxGetData (sect);

  ltrisect(v1, v2, v3, l1, l2, sect_p);

  /* assign the output parameters */
  plhs[0] = sect;

  return;
}

