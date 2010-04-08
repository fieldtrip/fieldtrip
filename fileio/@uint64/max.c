/*
 *
 * Copyright (C) 2007, Robert Oostenveld
 *
 * This implements [y, i] = max(x)
 *
 * This function does not check the occurence of integer overflowing.
 *
 */

#include <math.h>
#include "mex.h"
#include "matrix.h"

void
mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
  size_t M, N, n;
  UINT64_T *x, *m;
  double *i;
  
  if (nlhs > 2)
    mexErrMsgTxt ("Invalid number of output arguments");
  if (nrhs > 1)
    mexErrMsgTxt ("Invalid number of input arguments");
  
  if (!mxIsUint64(prhs[0]))
    mexErrMsgTxt ("Invalid type of input arguments (should be uint64)");
  
  M = mxGetM(prhs[0]);
  N = mxGetN(prhs[0]);
  
  if (N>1 & M>1)
    mexErrMsgTxt ("This works only on 1-D arrays (i.e. vectors)");
  
  if (N==0 & M==0)
    mexErrMsgTxt ("This works only on 1-D arrays (i.e. vectors)");
  
  
  plhs[0] = (mxArray*)mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
  plhs[1] = (mxArray*)mxCreateDoubleMatrix(1, 1, mxREAL);
  x = mxGetData(prhs[0]);
  m = mxGetData(plhs[0]);
  i = mxGetData(plhs[1]);
  
  n = 0;
  m[0] = x[n];
  i[0] = n+1; /* offset by one */
  
  for (n=1; n<(M*N); n++)
    if (x[n]>m[0])
  {
    m[0] = x[n];
    i[0] = n+1; /* offset by one */
    }
  
  return;
}

