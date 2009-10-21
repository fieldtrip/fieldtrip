/*
 *
 * Copyright (C) 2006-2007, Robert Oostenveld
 *
 * This implements [y] = abs(x)
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
  size_t M, N;
  UINT64_T *x, *y;
  int i;
  
  if (nlhs > 2)
    mexErrMsgTxt ("Invalid number of output arguments");
  if (nrhs > 1)
    mexErrMsgTxt ("Invalid number of input arguments");
  
  if (!mxIsUint64(prhs[0]))
    mexErrMsgTxt ("Invalid type of input arguments (should be uint64)");
  
  M = mxGetM(prhs[0]);
  N = mxGetN(prhs[0]);
  
  plhs[0] = (mxArray*)mxCreateNumericMatrix(M, N, mxUINT64_CLASS, mxREAL);
  x = mxGetData(prhs[0]);
  y = mxGetData(plhs[0]);
  
  for (i=0; i<(M*N); i++)
    y[i] = x[i];  /* the absolute value of an unsigned integer is identical to the original value */
  
  return;
}

