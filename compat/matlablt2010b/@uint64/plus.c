/*
 *
 * Copyright (C) 2006-2007, Robert Oostenveld
 *
 * This implements C = A + B for uint64 data types
 * and either A and B have the same size or one of them is scalar.
 *
 * This function does not check the occurence of integer overflowing.
 *
 * Note that UINT64_T is used here because the Borland C++ compiler (free version 5.5) does not
 * understand "unsigned long long", and the gcc compiler does not understand "unsigned __int64".
 * Matlab solves that with UINT64_T
 *
 * Regardless of the declaration of UINT64_T variables I did not manage to compile this code using
 * the LCC 2.4.1 compiler, which gives the following error
 *    >> mex plus.c
 *    Error plus.c: 70  compiler error in _kids--Bad rule number 0
 *    C:\PROGRA~1\MATLAB\R2007B\BIN\MEX.PL: Error: Compile of 'plus.c' failed.
 *
 */

#include <math.h>
#include "mex.h"

void
mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
  size_t M, N, i;
  UINT64_T *a, *b, *c;
  bool aScalar, bScalar;
  
  if (nlhs > 1)
    mexErrMsgTxt ("Invalid number of output arguments");
  if (nrhs != 2)
    mexErrMsgTxt ("Invalid number of input arguments");
  
  if (!mxIsUint64(prhs[0]))
    mexErrMsgTxt ("Invalid type of input arguments (should be uint64)");
  if (!mxIsUint64(prhs[1]))
    mexErrMsgTxt ("Invalid type of input arguments (should be uint64)");
  
  aScalar = ((mxGetN(prhs[0])==1) && (mxGetM(prhs[0])==1));
  bScalar = ((mxGetN(prhs[1])==1) && (mxGetM(prhs[1])==1));
  
  if (!aScalar && !bScalar)
  {
    if (mxGetM(prhs[0]) != mxGetM(prhs[1]))
      mexErrMsgTxt ("Invalid size of input arguments");
    if (mxGetN(prhs[0]) != mxGetN(prhs[1]))
      mexErrMsgTxt ("Invalid size of input arguments");
  }
  
  if (aScalar)
  {
    M = mxGetM(prhs[1]);
    N = mxGetN(prhs[1]);
  }
  else
  {
    M = mxGetM(prhs[0]);
    N = mxGetN(prhs[0]);
  }
  
  a = mxGetData(prhs[0]);
  b = mxGetData(prhs[1]);
  
  plhs[0] = (mxArray*)mxCreateNumericMatrix(M, N, mxUINT64_CLASS, mxREAL);
  c = mxGetData(plhs[0]);
  
  if (aScalar)
  {
    for (i=0; i<(M*N); i++)
      c[i] = a[0]+b[i];
  }
  else if (bScalar)
  {
    for (i=0; i<(M*N); i++)
      c[i] = a[i]+b[0];
  }
  else
  {
    for (i=0; i<(M*N); i++)
      c[i] = a[i]+b[i];
  }
  
  return;
}

