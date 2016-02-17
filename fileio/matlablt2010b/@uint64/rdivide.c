/*
 *
 * Copyright (C) 2006-2007, Robert Oostenveld
 *
 * This implements C = A / B for uint64 data types
 * and either A and B have the same size or one of them is scalar.
 *
 * This function does not check the occurence of integer overflowing.
 *
 */

#include "mex.h"
#include "matrix.h"
#include <math.h>

#if defined(_WIN32) || defined(_WIN64)
#include <limits.h>
#else
#include <stdint.h>
#endif

#ifndef UINT64_MAX
#define UINT64_MAX  _UI64_MAX
#endif

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
      if (b[i]==0)
    {
      mexWarnMsgTxt("Divide by zero.");
      c[i] = UINT64_MAX;
      }
      else
        c[i] = a[0]/b[i];
  }
  else if (bScalar)
  {
    for (i=0; i<(M*N); i++)
      if (b[0]==0)
    {
      mexWarnMsgTxt("Divide by zero.");
      c[i] = UINT64_MAX;
      }
      else
        c[i] = a[i]/b[0];
  }
  else
  {
    for (i=0; i<(M*N); i++)
      if (b[i]==0)
    {
      mexWarnMsgTxt("Divide by zero.");
      c[i] = UINT64_MAX;
      }
      else
        c[i] = a[i]/b[i];
  }
  
  return;
}

