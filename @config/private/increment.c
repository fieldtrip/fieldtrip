/*
 *
 * Copyright (C) 2007, Robert Oostenveld
 *
 * This implements X+=1, i.e. it increments the scalar X by one.
 *
 */

#include <math.h>
#include "mex.h"

void
mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
  double *x;
  
  if (nlhs > 0)
    mexErrMsgTxt ("Invalid number of output arguments");
  if (nrhs > 1)
    mexErrMsgTxt ("Invalid number of input arguments");
  
  if (!mxIsDouble(prhs[0]))
    mexErrMsgTxt ("Invalid type of input arguments (should be double)");
  
  x = mxGetData(prhs[0]);
  (*x) += 1;
  
  return;
}

