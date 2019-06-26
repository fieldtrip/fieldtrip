/*
 *
 * Copyright (C) 2007, Robert Oostenveld
 *
 * This implements Y=X by making a deep copy. A deep copy refers
 * to a copy in which all levels of data are copied. For example, a
 * deep copy of a cell-array copies each cell, and the contents of
 * the each cell (if any), and so on.
 *
 */

#include <math.h>
#include "mex.h"

void
mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
  double *x;
  
  if (nlhs > 1)
    mexErrMsgTxt ("Invalid number of output arguments");
  if (nrhs > 1)
    mexErrMsgTxt ("Invalid number of input arguments");
  
  plhs[0] = mxDuplicateArray(prhs[0]);
  
  return;
}

