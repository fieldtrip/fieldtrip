#ifndef mxInfoMexH
#define mxInfoMexH
/*
  matlab specific helper functions 
*/

#include "mex.h"
#include "matrix.h"
#include "mxInfo.h"

mxArray* mkmxArrayCopy(const MxInfo info);
/* convert back to mxArray */
mxArray* mkmxArray(const MxInfo info);
MxInfo mkmxInfoMxArray(const mxArray *mat, int nd);

#endif
