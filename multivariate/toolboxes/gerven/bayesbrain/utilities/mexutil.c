#include "mexutil.h"

/* Functions to create uninitialized arrays. */

mxArray *mxCreateNumericArrayE(mwSize ndim, const mwSize *dims, 
         mxClassID classid, mxComplexity ComplexFlag)
{
  mxArray *a;
  mwSize i;
  mwSize *dims1 = mxMalloc(ndim*sizeof(mwSize));
  mwSize sz = 1;
  for(i=0;i<ndim;i++) {
    sz *= dims[i];
    dims1[i] = 1;
  }
  a = mxCreateNumericArray(ndim,dims1,classid,ComplexFlag);
  sz *= mxGetElementSize(a);
  mxSetDimensions(a, dims, ndim);
  mxFree(dims1);
  mxSetData(a, mxRealloc(mxGetData(a), sz));
  if(ComplexFlag == mxCOMPLEX) {
    mxSetPi(a, mxRealloc(mxGetPi(a),sz));
  }
  return a;
}
mxArray *mxCreateNumericMatrixE(mwSize m, mwSize n, mxClassID classid, 
				mxComplexity ComplexFlag)
{
  mwSize sz = m*n*sizeof(double);
  mxArray *a = mxCreateNumericMatrix(1, 1, classid, ComplexFlag);
  mxSetM(a,m);
  mxSetN(a,n);
  mxSetPr(a, mxRealloc(mxGetPr(a),sz));
  if(ComplexFlag == mxCOMPLEX) {
    mxSetPi(a, mxRealloc(mxGetPi(a),sz));
  }
  return a;
}
mxArray *mxCreateDoubleMatrixE(mwSize m, mwSize n, 
			       mxComplexity ComplexFlag)
{
  return mxCreateNumericMatrixE(m,n,mxDOUBLE_CLASS,ComplexFlag);
}

