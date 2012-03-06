#include "mxInfo_mex.h"

/* memory management functions */
void *CALLOC(size_t nmemb, size_t size){  return mxCalloc(nmemb,size); }
void *MALLOC(size_t size) { return mxMalloc(size); }
void FREE(void *ptr) {  mxFree(ptr); }
void ERROR(const char *msg) { mexErrMsgTxt(msg); }
void WARNING(const char *msg) { mexWarnMsgTxt(msg); }

MxInfo mkmxInfoMxArray(const mxArray *mat, int nd){
  MxInfo minfo;
  int i;
  int rnd = mxGetNumberOfDimensions(mat);
  const int *dim = mxGetDimensions(mat);
  int dtype = mxGetClassID(mat);   
  if( !(dtype == DOUBLE_DTYPE || dtype==SINGLE_DTYPE || dtype==INT32_DTYPE) )
	 mexErrMsgTxt("Only defined for double/single/int32 matrices -- sorry!");
  minfo.dtype = dtype;
  if ( nd==0 ) nd=rnd; /* 0 indicates default dim size */
  if ( nd==0 ) nd=1;   /* always at least 1 dim */
  minfo.nd = nd;
  /* use a single malloc call to get all the mem we need */
  minfo.sz     = (int *)MALLOC(sizeof(int)*(minfo.nd*2+1));
  minfo.stride = minfo.sz+minfo.nd;
  minfo.numel     = 1;
  minfo.stride[0] = 1; 
  for (i = 0; i < minfo.nd; i++) {  /* copy dims to get result size */
	 minfo.sz[i]= (i<rnd) ? dim[i] : 1;/* comp the max idx, 1 pad extras */
	 /* compute the x/y strides for this dim */
	 minfo.stride[i+1] = minfo.stride[i]*minfo.sz[i]; /* nd+1 strides */
	 minfo.numel  *= minfo.sz[i]; /* total matrix size */
  }
  minfo.rp = mxGetPr(mat);
  if ( mxIsComplex(mat) ) minfo.ip = mxGetPi(mat); else minfo.ip=0;
  return minfo;
}

mxArray* mkmxArray(const MxInfo info){
  mxArray *mx = 
	 mxCreateNumericMatrix(0,0,info.dtype,(info.ip==0)?mxREAL:mxCOMPLEX);
  /* only valid if strides are consistent with single matrix */
  if ( !isContiguous(info) ) 
	 mexErrMsgTxt("Can't make mxarray from this mxInfo -- no uniform memory");
  
  mxSetDimensions(mx,info.sz,info.nd);
  mxSetPr(mx,info.rp); mxSetPi(mx,info.ip);
  return mx;
  
}

mxArray* mkmxArrayCopy(const MxInfo info){
  mxArray *mx = mxCreateNumericArray(info.nd,info.sz,info.dtype,
												 (info.ip==0)?mxREAL:mxCOMPLEX);
  mxSetDimensions(mx,info.sz,info.nd);
  copyData(info,info.rp,mxGetPr(mx));
  if ( info.ip ) copyData(info,info.ip,mxGetPi(mx));
  return mx;
}
