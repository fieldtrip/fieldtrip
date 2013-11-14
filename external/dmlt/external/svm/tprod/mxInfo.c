/*

Main code for wrapping matlab's matrices in something else.

$Id$



Copyright 2006-     by Jason D.R. Farquhar (jdrf@zepler.org)
Permission is granted for anyone to copy, use, or modify this
software and accompanying documents for any uncommercial
purposes, provided this copyright notice is retained, and note is
made of any changes that have been made. This software and
documents are distributed without any warranty, express or
implied


*/

#include "mxInfo.h"

/*-------------------------------------------------------------------------*/
/* struct to hold useful info for iterating over a n-d matrix              */
/* e.g. for 3 x 3 x 3 matrix:
   ndim=3, numel=27, sz=[2 2 2], stride=[1 3 9] */

MxInfo mkmxInfo(int nd, const int *dim, double *rp, double *ip, MxInfoDTypes dtype){
  MxInfo minfo;
  int i;
  if( !(dtype == DOUBLE_DTYPE || dtype==SINGLE_DTYPE || dtype==INT32_DTYPE) ){
	 ERROR("Only defined for double/single/int32 matrices -- sorry!");
  }
  minfo.dtype = dtype;
  if ( nd==0 ) nd=1;   /* always at least 1 dim */
  minfo.nd = nd;
  /* use a single malloc call to get all the mem we need */
  minfo.sz     = (int *)MALLOC(sizeof(int)*(minfo.nd*2+1));
  minfo.stride = minfo.sz+minfo.nd;
  minfo.numel     = 1;
  minfo.stride[0] = 1; 
  for (i = 0; i < minfo.nd; i++) {  /* copy dims to get result size */
	 minfo.sz[i]= (i<nd) ? dim[i] : 1;/* comp the max idx, 1 pad extras */
	 /* compute the x/y strides for this dim */
	 minfo.stride[i+1] = minfo.stride[i]*minfo.sz[i]; /* nd+1 strides */
	 minfo.numel  *= minfo.sz[i]; /* total matrix size */
  }
  minfo.rp = rp;
  minfo.ip = ip; 
  return minfo;
}

void printMxInfoSummary(FILE* fd, const MxInfo info){
  int i;
  fprintf(fd,"[");
  if ( info.nd!=1 ){
  for ( i=0; i < info.nd-1; i++ ) { fprintf(fd,"%dx",info.sz[i]); }
  fprintf(fd,"%d",info.sz[i]);
  } else {
	 fprintf(fd,"%dx1",info.sz[0]);
  }
  fprintf(fd,"]");
  if ( info.dtype==SINGLE_DTYPE ) fprintf(fd," (single"); else fprintf(fd," (double");
  if ( info.ip==0 ) fprintf(fd,")"); else fprintf(fd," complex)");
}

void printMxInfoData(FILE* fd, const MxInfo info){
  int i,j,idx;
  int nCol = (info.nd>1?info.sz[1]:1);
  for ( i=0; i < nCol; i++){
	 for ( j=0; j < info.sz[0]; j++){
		idx = i*info.stride[1]+j*info.stride[0];
		if( info.dtype==DOUBLE_DTYPE ) { 
		  fprintf(fd,"%lf ",info.rp[idx]); 
		} else {                         
		  fprintf(fd,"%f ",(((float*)info.rp)[idx]));  
		}
		if ( info.ip!=0 ) {
		  if( info.dtype==DOUBLE_DTYPE ) { 
			 fprintf(fd,"+ i %lf ",info.ip[idx]); 
		  } else {                         
			 fprintf(fd,"+ i %f ",(((float*)info.ip)[idx])); 
		  }
		}
	 }
  }
}

void printMxInfo(FILE* fd, const MxInfo info){
  printMxInfoSummary(fd,info); fprintf(fd,"\n"); /* summary */
  printMxInfoData(fd,info); fprintf(fd,"\n");/* data */
}


MxInfo mkemptymxInfo(int nd){
  MxInfo minfo;
  nd = (nd==0)?1:nd; /* always at least 1 dim */
  minfo.nd = nd;
  minfo.numel=0;
  /* use a single malloc call to get all the mem we need */
  minfo.sz     = (int *)CALLOC(nd*2+1,sizeof(int));
  minfo.stride = minfo.sz+nd;
  minfo.rp = 0;
  minfo.ip = 0; 
  minfo.dtype = DOUBLE_DTYPE;
  return minfo;
}

MxInfo copymxInfo(const MxInfo inf){
  MxInfo minfo;
  int i;
  minfo.nd = inf.nd;
  minfo.numel=inf.numel;
  /* use a single malloc call to get all the mem we need */
  minfo.sz     = (int *)CALLOC(inf.nd*2+1,sizeof(int));
  minfo.stride = minfo.sz+inf.nd;
  for(i=0;i<inf.nd;i++){ 
	 minfo.sz[i]=inf.sz[i]; 
	 minfo.stride[i]=inf.stride[i];
  }
  minfo.stride[i]=inf.stride[i]; /* nd+1 valid strides */
  minfo.rp = inf.rp;
  minfo.ip = inf.ip; 
  minfo.dtype = inf.dtype;
  return minfo;
}

void initmxInfo(MxInfo *minfo) {
  if ( minfo ) minfo->sz=0;
}

void delmxInfo(MxInfo *minfo) {
  if ( minfo && minfo->sz) FREE(minfo->sz);
}

/* compute if the input info array is contiguous in memory */
char isContiguous(const MxInfo info){
  int d;
  for(d=0; d<info.nd; d++) 
	 if ( info.stride[d+1] != info.stride[d]*info.sz[d] ) break; 
  return info.stride[0]==1 && d==info.nd;
}

/* inlined away ! --- or not because of silly MACS */
int sz(const MxInfo info,int i){
  return (i<info.nd)?info.sz[i]:1;
}
int stride(const MxInfo info, int i){
  return (i<info.nd+1)?info.stride[i]:info.stride[info.nd];
}

int dsz_bytes(const MxInfo info){
  int dsz=0;
  switch( info.dtype ){
  case LOGICAL_DTYPE: dsz=1; break;
  case CHAR_DTYPE   : dsz=1; break;
  case DOUBLE_DTYPE : dsz=sizeof(double); break;
  case SINGLE_DTYPE : dsz=sizeof(float); break;
 /*  case INT8_DTYPE   : dsz=1; */
  case UINT8_DTYPE  : dsz=1; break;
  case INT16_DTYPE  : dsz=2; break;
  case UINT16_DTYPE : dsz=2; break;
  case INT32_DTYPE  : dsz=4; break;
  case UINT32_DTYPE : dsz=4; break;
  case INT64_DTYPE  : dsz=8; break;
  case UINT64_DTYPE : dsz=8; break;
  }
  return dsz;
}

void copyData(const MxInfo info, const double *from, double *to){
  /* copy the real part */
  const char *xp=(char*)from;
  char *zp=(char *)to;
  int dsz=dsz_bytes(info);
  const char *zendp=(char *)zp+(info.numel)*dsz;
  if ( isContiguous(info) ){ /* simple linear copy */
	 while ( zp < zendp ) *zp++ = *xp++;
  } else {
	 int i;
	 int *subs =(int*)CALLOC(info.nd,sizeof(int));
	 while ( zp < zendp ) {
		for( i=0; i < dsz; i++ ) zp[i] = xp[i];
		zp += dsz;
		for( i=0; i < info.nd; i++ ){
		  /* if reached the last element of this dim */
		  xp += info.stride[i]*dsz; 
		  subs[i]++;        /* move on to the next element */
		  if( subs[i] < info.sz[i] ){/*move this dim on by one and stop!*/
			 break;
		  } else {
			 subs[i] = 0; /* reset to the start again! */
			 xp -= info.stride[i]*info.sz[i]*dsz; 
		  } 
		}
	 }
	 FREE(subs);
  }
}


