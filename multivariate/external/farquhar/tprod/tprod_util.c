/*

tprod.C

Generalised multiply accumulate operation. == tensor product with repeated
indicies

TODO:
  1) Would be nice to fall back on BLAS calls for the innermost cases when we
  have a standard 2d x 2d matrix multiply.

$Id: tprod_util.c,v 1.1.1.1 2008/02/27 14:42:55 roboos Exp $

 Copyright 2006-     by Jason Farquhar (jdrf@zepler.org)
 Permission is granted for anyone to copy, use, or modify this
 software and accompanying documents for any uncommercial purposes,
 provided this copyright notice is retained, and note is made of
 any changes that have been made. This software and documents are
 distributed without any warranty, express or implied.

 This code is based upon ideas by Douglas M. Schwarz (schwarz@servtech.com)
 and Adi Vehtari -- but heavily modifed to make it more general, easier to
 understand and faster!


*/

#include "mxInfo.h"
#include "tprod.h"

/*-------------------------------------------------------------------------*/
/* initialise the results matrix */
MxInfo initzmxInfo(int znd, 
						 const MxInfo xinfo, const MxInfo yinfo,
						 const int *x2yIdx, int xnidx, int ynidx){
  int i;
  MxInfo zinfo      = mkemptymxInfo(znd);
  zinfo.numel=1;
  zinfo.stride[0]=1;
  /* first ensure all output dims start with unit size -- so if we specify to
	  map input dim beyond the end we pad the gaps out with unit dims */
  for (i=0;i<znd;i++) zinfo.sz[i]=1;
  /* Now fill in the sizes specified in x2yIdx */
  for( i=0; i<xnidx+ynidx; i++){ /* loop over x2yIdx */
	 if ( x2yIdx[i]>0 ) { /* output dim given */
		if ( i<xnidx ) zinfo.sz[x2yIdx[i]-1]=sz(xinfo,i);
		else           zinfo.sz[x2yIdx[i]-1]=sz(yinfo,i-xnidx);
	 }
  }
  /* then compute the strides/numel & fill in the blanks*/
  for( i=0; i<zinfo.nd; i++){
	 /*if ( zinfo.sz[i]==0 ) { zinfo.sz[i]=1; }*/ /* Why did I do this? */
	 zinfo.stride[i+1] = zinfo.sz[i]*zinfo.stride[i];
	 zinfo.numel      *= zinfo.sz[i];
  }

  /* set the result type to the least precise of the input types */
  if ( xinfo.dtype == DOUBLE_DTYPE && yinfo.dtype == DOUBLE_DTYPE ) {
	 zinfo.dtype    = DOUBLE_DTYPE ;
  } else if ( xinfo.dtype == SINGLE_DTYPE || yinfo.dtype == SINGLE_DTYPE ){
	 zinfo.dtype    = SINGLE_DTYPE;
  } else { /* set to invalid type */
	 zinfo.dtype = 0 ; 
  }

  return zinfo;
}

/*-------------------------------------------------------------------------*/
/* initialize the size of the accumulated sub-matrices and update the x/y
	matrices */
void initrestmaccmxInfo(int znd,
							  const MxInfo xinfo, const MxInfo yinfo, 
							  const int x2yIdx[], int xnidx, int ynidx,
							  MxInfo *xrestinfo,  MxInfo *yrestinfo,
							  MxInfo *xmaccinfo,  MxInfo *ymaccinfo){
  int maccnd=0;
  int i;
  if( xrestinfo->nd < znd || yrestinfo->nd < znd ) 
	 mexErrMsgTxt("mxInfo structures size too small");

  /* ensure all output dims start with unit size -- so if we specify to map
	  input dim beyond the end we pad the gaps out with unit dims */
  for (i=0;i<xrestinfo->nd;i++) xrestinfo->sz[i]=1;
  for (i=0;i<yrestinfo->nd;i++) yrestinfo->sz[i]=1;

  /* ensure outputs have the same types as the inputs */
  xrestinfo->dtype = xinfo.dtype; xmaccinfo->dtype = xinfo.dtype;
  yrestinfo->dtype = yinfo.dtype; ymaccinfo->dtype = yinfo.dtype;

  /* this works in 2 stages -- first it sets up ?maccinfo/?restinfo and then
     squeezes the restinfo to minimum size */
  /* do the x part of the input */
  for( i=0; i<xnidx; i++ ) {

	 if ( x2yIdx[i] < 0 ) { /* accumulated dimension */
		/* compute the accumulated sub-matrices sizes */
		/* in the accumulated matrix the x size is the dim size */
		xmaccinfo->stride[maccnd]= stride(xinfo,i);
		xmaccinfo->sz[maccnd]    = sz(xinfo,i);
		/* in the accumulated matrix the y size is the dim size */
		ymaccinfo->stride[maccnd]= stride(yinfo,-x2yIdx[i]-1);
		ymaccinfo->sz[maccnd]    = sz(yinfo,-x2yIdx[i]-1);
		
		if ( sz(*xmaccinfo,maccnd) != sz(*ymaccinfo,maccnd) ) { 
		  mexErrMsgTxt("Accumulated dimensions sizes must match!");		  
		}
		maccnd++; /* 1 more accumulated dimension */
		
	 } else { /* non-accumulated dimension -- just copy to rest */
		/* in the x rest set, y's size is 1 */
		xrestinfo->stride[x2yIdx[i]-1] = stride(xinfo,i); 
		xrestinfo->sz[x2yIdx[i]-1]     = sz(xinfo,i);		
		yrestinfo->stride[x2yIdx[i]-1] = 0; 
		yrestinfo->sz[x2yIdx[i]-1]     = 1;
	 }
  }
  xmaccinfo->nd=maccnd; ymaccinfo->nd=maccnd; /* set actual number macc dims*/

  for( ; i<xnidx+ynidx; i++){ /* do the y part of x2yIdx */
	 if ( x2yIdx[i] < 0 ) { /* acc dimension -- allready done */
	
	 } else if ( x2yIdx[i] > 0){/* non-acc dimension -- just copy to rest */
 		if ( xrestinfo->sz[x2yIdx[i]-1] == 0 ){ /* if not allready set */
		  xrestinfo->stride[x2yIdx[i]-1] = 0; 
		  xrestinfo->sz[x2yIdx[i]-1]     = 1;
		}
		yrestinfo->stride[x2yIdx[i]-1] = stride(yinfo,i-xnidx); 
		yrestinfo->sz[x2yIdx[i]-1]     = sz(yinfo,i-xnidx);		
	 }
  }

  /* deal with size edge cases */
  if ( xrestinfo->nd==0 )    xrestinfo->nd=1;    /* everything macc'd away */
  if ( xrestinfo->sz[0]==0 ) xrestinfo->sz[0]=1; /* everything macc'd away */
  if ( xrestinfo->sz[1]==0 ) xrestinfo->sz[1]=1; 
  if ( xmaccinfo->nd == 0 )  xmaccinfo->sz[0]=1; /* nothing macc'd away */
  if ( ymaccinfo->nd == 0 )  ymaccinfo->sz[0]=1; /* nothing macc'd away */

  /* set the numels & the final strides -- to simplify debuggin */
  xrestinfo->numel=xrestinfo->sz[0];
  for(i=1;i<xrestinfo->nd;i++)xrestinfo->numel*=xrestinfo->sz[i];
  xrestinfo->stride[xrestinfo->nd]=xinfo.stride[xinfo.nd];
  yrestinfo->numel=yrestinfo->sz[0];
  for(i=1;i<yrestinfo->nd;i++)yrestinfo->numel*=yrestinfo->sz[i];
  yrestinfo->stride[yrestinfo->nd]=yinfo.stride[yinfo.nd];
  xmaccinfo->numel=xmaccinfo->sz[0];
  for(i=1;i<xmaccinfo->nd;i++)xmaccinfo->numel*=xmaccinfo->sz[i];
  xmaccinfo->stride[xmaccinfo->nd]=xinfo.stride[xinfo.nd];
  ymaccinfo->numel=ymaccinfo->sz[0];
  for(i=1;i<ymaccinfo->nd;i++)ymaccinfo->numel*=ymaccinfo->sz[i];
  ymaccinfo->stride[ymaccinfo->nd]=yinfo.stride[yinfo.nd];

  /* set the data pointers */
  xrestinfo->rp=xinfo.rp; xrestinfo->ip=xinfo.ip;
  yrestinfo->rp=yinfo.rp; yrestinfo->ip=yinfo.ip;
  xmaccinfo->rp=xinfo.rp; xmaccinfo->ip=xinfo.ip;
  ymaccinfo->rp=yinfo.rp; ymaccinfo->ip=yinfo.ip;
  
}

/* take as input a pair of mxInfo sets and squeeze them to remove redundant
	extra dimensions which either:
	a) are contiguous in both inputs
	b) are singleton's in both inputs
	c) are consistently contiguous in 1 and stride=0 in the other
*/
void squeezemxInfoPair(MxInfo *xinf, MxInfo *yinf){
  int i=0,nd=0;
  for(; i<xinf->nd; i++){
	 if (xinf->sz[i]<=1 && yinf->sz[i]<=1){
		/* both unset or unit size, ignore */
		
	 } else if( nd>0 && /* normal,and both prev sizes are contiguous */
					xinf->stride[i]==xinf->stride[nd-1]*xinf->sz[nd-1] &&
					yinf->stride[i]==yinf->stride[nd-1]*yinf->sz[nd-1] ) {
		/* merge this dim into the previous one */
		nd--;
		xinf->sz[nd]     = xinf->sz[nd]*xinf->sz[i];
		yinf->sz[nd]     = yinf->sz[nd]*yinf->sz[i];
		nd++;

	 } else { /* normal and non-contiguous */
		xinf->sz[nd]    = xinf->sz[i];	xinf->stride[nd]= xinf->stride[i];
		yinf->sz[nd]    = yinf->sz[i];	yinf->stride[nd]= yinf->stride[i];
		nd++;
	 }
  }
  if ( nd==0 ){ /* deal with all macc/op'd edge cases */
	 nd=1; /* make us look like a pair of scalars */
	 xinf->sz[0]=1; xinf->stride[0]=1; xinf->stride[1]=1;
	 yinf->sz[0]=1; yinf->stride[0]=1; yinf->stride[1]=1;
  } else {
	 xinf->stride[nd]=xinf->stride[i]; /* set the final strides */
	 yinf->stride[nd]=yinf->stride[i];
  }
  xinf->nd=nd; yinf->nd=nd;         /* set the new size */
}

/* optimise the order of dimension in the query so the tprod code works most
	efficiently */
void optimisetprodQuery(MxInfo *zrest, MxInfo *xrest, MxInfo *yrest, 
								MxInfo *xmacc, MxInfo *ymacc){
  int i;
  /* first just squeeze out redundant dimesions */
  squeezemxInfoPair(xrest,yrest);
  squeezemxInfoPair(xmacc,ymacc);

  /* validate that the input is of the right type */
  if ( zrest->nd == xrest->nd ){ 
	 /* then xrestinfo can't have been squashed so zrest is just zinfo */
	 
  } else { /* xrest is squashed, so compute the squashed zrest if possible */
	 if ( !isContiguous(*zrest) ) 
		mexErrMsgTxt("Z must be contiguous for tprod to work!");

	 /* we need to compute zrestinfo now to know how to stride over z */
	 zrest->numel=1;
	 zrest->stride[0]=zrest->stride[0]; /* should be 1? */
	 for (i=0;i<xrest->nd;i++){ /* N.B. this assumes z is contiguous! */
		zrest->sz[i]       = MAX(xrest->sz[i],yrest->sz[i]);
		zrest->stride[i+1] = zrest->sz[i]*zrest->stride[i];
		zrest->numel      *= zrest->sz[i];
	 }
	 zrest->nd=xrest->nd;
  }
  /*  if ( zrest->nd==0 ) zrest->nd=1; */ 


  /*move any aligned dimensions out of the first 2 dims so we can use the 2x2
	 step code*/
/*   for(i=0; i< MAX(zrest->nd,2); i++){ */
/* 	 if( xrest->sz[i]==yrest->sz[i] && */
/* 		  xrest->stride[i]>0 && yrest->stride[i]>0 ){ /\* aligned dimension! *\/ */
/* 		/\* see if we can find an op dim to swap with? *\/ */
/* 		int j,k; */
/* 		for (j=i+1;j<zrest->nd;j++){ */
/* 		  if ( xrest->sz[j]!=yrest->sz[j]  /\* not aligned dim *\/ */
/* 				 || xrest->stride[j]==0 || yrest->stride[j]==0 ) break; */
/* 		} */
/* 		/\* insert this op dim in this place and move everything else up by 1 *\/ */
/* 		if ( j < zrest->nd ) { */
/* 		  int zstride=zrest->stride[j], zsz=zrest->sz[j]; */
/* 		  int xstride=xrest->stride[j], xsz=xrest->sz[j]; */
/* 		  int ystride=yrest->stride[j], ysz=yrest->sz[j]; */
/* 		  for(k=j;k>i;k--){ /\* move down by 1 *\/ */
/* 			 zrest->stride[k]=zrest->stride[k-1]; zrest->sz[k]=zrest->sz[k-1]; */
/* 			 xrest->stride[k]=xrest->stride[k-1]; xrest->sz[k]=xrest->sz[k-1]; */
/* 			 yrest->stride[k]=yrest->stride[k-1]; yrest->sz[k]=yrest->sz[k-1]; */
/* 		  } */
/* 		  zrest->stride[k]=zstride; zrest->sz[k]=zsz;/\* insert in new place *\/ */
/* 		  xrest->stride[k]=xstride; xrest->sz[k]=xsz; */
/* 		  yrest->stride[k]=ystride; yrest->sz[k]=ysz; */
									 
/* 		} else { /\* no OP dims so nothing we can do ! *\/ */
/* 		} */
/* 	 } */
/*   } */

  /* N.B. In the inner loops we ASSUME that X comes first in the output.
	  Hence we re-arrange X and Y to ensure this is so */
  if( xrest->stride[0]==0 || 
		( yrest->stride[0]>0 && xrest->stride[0]>yrest->stride[0] )){
	 int tmp;
	 int *tmpip;
	 double *tmpdp;
	 /* swap the rest  */
	 tmp  =yrest->nd;     yrest->nd = xrest->nd;       xrest->nd = tmp;
	 tmp  =yrest->numel;  yrest->numel=xrest->numel;   xrest->numel = tmp;
	 tmpip=yrest->stride; yrest->stride=xrest->stride; xrest->stride=tmpip;
	 tmpip=yrest->sz;     yrest->sz=xrest->sz;         xrest->sz=tmpip;
	 tmpdp=yrest->rp;     yrest->rp=xrest->rp;         xrest->rp=tmpdp;
	 tmpdp=yrest->ip;     yrest->ip=xrest->ip;         xrest->ip=tmpdp;
	 tmp  =yrest->dtype;  yrest->dtype=xrest->dtype;   xrest->dtype=tmp;

	 /* swap the macc  */
	 tmp  =ymacc->nd;     ymacc->nd = xmacc->nd;       xmacc->nd = tmp;
	 tmp  =ymacc->numel;  ymacc->numel=xmacc->numel;   xmacc->numel = tmp;
	 tmpip=ymacc->stride; ymacc->stride=xmacc->stride; xmacc->stride=tmpip;
	 tmpip=ymacc->sz;     ymacc->sz=xmacc->sz;         xmacc->sz=tmpip;
	 tmpdp=ymacc->rp;     ymacc->rp=xmacc->rp;         xmacc->rp=tmpdp;
	 tmpdp=ymacc->ip;     ymacc->ip=xmacc->ip;         xmacc->ip=tmpdp;
	 tmp  =ymacc->dtype;  ymacc->dtype=xmacc->dtype;   xmacc->dtype=tmp;
  }
}
