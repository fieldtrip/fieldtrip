/*

tprod.C

Generalised multiply accumulate operation. == tensor product with repeated
indicies

  Z = tprod(X,Y,xdimspec,ydimspec,flags)

where the value in dimspec in determines what type of operation to perform
over that dimension, case are:
  0 -- must be singlento dimension, is squeezed out of the result

Examples:
  X = randn(100,100,100); % dim x samp x trials set of timeseries
  sf= randn(100,1);       % dim x 1 spatial filter
  tf= randn(100,1);       % dim x 1 temporal filter
  sfX=tprod(X,sf,1);      % spatially filter sfX -> [ 1 x samp x trials ]
  tfX=tprod(X,tf,[0 -1]); % temporaly fitler tfX -> [ dim x 1 x trials ]

  % with aligned/sequential dimensions
  Z=tprod(X,X,[1 -2 -3],[1 2 3]); % per trial data variance over dims
  Z=tprod(X,X,[-1 2 -3],[1 2 3]); % per trial data variance over time
  C=tprod(X,X,[1 2 -3],[2 1 3]);  % per trail dim covariance

result has size: [size(X) size(Y)] where the accumulated dims are now size==1

TODO:
  1) Would be nice to fall back on BLAS calls for the innermost cases when we
  have a standard 2d x 2d matrix multiply.

N.B. compile with TPRODLIB defined to get an version without mexFunction
defined to use in other mex-files
  

$Id$

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

/* process the input type information */
#ifndef XDTYPE 
#define XDTYPE DOUBLE_DTYPE
#endif
#ifndef YDTYPE 
#define YDTYPE DOUBLE_DTYPE
#endif

#if (XDTYPE == DOUBLE_DTYPE) && (YDTYPE == DOUBLE_DTYPE)
#define TYPESTR dd
#define XTYPE double
#define YTYPE double
#define ZDTYPE DOUBLE_DTYPE
#define ZTYPE double
#elif (XDTYPE == DOUBLE_DTYPE) && (YDTYPE == SINGLE_DTYPE)
#define TYPESTR ds
#define XTYPE double
#define YTYPE float
#define ZDTYPE SINGLE_DTYPE
#define ZTYPE float
#elif (XDTYPE == SINGLE_DTYPE) && (YDTYPE == DOUBLE_DTYPE)
#define TYPESTR sd
#define XTYPE float
#define YTYPE double
#define ZDTYPE SINGLE_DTYPE
#define ZTYPE float
#elif (XDTYPE == SINGLE_DTYPE) && (YDTYPE == SINGLE_DTYPE)
#define TYPESTR ss
#define XTYPE float
#define YTYPE float
#define ZDTYPE SINGLE_DTYPE
#define ZTYPE float
#elif (XDTYPE== INT32_DTYPE ) && (YDTYPE== INT32_DTYPE)
#define TYPESTR ii
#define XTYPE int
#define YTYPE int
#define ZDTYPE INT32_DTYPE
#define ZTYPE int
#else
#error Unsupported X,Y type combination
#endif

/* utility defines for constructing the function to call */
#define TYPESTRRR CAT(TYPESTR,rxry)
#define TYPESTRRC CAT(TYPESTR,rxcy)
#define TYPESTRCR CAT(TYPESTR,cxry)
#define TYPESTRCC CAT(TYPESTR,cxcy)

/*---------------------------------------------------------------------------*/
/* Use the macro-preprocessor to automatically generate code for each of the 
 * posible combinations of inputs */
/*---------------------------------------------------------------------------*/
#define XISCOMPLEX
#define YISCOMPLEX
/* #define FNNMEXT cxcy */
#include "tprod.def"   /* complex x, complex y */

#undef YISCOMPLEX
#define XISCOMPLEX
/* #define FNNMEXT cxry */
#include "tprod.def"   /* complex x, real y */

#undef XISCOMPLEX
#define YISCOMPLEX
/* #define FNNMEXT rxcy */
#include "tprod.def"   /* real x, complex y */

#undef XISCOMPLEX
#undef YISCOMPLEX
/* #define FNNMEXT rxry */
#include "tprod.def"   /* real x, real y */

/* function to compute the generalised tensor product */
TprodErrorCode CAT(TYPESTR,tprod)(const MxInfo zinfo, 
											 const MxInfo xrestin, const MxInfo yrestin,
											 const MxInfo xmaccin, const MxInfo ymaccin,
											 int blksz){
  /* make a copy of the inputs so we can modify them locally */
  MxInfo zrest=copymxInfo(zinfo);
  MxInfo xrest=copymxInfo(xrestin);
  MxInfo yrest=copymxInfo(yrestin);
  MxInfo xmacc=copymxInfo(xmaccin);
  MxInfo ymacc=copymxInfo(ymaccin);
  
  int intype, retVal = OK ; 
  const int BLKTHRESH=128*128*128; /* thresh size to switch to blk'd code*/
  int i,j;
  int *maccsubs=(int*)CALLOC(xmacc.nd+xrest.nd,sizeof(int));
  int *subs    =maccsubs+xmacc.nd;

  /* check the types match */
  if ( zinfo.dtype != ZDTYPE ) {
	 retVal = ZTYPEMISMATCH;
  } else if ( xrestin.dtype != XDTYPE || xmaccin.dtype != XDTYPE) {
	 retVal = XTYPEMISMATCH;
  } else if ( yrestin.dtype != YDTYPE || ymaccin.dtype != YDTYPE ) {
	 retVal = YTYPEMISMATCH; 
  } 
  if ( retVal != OK ) { /* something didn't match return error code */
	 delmxInfo(&zrest);delmxInfo(&xrest);delmxInfo(&yrest);
	 delmxInfo(&xmacc);delmxInfo(&ymacc);FREE(maccsubs); 
	 return retVal ;
  }


  /* query optimisation -- re-order and squeeze dims for max efficiency */
  /* N.B. we assume this has been done before we're called from now on! */
  /*optimisetprodQuery(&zrest,&xrest,&yrest,&xmacc,&ymacc); */

  /* N.B. opt-query may reverse the inputs so on comp this here */
  /* enum for the input type, rxry=0,cxry=1,rxcy=2,cxcy=3 */
  intype = (xrest.ip==0?0:1) + (yrest.ip==0?0:2); /* now compute type call */
  
  /*-------------------------------------------------------------------------*/
  /* use the fast b22 code when we have op over x/y in first 2 dims */
  /* 2 dims only if 2 output dims and y stride 0 in dim 0 and x stride 0 in
	  dim 1 (and the maccs aren't linear in both x and y) */
  if( zrest.nd>1 && zrest.sz[0]>1 && zrest.sz[1]>1 &&
		yrest.stride[0]==0 && xrest.stride[1]==0 ) {

	 if ( blksz<0 ||  /* use special 2x2x1 code */
			(zrest.nd==2 && xmacc.nd==1 && /* small enough 2x2x1 */
			 zrest.sz[0]*zrest.sz[1]*xmacc.sz[0] < BLKTHRESH) ){ 
		switch (intype){       
		case RXRY: CAT(dgemm22_,TYPESTRRR)(&zrest,&xrest,&yrest,&xmacc,&ymacc); break;
		case RXCY: CAT(dgemm22_,TYPESTRRC)(&zrest,&xrest,&yrest,&xmacc,&ymacc); break;
		case CXRY: CAT(dgemm22_,TYPESTRCR)(&zrest,&xrest,&yrest,&xmacc,&ymacc); break;
		case CXCY: CAT(dgemm22_,TYPESTRCC)(&zrest,&xrest,&yrest,&xmacc,&ymacc); break;
		}

	 } else if (blksz==0 || /* use non-blocked when input small enough */
		  zrest.sz[0]*zrest.sz[1]*xmacc.sz[0] < BLKTHRESH ) { 

		switch (intype){       
		case RXRY: 
		  CAT(b22XY_tprod_,TYPESTRRR)
			 ((ZTYPE*)zrest.rp, (ZTYPE*)zrest.ip, 
			  (XTYPE*)xrest.rp, (XTYPE*)xrest.ip, 
			  (YTYPE*)yrest.rp, (YTYPE*)yrest.ip,
			  &zrest, &xrest, &yrest, &xmacc, &ymacc, subs, maccsubs);
		  break;
		case RXCY: 
		  CAT(b22XY_tprod_,TYPESTRRC)
			 ((ZTYPE*)zrest.rp, (ZTYPE*)zrest.ip, 
			  (XTYPE*)xrest.rp, (XTYPE*)xrest.ip, 
			  (YTYPE*)yrest.rp, (YTYPE*)yrest.ip,
			  &zrest, &xrest, &yrest, &xmacc, &ymacc, subs, maccsubs);
		  break;
		case CXRY: 
		  CAT(b22XY_tprod_,TYPESTRCR)
			 ((ZTYPE*)zrest.rp, (ZTYPE*)zrest.ip, 
			  (XTYPE*)xrest.rp, (XTYPE*)xrest.ip, 
			  (YTYPE*)yrest.rp, (YTYPE*)yrest.ip,
			  &zrest, &xrest, &yrest, &xmacc, &ymacc, subs, maccsubs);
		  break;
		case CXCY: 
		  CAT(b22XY_tprod_,TYPESTRCC)
			 ((ZTYPE*)zrest.rp, (ZTYPE*)zrest.ip, 
			  (XTYPE*)xrest.rp, (XTYPE*)xrest.ip, 
			  (YTYPE*)yrest.rp, (YTYPE*)yrest.ip,
			  &zrest, &xrest, &yrest, &xmacc, &ymacc, subs, maccsubs);
		  break;
		}
 		
	 } else { /* use the blocked code */
        MxInfo zrestIn, xrestIn, yrestIn;
        int maccblksz;
		MxInfo xmaccIn, ymaccIn;   /* inner loop mxInfo description */
		MxInfo xmaccOut, ymaccOut; /* outer macc loop mxInfo description */
		int maccsz, maccblkdim, blkmaccstride, maccOutRem;
        
		/* compute the block size for the 1st 2 dims */
		int restOutrem0, restOutrem1, blksz0, blksz1;
		if ( zrest.sz[0]*zrest.sz[1] < blksz*blksz ) { 
		  /* first 2 dims are less than blksz so use them all in 1 step */
		  /* hence we don't need to step over them at all, so remove them. */
		  blksz0=zrest.sz[0]; blksz1=zrest.sz[1];
		} else { /* limit them to blksz */
		  blksz0=blksz;      blksz1=blksz;
		}
		restOutrem0       =zinfo.sz[0]%blksz0;        /* compute remainder 0 */
		restOutrem1       =zinfo.sz[1]%blksz1;        /* compute remainder 1 */

		zrestIn=mkemptymxInfo(2);xrestIn=mkemptymxInfo(2);yrestIn=mkemptymxInfo(2);

		zrestIn.stride[0]=zrest.stride[0];	zrestIn.sz[0]=blksz0;
		xrestIn.stride[0]=xrest.stride[0];	xrestIn.sz[0]=blksz0;
		yrestIn.stride[0]=0;	               yrestIn.sz[0]=1;		
		zrestIn.stride[1]=zrest.stride[1];	zrestIn.sz[1]=blksz1;
		xrestIn.stride[1]=0;	               xrestIn.sz[1]=1;
		yrestIn.stride[1]=yrest.stride[1];	yrestIn.sz[1]=blksz1;		
	
		/* N.B. this is more complex as we may need to go a few dims down the
			macc to get a sufficiently large block */
		/* modify the macc to work in blocks of blksz */
		/* size to fill cache */
		maccblksz=(3*blksz*blksz-blksz0*blksz1)/(blksz0+blksz1); 

		/*identify the macc dim which gives blksz sum*/
		if ( xmacc.nd==0 ) { /*check for empty macc input*/
		  /* empty macc input */
		  maccsz=0;
		  maccblkdim=0;
		  blkmaccstride=0;
		  maccOutRem=0;
		  xmaccIn=mkemptymxInfo(0);  ymaccIn=mkemptymxInfo(0);
		  xmaccOut=mkemptymxInfo(0);	ymaccOut=mkemptymxInfo(0);

		} else { /* compute the inner and outer macc descriptions */
		  maccsz=xmacc.sz[0];
		  for(i=0;i<xmacc.nd-1&& maccsz<maccblksz;i++,maccsz*=xmacc.sz[i]);
		  if( maccsz < maccblksz ) maccblksz=maccsz;/*take all if macc to small*/
		  maccsz=maccsz/xmacc.sz[i]; /* get num-el up to this dim */
		  maccblkdim=i;                  /* record which dim we blk at */
		  blkmaccstride=maccblksz/maccsz;/*num steps in dim per blksz elm*/
		  /* now we need to cut this dimension evenly */
		  /*attempt to find an int multiple of this stride which matches blksz*/
		  /* 	 for(i=blkmaccstride/2+1;  */
		  /* 		 i<xmacc.sz[maccblkdim]/2 && i<maccblksz*2/maccsz-1; i++){*/
		  /* 		if( xmacc.sz[maccblkdim]%i==0 ) break; */
		  /* 	 } */
		  /* 	 blkmaccstride=i; */

		  /* create an mxInfo for use in the macc inner loop */
		  xmaccIn=mkemptymxInfo(maccblkdim+1);ymaccIn=mkemptymxInfo(maccblkdim+1);
		  for ( i=0; i < maccblkdim; i++){
			 xmaccIn.stride[i]=xmacc.stride[i]; xmaccIn.sz[i]=xmacc.sz[i];
			 ymaccIn.stride[i]=ymacc.stride[i]; ymaccIn.sz[i]=ymacc.sz[i];
		  }
		  xmaccIn.stride[i]=xmacc.stride[i]; xmaccIn.sz[i]=blkmaccstride; 
		  ymaccIn.stride[i]=ymacc.stride[i]; ymaccIn.sz[i]=blkmaccstride;	 
		  maccOutRem = xmacc.sz[i]%blkmaccstride;

		  /* create a mxInfo for the outer macc loop */
		  xmaccOut=mkemptymxInfo(xmacc.nd-maccblkdim);
		  ymaccOut=mkemptymxInfo(ymacc.nd-maccblkdim);
		  if ( blkmaccstride != xmacc.sz[maccblkdim] ){
			 xmaccOut.stride[0]=  xmacc.stride[maccblkdim]*xmaccIn.sz[maccblkdim];
			 xmaccOut.sz[0]    =  xmacc.sz[maccblkdim]    /xmaccIn.sz[maccblkdim];
 			 ymaccOut.stride[0]=  ymacc.stride[maccblkdim]*ymaccIn.sz[maccblkdim];
			 ymaccOut.sz[0]    =  ymacc.sz[maccblkdim]    /ymaccIn.sz[maccblkdim];
			 j=1;
		  } else { /* just copy from 1 dim further on */
			 j=0;
		  }
		  /* copy the sizes up to start of macc info */
		  for(i=maccblkdim+1; i < xmaccOut.nd; i++,j++){
			 xmaccOut.stride[j]=xmacc.stride[i+maccblkdim];
			 xmaccOut.sz[j]    =xmacc.sz[i];
			 ymaccOut.stride[j]=ymacc.stride[i+maccblkdim];
			 ymaccOut.sz[j]    =ymacc.sz[i];
		  }
		  xmaccOut.nd=j; ymaccOut.nd=j;
		} 	 

		switch (intype){       
		case RXRY: 
		  CAT(blk_b22XY_tprod_,TYPESTRRR)
			 ((ZTYPE*)zrest.rp, (ZTYPE*)zrest.ip, 
			  (XTYPE*)xrest.rp, (XTYPE*)xrest.ip, 
			  (YTYPE*)yrest.rp, (YTYPE*)yrest.ip,
			  &zrest, &xrest, &yrest, &xmaccIn,  &ymaccIn,  subs, maccsubs,
			  &xmaccOut, &ymaccOut, maccOutRem,
			  &zrestIn, &xrestIn, &yrestIn,blksz0,blksz1,restOutrem0,restOutrem1);
		  break;
		case RXCY: 
		  CAT(blk_b22XY_tprod_,TYPESTRRC)
			 ((ZTYPE*)zrest.rp, (ZTYPE*)zrest.ip, 
			  (XTYPE*)xrest.rp, (XTYPE*)xrest.ip, 
			  (YTYPE*)yrest.rp, (YTYPE*)yrest.ip,
			  &zrest, &xrest, &yrest, &xmaccIn,  &ymaccIn,  subs, maccsubs,
			  &xmaccOut, &ymaccOut, maccOutRem,
			  &zrestIn, &xrestIn, &yrestIn,blksz0,blksz1,restOutrem0,restOutrem1);
		  break;
		case CXRY: 
		  CAT(blk_b22XY_tprod_,TYPESTRCR)
			 ((ZTYPE*)zrest.rp, (ZTYPE*)zrest.ip, 
			  (XTYPE*)xrest.rp, (XTYPE*)xrest.ip, 
			  (YTYPE*)yrest.rp, (YTYPE*)yrest.ip,
			  &zrest, &xrest, &yrest, &xmaccIn,  &ymaccIn,  subs, maccsubs,
			  &xmaccOut, &ymaccOut, maccOutRem,
			  &zrestIn, &xrestIn, &yrestIn,blksz0,blksz1,restOutrem0,restOutrem1);
		  break;
		case CXCY: 
		  CAT(blk_b22XY_tprod_,TYPESTRCC) 
			 ((ZTYPE*)zrest.rp, (ZTYPE*)zrest.ip, 
			  (XTYPE*)xrest.rp, (XTYPE*)xrest.ip, 
			  (YTYPE*)yrest.rp, (YTYPE*)yrest.ip,
			  &zrest, &xrest, &yrest, &xmaccIn,  &ymaccIn,  subs, maccsubs,
			  &xmaccOut, &ymaccOut, maccOutRem,
			  &zrestIn, &xrestIn, &yrestIn,blksz0,blksz1,restOutrem0,restOutrem1);
		  break;
		}
		

		delmxInfo(&zrestIn);  delmxInfo(&xrestIn); delmxInfo(&yrestIn);
		delmxInfo(&xmaccIn);  delmxInfo(&ymaccIn);
		delmxInfo(&xmaccOut); delmxInfo(&ymaccOut);
	 }
		
	 /*---------------------------------------------------------------------*/
	 /*couldn't use the fast code, can we use the re-order inner loop code? */
  } else if ( zrest.stride[0]==1 && xrest.stride[0]==1 && 
				  yrest.stride[0]==0 ) { 
		switch (intype){       
		case RXRY: 
		  CAT(ro_tprod_,TYPESTRRR)((ZTYPE*)zrest.rp, (ZTYPE*)zrest.ip, 
											(XTYPE*)xrest.rp, (XTYPE*)xrest.ip, 
											(YTYPE*)yrest.rp, (YTYPE*)yrest.ip,
											&zrest, &xrest, &yrest,
											&xmacc, &ymacc, subs, maccsubs);
		  break;
		case RXCY: 
		  CAT(ro_tprod_,TYPESTRRC)((ZTYPE*)zrest.rp, (ZTYPE*)zrest.ip, 
											(XTYPE*)xrest.rp, (XTYPE*)xrest.ip, 
											(YTYPE*)yrest.rp, (YTYPE*)yrest.ip,
											&zrest, &xrest, &yrest,
											&xmacc, &ymacc, subs, maccsubs);
		  break;
		case CXRY: 
		  CAT(ro_tprod_,TYPESTRCR)((ZTYPE*)zrest.rp, (ZTYPE*)zrest.ip, 
											(XTYPE*)xrest.rp, (XTYPE*)xrest.ip, 
											(YTYPE*)yrest.rp, (YTYPE*)yrest.ip,
											&zrest, &xrest, &yrest,
											&xmacc, &ymacc, subs, maccsubs);
		  break;
		case CXCY: 
		  CAT(ro_tprod_,TYPESTRCC)((ZTYPE*)zrest.rp, (ZTYPE*)zrest.ip, 
											(XTYPE*)xrest.rp, (XTYPE*)xrest.ip, 
											(YTYPE*)yrest.rp, (YTYPE*)yrest.ip,
											&zrest, &xrest, &yrest,
											&xmacc, &ymacc, subs, maccsubs);
		  break;
		}

	 /*-----------------------------------------------------------------------*/
  } else {/* Don't know what else to do so fall back on conventional macc */
	 /* then the normal loop order is the best */
		switch (intype){       
		case RXRY: 
		  CAT(simp_tprod_,TYPESTRRR)((ZTYPE*)zrest.rp, (ZTYPE*)zrest.ip, 
											(XTYPE*)xrest.rp, (XTYPE*)xrest.ip, 
											(YTYPE*)yrest.rp, (YTYPE*)yrest.ip,
											  &zrest, &xrest, &yrest,
											  &xmacc, &ymacc, subs, maccsubs);
		  break;
		case RXCY: 
		  CAT(simp_tprod_,TYPESTRRC)((ZTYPE*)zrest.rp, (ZTYPE*)zrest.ip, 
											  (XTYPE*)xrest.rp, (XTYPE*)xrest.ip, 
											  (YTYPE*)yrest.rp, (YTYPE*)yrest.ip,
											  &zrest, &xrest, &yrest,
											  &xmacc, &ymacc, subs, maccsubs);
		  break;
		case CXRY: 
		  CAT(simp_tprod_,TYPESTRCR)((ZTYPE*)zrest.rp, (ZTYPE*)zrest.ip, 
											  (XTYPE*)xrest.rp, (XTYPE*)xrest.ip, 
											  (YTYPE*)yrest.rp, (YTYPE*)yrest.ip,
											  &zrest, &xrest, &yrest,
											  &xmacc, &ymacc, subs, maccsubs);
		  break;
		case CXCY: 
		  CAT(simp_tprod_,TYPESTRCC)((ZTYPE*)zrest.rp, (ZTYPE*)zrest.ip, 
											  (XTYPE*)xrest.rp, (XTYPE*)xrest.ip, 
											  (YTYPE*)yrest.rp, (YTYPE*)yrest.ip,
											  &zrest, &xrest, &yrest,
											  &xmacc, &ymacc, subs, maccsubs);
		  break;
		}

  }   
  delmxInfo(&zrest);delmxInfo(&xrest);delmxInfo(&yrest);
  delmxInfo(&xmacc);delmxInfo(&ymacc);FREE(maccsubs); 
  return retVal;
}
