/*

 REPOP.C
 Generalized arithmetic operators which implicity replicate their inputs when
 they are not large enough for the input dimensions.

 This file contains the generic c-code which can be used stand-alone
 (i.e. w/o reference to matlab)
 
 * When compilied without a specific operator selected it generates all
 operators, and the associated function to call the selected one.

 * When compilied with a specific operator selected (via -D_OPERATORNM_ ) it
 generates only the appropriate code and utility functions.

 Copyright 2006-     by Jason Farquhar (jdrf@zepler.org)
 Permission is granted for anyone to copy, use, or modify this
 software and accompanying documents for any uncommercial purposes,
 provided this copyright notice is retained, and note is made of
 any changes that have been made. This software and documents are
 distributed without any warranty, express or implied.

 This code is based upon ideas by Douglas M. Schwarz (schwarz@servtech.com)
 and Adi Vehtari -- but heavily modifed to make it more general, easier to
 understand and faster!

TESTCASES

X=randn(100,1000); Y=randn(size(X));
X=complex(randn(100,1000),randn(100,1000)); Y=complex(randn(size(X)),randn(size(X)));
mimage(Z,T,Z-T);
% PLUS
Z=repop(X,Y(:,1),'+');   T=X+repmat(Y(:,1),[1,size(X,2)]);  
Z=repop(X,Y(:,1:2),'+'); T=X+repmat(Y(:,1:2),[1,size(X,2)/2]);
Z=repop(X,Y(1,:),'+');   T=X+repmat(Y(1,:),[size(X,1),1]);
Z=repop(X,Y(:,:),'+');   T=X+repmat(Y(:,:),[1,1]); 
% TIMES
Z=repop(X,Y(:,1),'*');   T=X.*repmat(Y(:,1),[1,size(X,2)]);
Z=repop(X,Y(:,1:2),'*'); T=X.*repmat(Y(:,1:2),[1,size(X,2)/2]);
Z=repop(X,Y(1,:),'*');   T=X.*repmat(Y(1,:),[size(X,1),1]);
Z=repop(X,Y(:,:),'*');   T=X.*repmat(Y(:,:),[1,1]);  

*/

#include <string.h>

#include "mxInfo.h"
#include "repop.h"


/*---------------------------------------------------------------------------*
 * use the macro-preprocessor and the repops.template file to automatically *
 * generate the functions to call for the different combinations of real and
 * complex inputs */
/* decide if we're building 1 or all operators */
#if defined(_PLUS_) ||defined(_MINUS_) || defined(_TIMES_) \
    || defined(_RDIVIDE_) || defined(_LDIVIDE_) || defined(_POWER_)\
    || defined(_EQ_) || defined(_NE_) || defined(_LT_) || defined(_GT_) \
    || defined(_LE_) || defined(_GE_) \
    || defined(_MAX_) || defined(_MIN_) || defined(_MOD_)
#undef ALLOPS
#else 
#define ALLOPS
#endif
/* --------------------------------------------------------------------*/

/* we default to double matrix types */
/*---------------------------------------------------------------------------*/
/* Set the default input type if not defined */
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
#define OZDTYPE DOUBLE_DTYPE
#define OZTYPE double
#elif (XDTYPE == DOUBLE_DTYPE) && (YDTYPE == SINGLE_DTYPE)
#define TYPESTR ds
#define XTYPE double
#define YTYPE float
#define OZDTYPE SINGLE_DTYPE
#define OZTYPE float
#elif (XDTYPE == SINGLE_DTYPE) && (YDTYPE == DOUBLE_DTYPE)
#define TYPESTR sd
#define XTYPE float
#define YTYPE double
#define OZDTYPE SINGLE_DTYPE
#define OZTYPE float
#elif (XDTYPE == SINGLE_DTYPE) && (YDTYPE == SINGLE_DTYPE)
#define TYPESTR ss
#define XTYPE float
#define YTYPE float
#define OZDTYPE SINGLE_DTYPE
#define OZTYPE float
#elif (XDTYPE== INT32_DTYPE ) && (YDTYPE== INT32_DTYPE)
#define TYPESTR ii
#define XTYPE int
#define YTYPE int
#define OZDTYPE INT32_DTYPE
#define OZTYPE int
#else
#error Unsupported X,Y type combination
#endif

#define ZTYPE  OZTYPE
#define ZDTYPE OZDTYPE

/* utility defines for constructing the function to call */
#define TYPESTRRR CAT(TYPESTR,rxry)
#define TYPESTRRC CAT(TYPESTR,rxcy)
#define TYPESTRCR CAT(TYPESTR,cxry)
#define TYPESTRCC CAT(TYPESTR,cxcy)

/* Generate the code for each of the different operators */
#if defined(_PLUS_) || defined(ALLOPS)
#define OPID PLUS
#define OPNM plus
#define XISCOMPLEX
#define YISCOMPLEX
#include "repop.def"   /* complex x, complex y */
#undef YISCOMPLEX
#define XISCOMPLEX
#include "repop.def"   /* complex x, real y */
#undef XISCOMPLEX
#define YISCOMPLEX
#include "repop.def"   /* real x, complex y */
#undef XISCOMPLEX
#undef YISCOMPLEX
#include "repop.def"   /* real x, real y */
#endif

#if defined(_MINUS_) || defined(ALLOPS)
#undef OPNM
#undef OPID
#define OPID MINUS
#define OPNM minus
#define XISCOMPLEX
#define YISCOMPLEX
#include "repop.def"   /* complex x, complex y */
#undef YISCOMPLEX
#define XISCOMPLEX
#include "repop.def"   /* complex x, real y */
#undef XISCOMPLEX
#define YISCOMPLEX
#include "repop.def"   /* real x, complex y */
#undef XISCOMPLEX
#undef YISCOMPLEX
#include "repop.def"   /* real x, real y */
#endif

#if defined(_TIMES_) || defined(ALLOPS)
#undef OPNM
#undef OPID
#define OPID TIMES
#define OPNM times
#define XISCOMPLEX
#define YISCOMPLEX
#include "repop.def"   /* complex x, complex y */
#undef YISCOMPLEX
#define XISCOMPLEX
#include "repop.def"   /* complex x, real y */
#undef XISCOMPLEX
#define YISCOMPLEX
#include "repop.def"   /* real x, complex y */
#undef XISCOMPLEX
#undef YISCOMPLEX
#include "repop.def"   /* real x, real y */
#endif

#if defined(_RDIVIDE_) || defined(ALLOPS)
#undef OPNM
#undef OPID
#define OPID RDIVIDE
#define OPNM rdivide
#define XISCOMPLEX
#define YISCOMPLEX
#include "repop.def"   /* complex x, complex y */
#undef YISCOMPLEX
#define XISCOMPLEX
#include "repop.def"   /* complex x, real y */
#undef XISCOMPLEX
#define YISCOMPLEX
#include "repop.def"   /* real x, complex y */
#undef XISCOMPLEX
#undef YISCOMPLEX
#include "repop.def"   /* real x, real y */
#endif

#if defined(_LDIVIDE_) || defined(ALLOPS)
#undef OPNM
#undef OPID
#define OPID LDIVIDE
#define OPNM ldivide
#define XISCOMPLEX
#define YISCOMPLEX
#include "repop.def"   /* complex x, complex y */
#undef YISCOMPLEX
#define XISCOMPLEX
#include "repop.def"   /* complex x, real y */
#undef XISCOMPLEX
#define YISCOMPLEX
#include "repop.def"   /* real x, complex y */
#undef XISCOMPLEX
#undef YISCOMPLEX
#include "repop.def"   /* real x, real y */
#endif

#if defined(_POWER_) || defined(ALLOPS)
#undef OPNM
#undef OPID
#define OPID POWER
#define OPNM power
#define XISCOMPLEX
#define YISCOMPLEX
#include "repop.def"   /* complex x, complex y */
#undef YISCOMPLEX
#define XISCOMPLEX
#include "repop.def"   /* complex x, real y */
#undef XISCOMPLEX
#define YISCOMPLEX
#include "repop.def"   /* real x, complex y */
#undef XISCOMPLEX
#undef YISCOMPLEX
#include "repop.def"   /* real x, real y */
#endif

#if defined(_EQ_) || defined(ALLOPS)
#undef OPNM
#undef OPID
#define OPID EQ
#define OPNM eq
/* re-define the ZTYPE to logical */
#undef ZTYPE
#define ZTYPE  BOOLTYPE
#undef ZDTYPE
#define ZDTYPE LOGICAL_DTYPE

#define XISCOMPLEX
#define YISCOMPLEX
#include "repop.def"   /* complex x, complex y */
#undef YISCOMPLEX
#define XISCOMPLEX
#include "repop.def"   /* complex x, real y */
#undef XISCOMPLEX
#define YISCOMPLEX
#include "repop.def"   /* real x, complex y */
#undef XISCOMPLEX
#undef YISCOMPLEX
#include "repop.def"   /* real x, real y */
/* reset the z-type */
#undef ZTYPE
#define ZTYPE  OZTYPE
#undef ZDTYPE
#define ZDTYPE OZDTYPE
#endif

#if defined(_NE_) || defined(ALLOPS)
#undef OPNM
#undef OPID
#define OPID NE
#define OPNM ne
/* re-define the ZTYPE to logical */
#undef ZTYPE
#define ZTYPE  BOOLTYPE
#undef ZDTYPE
#define ZDTYPE LOGICAL_DTYPE

#define XISCOMPLEX
#define YISCOMPLEX
#include "repop.def"   /* complex x, complex y */
#undef YISCOMPLEX
#define XISCOMPLEX
#include "repop.def"   /* complex x, real y */
#undef XISCOMPLEX
#define YISCOMPLEX
#include "repop.def"   /* real x, complex y */
#undef XISCOMPLEX
#undef YISCOMPLEX
#include "repop.def"   /* real x, real y */
/* reset the z-type */
#undef ZTYPE
#define ZTYPE  OZTYPE
#undef ZDTYPE
#define ZDTYPE OZDTYPE
#endif


#if defined(_LT_) || defined(ALLOPS)
#undef OPNM
#undef OPID
#define OPID LT
#define OPNM lt
/* re-define the ZTYPE to logical */
#undef ZTYPE
#define ZTYPE  BOOLTYPE
#undef ZDTYPE
#define ZDTYPE LOGICAL_DTYPE

#define XISCOMPLEX
#define YISCOMPLEX
#include "repop.def"   /* complex x, complex y */
#undef YISCOMPLEX
#define XISCOMPLEX
#include "repop.def"   /* complex x, real y */
#undef XISCOMPLEX
#define YISCOMPLEX
#include "repop.def"   /* real x, complex y */
#undef XISCOMPLEX
#undef YISCOMPLEX
#include "repop.def"   /* real x, real y */
/* reset the z-type */
#undef ZTYPE
#define ZTYPE  OZTYPE
#undef ZDTYPE
#define ZDTYPE OZDTYPE
#endif

#if defined(_GT_) || defined(ALLOPS)
#undef OPNM
#undef OPID
#define OPID GT
#define OPNM gt
/* re-define the ZTYPE to logical */
#undef ZTYPE
#define ZTYPE  BOOLTYPE
#undef ZDTYPE
#define ZDTYPE LOGICAL_DTYPE

#define XISCOMPLEX
#define YISCOMPLEX
#include "repop.def"   /* complex x, complex y */
#undef YISCOMPLEX
#define XISCOMPLEX
#include "repop.def"   /* complex x, real y */
#undef XISCOMPLEX
#define YISCOMPLEX
#include "repop.def"   /* real x, complex y */
#undef XISCOMPLEX
#undef YISCOMPLEX
#include "repop.def"   /* real x, real y */
/* reset the z-type */
#undef ZTYPE
#define ZTYPE  OZTYPE
#undef ZDTYPE
#define ZDTYPE OZDTYPE
#endif

#if defined(_LE_) || defined(ALLOPS)
#undef OPNM
#undef OPID
#define OPID LE
#define OPNM le
/* re-define the ZTYPE to logical */
#undef ZTYPE
#define ZTYPE  BOOLTYPE
#undef ZDTYPE
#define ZDTYPE LOGICAL_DTYPE

#define XISCOMPLEX
#define YISCOMPLEX
#include "repop.def"   /* complex x, complex y */
#undef YISCOMPLEX
#define XISCOMPLEX
#include "repop.def"   /* complex x, real y */
#undef XISCOMPLEX
#define YISCOMPLEX
#include "repop.def"   /* real x, complex y */
#undef XISCOMPLEX
#undef YISCOMPLEX
#include "repop.def"   /* real x, real y */
/* reset the z-type */
#undef ZTYPE
#define ZTYPE  OZTYPE
#undef ZDTYPE
#define ZDTYPE OZDTYPE
#endif

#if defined(_GE_) || defined(ALLOPS)
#undef OPNM
#undef OPID
#define OPID GE
#define OPNM ge
/* re-define the ZTYPE to logical */
#undef ZTYPE
#define ZTYPE  BOOLTYPE
#undef ZDTYPE
#define ZDTYPE LOGICAL_DTYPE

#define XISCOMPLEX
#define YISCOMPLEX
#include "repop.def"   /* complex x, complex y */
#undef YISCOMPLEX
#define XISCOMPLEX
#include "repop.def"   /* complex x, real y */
#undef XISCOMPLEX
#define YISCOMPLEX
#include "repop.def"   /* real x, complex y */
#undef XISCOMPLEX
#undef YISCOMPLEX
#include "repop.def"   /* real x, real y */
/* reset the z-type */
#undef ZTYPE
#define ZTYPE  OZTYPE
#undef ZDTYPE
#define ZDTYPE OZDTYPE
#endif

#if defined(_MINOP_) || defined(ALLOPS)
#undef OPNM
#undef OPID
#define OPID MINOP
#define OPNM min
#define XISCOMPLEX
#define YISCOMPLEX
#include "repop.def"   /* complex x, complex y */
#undef YISCOMPLEX
#define XISCOMPLEX
#include "repop.def"   /* complex x, real y */
#undef XISCOMPLEX
#define YISCOMPLEX
#include "repop.def"   /* real x, complex y */
#undef XISCOMPLEX
#undef YISCOMPLEX
#include "repop.def"   /* real x, real y */
#endif

#if defined(_MAXOP_) || defined(ALLOPS)
#undef OPNM
#undef OPID
#define OPID MAXOP
#define OPNM max
#define XISCOMPLEX
#define YISCOMPLEX
#include "repop.def"   /* complex x, complex y */
#undef YISCOMPLEX
#define XISCOMPLEX
#include "repop.def"   /* complex x, real y */
#undef XISCOMPLEX
#define YISCOMPLEX
#include "repop.def"   /* real x, complex y */
#undef XISCOMPLEX
#undef YISCOMPLEX
#include "repop.def"   /* real x, real y */
#endif

#if defined(_MODOP_) || defined(ALLOPS)
#undef OPNM
#undef OPID
#define OPID MOD
#define OPNM mod
#define XISCOMPLEX
#define YISCOMPLEX
#include "repop.def"   /* complex x, complex y */
#undef YISCOMPLEX
#define XISCOMPLEX
#include "repop.def"   /* complex x, real y */
#undef XISCOMPLEX
#define YISCOMPLEX
#include "repop.def"   /* real x, complex y */
#undef XISCOMPLEX
#undef YISCOMPLEX
#include "repop.def"   /* real x, real y */
#endif

#undef OPNM
#undef OPID


/****************************************************************************
 * the main work function 
 ***************************************************************************/
/* only compile if someone hasn't asked for a single operator */
#if defined(ALLOPS) 

RepopErrorCode CAT(TYPESTR,repop)(MxInfo zinfoorig, 
											 const MxInfo xinfoorig, 
											 const MxInfo yinfoorig, 
											 int opid){
  RepopErrorCode retVal=OK ;
  MxInfo xinfo, yinfo, zinfo; /* make copy of input x/y info to modify here */
  int i, intype;
  int *repxy   = MALLOC(zinfoorig.nd*sizeof(int));

  /* check the types match */
  if ( xinfoorig.dtype != XDTYPE ) {
	 retVal = XTYPEMISMATCH;
  } else if ( yinfoorig.dtype != YDTYPE ) {
	 retVal = YTYPEMISMATCH; 
  }
  if ( retVal != OK ) { /* something didn't match return error code */
	 FREE(repxy); return retVal ;
  }

  /* Types are OK so do the operation */
  zinfo=copymxInfo(zinfoorig);
  xinfo=copymxInfo(xinfoorig);
  yinfo=copymxInfo(yinfoorig);
  
  /* identify the rep dims */
  for( i=0; i<zinfo.nd; i++ ){
	 /* repxy indicates which of x/y to rep: 0=none, <0=rep X, >0 rep Y */
	 repxy[i]= xinfo.sz[i]==yinfo.sz[i]?0:(xinfo.sz[i]<yinfo.sz[i]?1:-1);
	 /* as an efficiency hack -- if rep-size == 1 then set its stride=0 
		 and mark as non-repd to save the test and step back overhead */
	 if ( repxy[i]>0 && xinfo.sz[i]==1 ) { repxy[i]=0; xinfo.stride[i]=0; }
	 if ( repxy[i]<0 && yinfo.sz[i]==1 ) { repxy[i]=0; yinfo.stride[i]=0; }	 
  }
  
  repopqueryOptimise(&zinfo,&xinfo,&yinfo,repxy);

  /* enum for the input type, rxry=0,cxry=1,rxcy=2,cxcy=3 */
  intype = (xinfo.ip==0?0:1) + (yinfo.ip==0?0:2); 
 
  /*-------------------------------------------------------------------*/
  /* call the appropriate function to compute the result for this
	  operator and type of inputs */
  switch (intype){ 
  case RXRY: /* real x, real y */
	 switch(opid){
	 case PLUS:    retVal=CAT(repplus_,TYPESTRRR)   (zinfo,xinfo,yinfo,repxy);  break;
	 case MINUS:   retVal=CAT(repminus_,TYPESTRRR)  (zinfo,xinfo,yinfo,repxy);  break;
	 case TIMES:   retVal=CAT(reptimes_,TYPESTRRR)  (zinfo,xinfo,yinfo,repxy);  break;
	 case RDIVIDE: retVal=CAT(reprdivide_,TYPESTRRR)(zinfo,xinfo,yinfo,repxy);  break;
	 case LDIVIDE: retVal=CAT(repldivide_,TYPESTRRR)(zinfo,xinfo,yinfo,repxy);  break;
	 case POWER:   retVal=CAT(reppower_,TYPESTRRR)  (zinfo,xinfo,yinfo,repxy);  break;
	 case MOD:     retVal=CAT(repmod_,TYPESTRRR)    (zinfo,xinfo,yinfo,repxy);  break;
	 case EQ:      retVal=CAT(repeq_,TYPESTRRR)     (zinfo,xinfo,yinfo,repxy);  break;
	 case NE:      retVal=CAT(repne_,TYPESTRRR)     (zinfo,xinfo,yinfo,repxy);  break;
	 case LT:      retVal=CAT(replt_,TYPESTRRR)     (zinfo,xinfo,yinfo,repxy);  break;
	 case GT:      retVal=CAT(repgt_,TYPESTRRR)     (zinfo,xinfo,yinfo,repxy);  break;
	 case LE:      retVal=CAT(reple_,TYPESTRRR)     (zinfo,xinfo,yinfo,repxy);  break;
	 case GE:      retVal=CAT(repge_,TYPESTRRR)     (zinfo,xinfo,yinfo,repxy);  break;
	 case MINOP:   retVal=CAT(repmin_,TYPESTRRR)    (zinfo,xinfo,yinfo,repxy);  break;
	 case MAXOP:   retVal=CAT(repmax_,TYPESTRRR)    (zinfo,xinfo,yinfo,repxy);  break;
	 default: retVal=UNDEFOPERATOR;
	 }
	 break;
  case CXRY: /* complex x, real y */
	 switch(opid){ 
		/* N.B. for mixed cases can treat as real/real as complex will be
			unaffected.  But only if we memcpy+replicate accross the complex
			part after! */ 
	 case PLUS:		retVal=CAT(repplus_,TYPESTRCR)   (zinfo,xinfo,yinfo,repxy);  break;
	 case MINUS:   retVal=CAT(repminus_,TYPESTRCR)  (zinfo,xinfo,yinfo,repxy);  break;
	 case TIMES: 
		/* 2 separate multi's is marginaly faster if cy is large (so r/im
			parts interfere in the cache?), but slower when its small... at best
			only 2% better to ignore */
	   retVal=CAT(reptimes_,TYPESTRCR)(zinfo,xinfo,yinfo,repxy);
		break;
	 case RDIVIDE: retVal=CAT(reprdivide_,TYPESTRCR)(zinfo,xinfo,yinfo,repxy);  break;
	 case LDIVIDE: retVal=CAT(repldivide_,TYPESTRCR)(zinfo,xinfo,yinfo,repxy);  break;
	 case POWER:   retVal=CAT(reppower_,TYPESTRCR)  (zinfo,xinfo,yinfo,repxy);  break;
	 case MOD:     retVal=CAT(repmod_,TYPESTRCR)  (zinfo,xinfo,yinfo,repxy);  break;
	 case EQ:      retVal=CAT(repeq_,TYPESTRCR)     (zinfo,xinfo,yinfo,repxy);  break;
	 case NE:      retVal=CAT(repne_,TYPESTRCR)     (zinfo,xinfo,yinfo,repxy);  break;
	 case LT:      retVal=CAT(replt_,TYPESTRCR)     (zinfo,xinfo,yinfo,repxy);  break;
	 case GT:      retVal=CAT(repgt_,TYPESTRCR)     (zinfo,xinfo,yinfo,repxy);  break;
	 case LE:      retVal=CAT(reple_,TYPESTRCR)     (zinfo,xinfo,yinfo,repxy);  break;
	 case GE:      retVal=CAT(repge_,TYPESTRCR)     (zinfo,xinfo,yinfo,repxy);  break;
	 case MINOP:   retVal=CAT(repmin_,TYPESTRCR)    (zinfo,xinfo,yinfo,repxy);  break;
	 case MAXOP:   retVal=CAT(repmax_,TYPESTRCR)    (zinfo,xinfo,yinfo,repxy);  break;
	 default: retVal=UNDEFOPERATOR;
	 }
	 break;
  case RXCY: /* real x, complex y */
	 switch(opid){
		/* N.B. for mixed cases can treat as real/real as complex will be
			unaffected. But only if we memcpy+replicate accross the complex part
			after! */ 
	 case PLUS:    retVal=CAT(repplus_,TYPESTRRC)   (zinfo,xinfo,yinfo,repxy);  break;
		/* N.B. must treat y as complex to change sign of complex part */
	 case MINUS:   retVal=CAT(repminus_,TYPESTRRC)  (zinfo,xinfo,yinfo,repxy);  break;
	 case TIMES:  
		/* 2 separate multi's is marginaly faster if cy is large (so r/im
			parts interfere in the cache?), but slower when its small... at best
			only 2% better to ignore */
	   retVal=CAT(reptimes_,TYPESTRRC)(zinfo,xinfo,yinfo,repxy);
		break;
	 case RDIVIDE: retVal=CAT(reprdivide_,TYPESTRRC)(zinfo,xinfo,yinfo,repxy);  break;
	 case LDIVIDE: retVal=CAT(repldivide_,TYPESTRRC)(zinfo,xinfo,yinfo,repxy);  break;
	 case POWER:   retVal=CAT(reppower_,TYPESTRRC)  (zinfo,xinfo,yinfo,repxy);  break;
	 case MOD:     retVal=CAT(repmod_,TYPESTRRC)  (zinfo,xinfo,yinfo,repxy);  break;
	 case EQ:      retVal=CAT(repeq_,TYPESTRRC)     (zinfo,xinfo,yinfo,repxy);  break;
	 case NE:      retVal=CAT(repne_,TYPESTRRC)     (zinfo,xinfo,yinfo,repxy);  break;
	 case LT:      retVal=CAT(replt_,TYPESTRRC)     (zinfo,xinfo,yinfo,repxy);  break;
	 case GT:      retVal=CAT(repgt_,TYPESTRRC)     (zinfo,xinfo,yinfo,repxy);  break;
	 case LE:      retVal=CAT(reple_,TYPESTRRC)     (zinfo,xinfo,yinfo,repxy);  break;
	 case GE:      retVal=CAT(repge_,TYPESTRRC)     (zinfo,xinfo,yinfo,repxy);  break;
	 case MINOP:   retVal=CAT(repmin_,TYPESTRRC)    (zinfo,xinfo,yinfo,repxy);  break;
	 case MAXOP:   retVal=CAT(repmax_,TYPESTRRC)    (zinfo,xinfo,yinfo,repxy);  break;
	 default: retVal=UNDEFOPERATOR;
	 }
	 break;
  case CXCY: /* complex x, complex y */
	 switch(opid){
	 case PLUS:  /* marginaly faster to do as 2 separate additions */
		retVal=CAT(repplus_,TYPESTRRR)   (zinfo,xinfo,yinfo,repxy);
		zinfo.rp=zinfo.ip; xinfo.rp=xinfo.ip; yinfo.rp=yinfo.ip;
		retVal=CAT(repplus_,TYPESTRRR)   (zinfo,xinfo,yinfo,repxy);  
		break;
	 case MINUS: /* marginaly faster to do as 2 separate subtractions */
		retVal=CAT(repminus_,TYPESTRRR)   (zinfo,xinfo,yinfo,repxy);
		zinfo.rp=zinfo.ip; xinfo.rp=xinfo.ip; yinfo.rp=yinfo.ip;
		retVal=CAT(repminus_,TYPESTRRR)   (zinfo,xinfo,yinfo,repxy);  
		break;		
	 case TIMES:   retVal=CAT(reptimes_,TYPESTRCC)  (zinfo,xinfo,yinfo,repxy);  break;
	 case RDIVIDE: retVal=CAT(reprdivide_,TYPESTRCC)(zinfo,xinfo,yinfo,repxy);  break;
	 case LDIVIDE: retVal=CAT(repldivide_,TYPESTRCC)(zinfo,xinfo,yinfo,repxy);  break;
	 case POWER:   retVal=CAT(reppower_,TYPESTRCC)  (zinfo,xinfo,yinfo,repxy);  break;
	 case MOD:     retVal=CAT(repmod_,TYPESTRCC)  (zinfo,xinfo,yinfo,repxy);  break;
	 case EQ:      retVal=CAT(repeq_,TYPESTRCC)     (zinfo,xinfo,yinfo,repxy);  break;
	 case NE:      retVal=CAT(repne_,TYPESTRCC)     (zinfo,xinfo,yinfo,repxy);  break;
	 case LT:      retVal=CAT(replt_,TYPESTRCC)     (zinfo,xinfo,yinfo,repxy);  break;
	 case GT:      retVal=CAT(repgt_,TYPESTRCC)     (zinfo,xinfo,yinfo,repxy);  break;
	 case LE:      retVal=CAT(reple_,TYPESTRCC)     (zinfo,xinfo,yinfo,repxy);  break;
	 case GE:      retVal=CAT(repge_,TYPESTRCC)     (zinfo,xinfo,yinfo,repxy);  break;
	 case MINOP:   retVal=CAT(repmin_,TYPESTRCC)    (zinfo,xinfo,yinfo,repxy);  break;
	 case MAXOP:   retVal=CAT(repmax_,TYPESTRCC)    (zinfo,xinfo,yinfo,repxy);  break;
	 default: retVal=UNDEFOPERATOR;
	 }
	 break;
  default:
	 retVal = INTYPEMISMATCH ;
  }

  FREE(repxy);
  delmxInfo(&xinfo); delmxInfo(&yinfo); delmxInfo(&zinfo);
  return retVal; 
}
#endif
