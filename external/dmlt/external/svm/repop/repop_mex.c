/*

 REPOP.C
 Generalized arithmetic operators which implicity replicate their inputs when
 they are not large enough for the input dimensions.

 This file contains the matlab wrapper code -- relying on repop.c to do the
 main work.

 Copyright 2006-     by Jason Farquhar (jdrf@zepler.org)

 Permission is granted for anyone to copy, use, or modify this
 software and accompanying documents for any uncommercial purposes,
 provided this copyright notice is retained, and note is made of
 any changes that have been made. This software and documents are
 distributed without any warranty, express or implied.

 This code is based upon ideas by Douglas M. Schwarz (schwarz@servtech.com)
 and Adi Vehtari -- but heavily modifed to make it more general, easier to
 understand and faster!

 $Id$

TESTCASES

X=randn(100,1000); Y=randn(size(X));
X=complex(randn(100,1000),randn(100,1000)); Y=complex(randn(size(X)),randn(size(X)));
mimage(Z,T,Z-T);
% PLUS
Z=repop(X,'+',Y(:,1));   T=X+repmat(Y(:,1),[1,size(X,2)]);  
Z=repop(X,'+',Y(:,1:2)); T=X+repmat(Y(:,1:2),[1,size(X,2)/2]);
Z=repop(X,'+',Y(1,:));   T=X+repmat(Y(1,:),[size(X,1),1]);
Z=repop(X,'+',Y(:,:));   T=X+repmat(Y(:,:),[1,1]); 
% TIMES
Z=repop(X,'*',Y(:,1));   T=X.*repmat(Y(:,1),[1,size(X,2)]);
Z=repop(X,'*',Y(:,1:2)); T=X.*repmat(Y(:,1:2),[1,size(X,2)/2]);
Z=repop(X,'*',Y(1,:));   T=X.*repmat(Y(1,:),[size(X,1),1]);
Z=repop(X,'*',Y(:,:));   T=X.*repmat(Y(:,:),[1,1]);  

*/

#include "mex.h"
#include "matrix.h"
#include "mxInfo.h"
#include "mxInfo_mex.h"
#include "string.h"
#include "repop.h"


/* Function prototype for undocumented in-place utiltity code */
extern mxArray *mxCreateSharedDataCopy(const mxArray *pr);
extern bool     mxUnshareArray(const mxArray *pr, const bool noDeepCopy);
extern mxArray *mxUnreference(const mxArray *pr);


/*---------------------------------------------------------------------------*
 * use the macro-preprocessor and the repops.template file to automatically *
 * generate the functions to call for the different combinations of real and
 * complex inputs */

#if defined(_PLUS_) ||defined(_MINUS_) || defined(_TIMES_) || defined(_POWER_)\
    || defined(_EQ_) || defined(_NE_) || defined(_LT_) || defined(_GT_) \
    || defined(_LE_) || defined(_GE_) || defined(_MINOP_) || defined(_MAXOP_)\
    || defined(_MOD_)
#error Only define the operator to use for repop_ind.c
#endif

/*-------------------------------------------------------------------------*/
/* the wrapper function */
void mexFunction(const int nlhs, mxArray *plhs[],
                 const int nrhs, const mxArray *prhs[]) {
  RepopErrorCode retVal=OK ;
  char *opnmStr;
  int xnd, ynd, znd, opid, xargi=-1, yargi=-1, opargi=-1, i=0;
  MxInfo xinfo, yinfo, zinfo;
  int repNonUnitDim=0, inPlace=false;
                     
  if (nrhs < 3) ERROR("repop: Incorrect number of inputs.");        
  if (nlhs > 1) ERROR("repop: Too many output arguments.");        

  if( nrhs>3 ){ /* parse the options */
	 if ( mxIsChar(prhs[3]) ){
		opnmStr=mxArrayToString(prhs[3]);
		for (i=0 ; opnmStr[i] != 0; i++ ) {
		  switch ( opnmStr[i] ) {
		  case 'i': case 'I':   inPlace=1;  break;
		  case 'm': case 'M':   repNonUnitDim=1;  break;
		  case 'n': case 'N':   repNonUnitDim=2;  break;
		  default: WARNING("repop: Unrecognised tprod option");
		  }
		}
		mxFree(opnmStr);
	 }
  }

  /* first -- find the operator position and decipher which operator we're
	  supposed to use */
  if( mxIsChar(prhs[0])) { xargi=1; yargi=2; opargi=0; }
  if( mxIsChar(prhs[1])) { xargi=0; yargi=2; opargi=1; }
  if( mxIsChar(prhs[2])) { xargi=0; yargi=1; opargi=2; }
  if( opargi<0 || !mxIsChar(prhs[opargi])) mexErrMsgTxt("repop: one arg must be operator type"); 
  opnmStr=mxArrayToString(prhs[opargi]);
  opid = getOpid(opnmStr);  
  if ( opid < 0 ) 
	 mexErrMsgTxt("repop:Unrecognised operator. Must be one of: + - * / \\ ^ % = < > ~= <= >= min max");

  xnd = mxGetNumberOfDimensions(prhs[xargi]);
  ynd = mxGetNumberOfDimensions(prhs[yargi]);
  znd = max(xnd,ynd);
  xinfo=mkmxInfoMxArray(prhs[xargi],znd);
  yinfo=mkmxInfoMxArray(prhs[yargi],znd);

  /* empty set as input for y means make it copy of x */
  if ( yinfo.numel==0 && yinfo.nd==2 && yinfo.sz[0]==0 && yinfo.sz[1]==0 ) { 
	 yinfo=copymxInfo(xinfo); 
  }

  /* check the types are what we can use */
  if ( !(xinfo.dtype == DOUBLE_DTYPE || xinfo.dtype == SINGLE_DTYPE ) ){
	 ERROR("tprod: X type unsuppored: only full double/single");
  }
  if ( mxIsSparse(prhs[xargi]) ){
	 ERROR("tprod: X is sparse, only full double/single supported");
  }
  if ( !(yinfo.dtype == DOUBLE_DTYPE || yinfo.dtype == SINGLE_DTYPE ) ){
	 ERROR("tprod: Y type unsuppored: only double, single");
  }
  if ( mxIsSparse(prhs[yargi]) ){
	 ERROR("tprod: X is sparse, only full double/single supported");
  }

  #ifdef LOGGING  
  FILE *fd = fopen("/tmp/repop.log","a");
  if( fd==0 ){WARNING("Couldn't open log file /tmp/repop.log\n"); fd=stderr;}
  fprintf(fd,"repop( "); printMxInfoSummary(fd,xinfo);
  fprintf(fd,", \'%s\', ",opnmStr); 
  printMxInfoSummary(fd,yinfo);fprintf(fd,")\n");
  if ( fd!=stderr ) fclose(fd);
  #endif

  /*MxInfo zinfo; */
  /* deal with empty inputs */  
  if ( xinfo.numel == 0 || yinfo.numel==0 ) {
	 zinfo=mkemptymxInfo(2);
	 plhs[0] = mxCreateNumericArray(zinfo.nd,zinfo.sz,mxDOUBLE_CLASS,mxREAL); 

  } else {

	 /* compute the z size */
	 zinfo=initzinfo(xinfo,yinfo,repNonUnitDim);

	 if ( inPlace ) { /* test if requested in-place operation is possible */
		if ( zinfo.nd != xinfo.nd ) {
		  inPlace=0;
		} else {
		  for ( i=0; i<xinfo.nd; i++ ) /* check all dimensions sizes match */
			 if ( zinfo.sz[i] != xinfo.sz[i] ) { inPlace=false; break; }
		}
		/* test if X's complex status matches the result */
		if ( (xinfo.ip==0 && yinfo.ip!=0) || /* y-complex & x-real */
			  (xinfo.ip==0 && opid==POWER ) ) /* x-real & power operator */
		  inPlace=0; 
		if ( xinfo.dtype!=zinfo.dtype ) inPlace = 0 ; /* not x,y differ */
		if ( !inPlace ) 
		  WARNING("repop: InPlace failure: size(Z)~=size(X) | isreal(Z)~=isreal(X) | isa(Z)~=isa(X)");
	 }

	 /* allocate the space for z */
	 if( opid >= EQ && opid <= GE ) {
		/*logical operator so make logical output array*/
		plhs[0]  = mxCreateLogicalArray(zinfo.nd, zinfo.sz);
		zinfo.rp = (double*)mxGetLogicals(plhs[0]);

	 } else { /* not a logical operator so Z is same type as X */
		if ( opid == POWER || xinfo.ip != 0 || yinfo.ip != 0) {
		  zinfo.ip = (double*)1 ; /* set ip to indicate need complex array */
		}
		if ( inPlace ) { /* use X inplace */

		  /* Simplier alternative? */
		  mxUnshareArray(prhs[xargi],0); /* ensure X isnt deep copy */
		  plhs[0] = mxCreateSharedDataCopy(prhs[xargi]);

/* 		  /\* create a empty temp array & set X to point to it  *\/ */
/* 		  plhs[0] = mxCreateNumericMatrix(0,0,zinfo.dtype, */
/* 													 zinfo.ip?mxCOMPLEX:mxREAL); */
/* 		  mxSetPr((mxArray*)prhs[xargi],mxGetPr(plhs[0]));  */
/* 		  mxSetPi((mxArray*)prhs[xargi],mxGetPi(plhs[0])); */
/* 		  mxSetDimensions((mxArray*)prhs[xargi],mxGetDimensions(plhs[0]),2); */
/* 		  /\* set Z to use Xs old memory instead *\/ */
/* 		  mxSetPr(plhs[0],xinfo.rp);		   */
/* 		  mxSetPi(plhs[0],xinfo.ip); */
/* 		  mxSetDimensions(plhs[0],zinfo.sz,zinfo.nd); */
		} else { /* allocate the output array info+mem */
		  plhs[0] = mxCreateNumericArray(zinfo.nd, zinfo.sz, zinfo.dtype, 
													zinfo.ip?mxCOMPLEX:mxREAL);
		}
		zinfo.rp = mxGetPr(plhs[0]);  
		zinfo.ip = mxGetPi(plhs[0]);
	 }
	 zinfo.dtype = mxGetClassID(plhs[0]); /* set the type indicator */

	 /* do the (appropriately typed) operation */
	 if(        xinfo.dtype==DOUBLE_DTYPE && yinfo.dtype==DOUBLE_DTYPE ) {/*dd*/
		retVal = ddrepop( zinfo,  xinfo,  yinfo, opid ); 
	 
	 } else if( xinfo.dtype==DOUBLE_DTYPE && yinfo.dtype==SINGLE_DTYPE ) {/*ds*/
		retVal = dsrepop( zinfo,  xinfo,  yinfo, opid ); 
	 
	 } else if( xinfo.dtype==SINGLE_DTYPE && yinfo.dtype==DOUBLE_DTYPE ) {/*sd*/
		retVal = sdrepop( zinfo,  xinfo,  yinfo, opid ); 
	 
	 } else if( xinfo.dtype==SINGLE_DTYPE && yinfo.dtype==SINGLE_DTYPE ){/*ss*/
		retVal = ssrepop( zinfo,  xinfo,  yinfo, opid ); 
	 
	 } else {
		ERROR("tprod: Inputs of unsupported type: only double/single");		
	 }
	 if ( opid == POWER ) removeZeroImag(&zinfo); /* pure real if poss */

	 /* check for errors */
	 switch ( retVal ) {
	 case ZTYPEMISMATCH : 
		ERROR("repop: Z is of unsupported type"); break;
	 case XTYPEMISMATCH :
		ERROR("repop: X is of unsupported type"); break;
	 case YTYPEMISMATCH :
		ERROR("repop: Y is of unsupported type"); break;
	 case INTYPEMISMATCH :
		ERROR("repop: input real/complex mix unsupported"); break;
	 case UNDEFOPERATOR :
		ERROR("repop: Unsupported operator requested"); break;
	 case OTHERERROR :
		ERROR("repop: Something else went wrong, sorry!"); break;
	 default: ;
	 }
  }

  /* free the allocated RAM */
  delmxInfo(&zinfo);delmxInfo(&xinfo);delmxInfo(&yinfo);
  mxFree(opnmStr); /* free this ram */
}

