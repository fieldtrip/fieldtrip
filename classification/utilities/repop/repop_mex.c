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

 $Id: repop_mex.c,v 1.1.1.1 2008/02/27 14:42:55 roboos Exp $

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
#include "string.h"
#include "repop.h"

/*---------------------------------------------------------------------------*
 * use the macro-preprocessor and the repops.template file to automatically *
 * generate the functions to call for the different combinations of real and
 * complex inputs */

#if defined(_PLUS_) ||defined(_MINUS_) || defined(_TIMES_) || defined(_POWER_)\
    || defined(_EQ_) || defined(_NE_) || defined(_LT_) || defined(_GT_) \
    || defined(_LE_) || defined(_GE_)
#error Only define the operator to use for repop_ind.c
#endif

/*-------------------------------------------------------------------------*/
/* the wrapper function */
void mexFunction(const int nlhs, mxArray *plhs[],
                 const int nrhs, const mxArray *prhs[]) {
  RepopErrorCode retVal=OK ;
  char *opnmStr;
  int xnd, ynd, znd, opid, xargi=0, yargi=1, opargi=2, i=0;
  int intype, indtype; /* real/complex input mix, and input datatype mix */
  MxInfo xinfo, yinfo, zinfo;
  int repNonUnitDim=0, inPlace=false;
                     
  if (nrhs < 3) mexErrMsgTxt("repop: Incorrect number of inputs.");        
  if (nlhs > 1) mexErrMsgTxt("repop: Too many output arguments.");        

  if( nrhs>3 ){ /* parse the options */
	 if ( mxIsChar(prhs[3]) ){
		opnmStr=mxArrayToString(prhs[3]);
		for (i=0 ; opnmStr[i] != 0; i++ ) {
		  switch ( opnmStr[i] ) {
		  case 'i': case 'I':   inPlace=1;  break;
		  case 'm': case 'M':   repNonUnitDim=1;  break;
		  case 'n': case 'N':   repNonUnitDim=2;  break;
		  default: 
			 mexWarnMsgTxt("repop: Unrecognised tprod option");
		  }
		}
		mxFree(opnmStr);
	 }
  }

  /* first -- find the operator position and decipher which operator we're
	  supposed to use */
  if( mxIsChar(prhs[0])) { xargi=1; yargi=2; opargi=0; }
  if( mxIsChar(prhs[1])) { yargi=2; opargi=1; }
  if(!mxIsChar(prhs[opargi])) mexErrMsgTxt("repop: one arg must be operator type"); 
  opnmStr=mxArrayToString(prhs[opargi]);
  opid = getOpid(opnmStr);  
  mxFree(opnmStr); /* free this ram */
  if ( opid < 0 || opid > GE ) 
	 mexErrMsgTxt("repop:Unrecognised operator. Must be one of: + - * / \\ = < > ~= <= >=");

  xnd = mxGetNumberOfDimensions(prhs[xargi]);
  ynd = mxGetNumberOfDimensions(prhs[yargi]);
  znd = max(xnd,ynd);
  xinfo=mkmxInfo(prhs[xargi],znd);
  yinfo=mkmxInfo(prhs[yargi],znd);

  /* empty set as input for y means make it copy of x */
  if ( yinfo.numel==0 && yinfo.nd==2 && yinfo.sz[0]==0 && yinfo.sz[1]==0 ) { 
	 yinfo=copymxInfo(xinfo); 
  }

  /* type info of the inputs */
  intype= (xinfo.ip==0?0:1) + (yinfo.ip==0?0:2);/*%1=x complex %2=y complex*/
  if ( xinfo.dtype == DOUBLE_DTYPE ){
	 if ( yinfo.dtype == DOUBLE_DTYPE ) { /* dd */
		indtype = 11; 
	 } else if ( yinfo.dtype == SINGLE_DTYPE ) { /* ds */
		indtype = 12; 
	 } else {
		mexErrMsgTxt("repop: Y is of unsupported type: only double or single");
	 }
  } else if ( xinfo.dtype == SINGLE_DTYPE ) {
	 if ( yinfo.dtype == DOUBLE_DTYPE ) { /* sd */
		indtype = 21; 
	 } else if ( yinfo.dtype == SINGLE_DTYPE ) { /* ss */
		indtype = 22; 
	 } else { 
		mexErrMsgTxt("repop: Y is of unsupported type: only double or single");
	 }
  } else {
	 mexErrMsgTxt("repop: X is of unsuppored type: only double or single");
  }

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
		if ( (intype&1)==0 && ( intype!=0 || opid==POWER ) ) inPlace=0; 
		if ( indtype==12 ) inPlace = 0 ; /* not if single, double call */
		if ( !inPlace ) 
		  mexWarnMsgTxt("repop: InPlace failure: size(Z)~=size(X) | isreal(Z)~=isreal(X) | isa(Z)~=isa(X)");
	 }

	 /* allocate the space for z */
	 if( opid >= EQ ) {/*logical operator so make logical output array*/
		plhs[0]  = mxCreateLogicalArray(zinfo.nd, zinfo.sz);
		zinfo.rp = (double*)mxGetLogicals(plhs[0]);
	 } else { /* not a logical operator so Z is same type as X */
		if ( opid == POWER || intype != 0 ) {
		  zinfo.ip = (double*)1 ; /* set ip to indicate need complex array */
		}
		if ( inPlace ) { /* use X inplace */
		  /* create a empty temp array & set X to point to it  */
		  plhs[0] = mxCreateNumericMatrix(0,0,zinfo.dtype,
													 zinfo.ip?mxCOMPLEX:mxREAL);
		  mxSetPr((mxArray*)prhs[xargi],mxGetPr(plhs[0])); 
		  mxSetPi((mxArray*)prhs[xargi],mxGetPi(plhs[0]));
		  mxSetDimensions((mxArray*)prhs[xargi],mxGetDimensions(plhs[0]),2);
		  /* set Z to use Xs old memory instead */
		  mxSetPr(plhs[0],xinfo.rp);		  mxSetPi(plhs[0],xinfo.ip);
		  mxSetDimensions(plhs[0],zinfo.sz,zinfo.nd);
		} else { /* allocate the output array info+mem */
		  plhs[0] = mxCreateNumericArray(zinfo.nd, zinfo.sz, zinfo.dtype, 
													zinfo.ip?mxCOMPLEX:mxREAL);
		}
		zinfo.rp = mxGetPr(plhs[0]);  
		zinfo.ip = mxGetPi(plhs[0]);
	 }
	 /* do the (appropriately typed) operation */
	 switch ( indtype ) {
	 case 11: /* dd */
		retVal = ddrepop( zinfo,  xinfo,  yinfo, opid ); break;
	 case 12: /* ds */
		retVal = dsrepop( zinfo,  xinfo,  yinfo, opid ); break;
	 case 21: /* sd */
		retVal = sdrepop( zinfo,  xinfo,  yinfo, opid ); break;
	 case 22: /* ss */
		retVal = ssrepop( zinfo,  xinfo,  yinfo, opid ); break;
	 }

	 /* check for errors */
	 switch ( retVal ) {
	 case ZTYPEMISMATCH : 
		mexErrMsgTxt("repop: Z is of unsupported type"); break;
	 case XTYPEMISMATCH :
		mexErrMsgTxt("repop: X is of unsupported type"); break;
	 case YTYPEMISMATCH :
		mexErrMsgTxt("repop: Y is of unsupported type"); break;
	 case INTYPEMISMATCH :
		mexErrMsgTxt("repop: input real/complex mix unsupported"); break;
	 case OTHERERROR :
		mexErrMsgTxt("repop: Something else went wrong, sorry!"); break;
	 default: ;
	 }
  }

  /* free the allocated RAM */
  delmxInfo(&zinfo);delmxInfo(&xinfo);delmxInfo(&yinfo);
}

