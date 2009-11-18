/*

tprod_mex.c

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
  Z=tprod(X,[],[1 -2 -3],[1 2 3]); % per trial data variance over dims
  Z=tprod(X,[],[-1 2 -3],[1 2 3]); % per trial data variance over time
  C=tprod(X,[],[1 2 -3],[2 1 3]);  % per trail dim covariance

The last argument determines wether to use the, fast but perhaps not memory
efficient MATLAB code when possible.  (default is true, set to 0 to only use
memory efficient code)

result has size: [size(X) size(Y)] where the accumulated dims are now size==1

N.B. compile with TPRODLIB defined to get an version without mexFunction
defined to use in other mex-files
  
$Id: tprod_mex.c,v 1.1.1.1 2008/02/27 14:42:55 roboos Exp $


*/

#include "mex.h"
#include "matrix.h"
#include "mxInfo.h"
#include "tprod.h"

#ifndef MAX
#define MAX(A,B)  ((A) > (B) ? (A) : (B))
#define MIN(A,B)  ((A) < (B) ? (A) : (B))
#endif

/* matlab call to do the matrix product when possible */
mxArray* MATLAB_mm(MxInfo zinfo, const MxInfo xinfo, const MxInfo yinfo,
						 const MxInfo xrestinfo, const MxInfo yrestinfo,
						 const MxInfo xmaccinfo, const MxInfo ymaccinfo);

/* size at which it becomes more efficient to use tprod code than calling
	back to matlab */
const int MATLABMACCTHRESHOLDSIZE=201;
const int DEFAULTBLKSZ=32;

/*
  This function computes the generalised matrix products.  Basicially it
  works by spliting the input X and Y matrices into 2 *virtual* sub-matrices,
  one for the "outer" product dimensions (x/y rest) (over which the cartesian
  product is computed) and one for the "inner" product matrices (over which
  the inner product - or multiply accumulate - is computed).  Once these 2
  matrices have been computed the result is simply an outer product in
  tprod and inner product in macc */
void mexFunction(const int nlhs, mxArray *plhs[], 
                 const int nrhs, const mxArray *prhs[]) {
  int i, j, maccnd=0, seqnd=0, BLKSZ=DEFAULTBLKSZ, intype, indtype, retVal;
  const int xarg=0, ydimspecarg=3;
  int xdimspecarg=0, yarg=0; /* possibly other way round? */
  bool useMATLAB=true;
  int callType=0;
  MxInfo xinfo, yinfo, zinfo;
  MxInfo xmaccinfo, ymaccinfo, xrestinfo, yrestinfo;
  int znd, xnidx=0, ynidx=0;
  int *x2yIdx=0;
  if (nrhs < 4 || nrhs >6)	 mexErrMsgTxt("tprod: Incorrect number of inputs.");
  if (nlhs > 1)	 mexErrMsgTxt("tprod: Too many output arguments.");
  /* parse the tprod options list */
  if( nrhs >= 5 && mxIsChar(prhs[4]) ) {
	 char *opnmStr=mxArrayToString(prhs[4]);
	 for (i=0 ; opnmStr[i] != 0; i++ ) {
		switch ( opnmStr[i] ) {
		case 'm': case 'M':   useMATLAB=false;  break;
		case 'n': case 'N':   callType=1;    break; /* new call type */
		case 'o': case 'O':   callType=-1;   break; /* old call type */
		default: 
		  mexWarnMsgTxt("tprod: Unrecognised tprod option");
		}
	 }
	 mxFree(opnmStr); opnmStr=0;
  }

  if ( nrhs==6 && mxGetNumberOfElements(prhs[5])==1 ){
	 BLKSZ=(int)mxGetScalar(prhs[5]);
  }

  /* Get X */
  xinfo = mkmxInfo(prhs[xarg],0);
  /* remove trailing singlenton dimensions from x and y */
  for( i=xinfo.nd; i>1; i-- ) if ( xinfo.sz[i-1]!=1 ) break;  xinfo.nd=i;

  /*----------------------------------------------------------------------*/
  /* Identify the calling convention used */
  if ( callType == 0 ) {
	 /* Try and identify the calling convention used, i.e. 
		 X,Y,xdimspec,ydimspec, or X,xdimspec,Y,ydimspec (deprecated) */
	 if (  mxGetNumberOfDimensions(prhs[1])==2 &&
			 ((mxGetN(prhs[1])==1 && mxGetM(prhs[1])>=xinfo.nd) /* size is OK */
			  || (mxGetN(prhs[1])>=xinfo.nd && mxGetM(prhs[1])==1)) ) {
		callType = 1 ;  /* new call type */
	 } 
	 if( mxGetNumberOfDimensions(prhs[2])==2 &&
		  ((mxGetN(prhs[2])==1 && mxGetM(prhs[2])>=xinfo.nd)/* size OK?*/
			|| (mxGetN(prhs[2])>=xinfo.nd && mxGetM(prhs[2])==1)) ){
		/* argument 3 is CONSISTENT */
		if ( callType == 1 ) { /* argument 2 consistent ALSO */
		  mexWarnMsgTxt("tprod: Could not unambigiously determine calling convention, tprod(X,xdimspec,Y,ydimspec) assumed. Use 'n' or 'o' to force new/old convention.");
		} else {
		  callType = -1; 
		  xdimspecarg=2; yarg=1; /* X,Y,xdimspec,ydimspec call */
		}
	 } 
  }
  switch ( callType ) {
  case 1:  xdimspecarg=1; yarg=2; break;/* new type: X,xdimspec,Y,xdimspec */
  case -1: xdimspecarg=2; yarg=1; break;/* old type: X,Y,xdimspec,xdimspec */
  otherwise: mexErrMsgTxt("tprod: Couldnt identify calling convention.");
  }

  /* Now we know where the Y is we can get it too! */
  yinfo = mkmxInfo(prhs[yarg],0); 
  intype = (xinfo.ip==0?0:1) + (yinfo.ip==0?0:2); /* now compute type call */
  /* empty set as input for y means make it copy of x */
  if ( yinfo.numel==0 && yinfo.nd==2 && yinfo.sz[0]==0 && yinfo.sz[1]==0 ) { 
	 yinfo=copymxInfo(xinfo); 
  }
  /* remove trailing singlenton dimensions from x and y */
  for( i=yinfo.nd; i>1; i-- ) if ( yinfo.sz[i-1]!=1 ) break;  yinfo.nd=i;

  /* check the types are what we can use */
  if ( xinfo.dtype == DOUBLE_DTYPE ){
	 if ( yinfo.dtype == DOUBLE_DTYPE ) { /* dd */
		indtype = 11; 
	 } else if ( yinfo.dtype == SINGLE_DTYPE ) { /* ds */
		indtype = 12;
	 } else {
		mexErrMsgTxt("tprod: Y is of unsupported type: only double");/* or single");*/
	 }
  } else if ( xinfo.dtype == SINGLE_DTYPE ) {
	 if ( yinfo.dtype == DOUBLE_DTYPE ) { /* sd */
		indtype = 21;
	 } else if ( yinfo.dtype == SINGLE_DTYPE ) { /* ss */
		indtype = 22;
	 } else {
		mexErrMsgTxt("tprod: Y is of unsupported type: only double or single");
	 }
  } else {
	 mexErrMsgTxt("tprod: X is of unsuppored type: only double");/* or single");*/
  }

     
  /*-------------------------------------------------------------------------
   * Initialise x2yIdx which maps from input dimensions to the type of output
	* we want.  The format of x2yIdx is:
	*   [Xdimlabels Ydimlabels] 
	* so it reflects the order of the dimensions in the output Z.
   * The value and sign of the integer label tells us what type of
   * operation to perform on this dimension.
	*    0   this is a squeezed dim which must be singlenton
	*   -ve  this is the index in X/Y of the matching dim to be accumulated
	*   +ve  this number is the position of the dimension at this location in
   *        the output matrix.
	*-------------------------------------------------------------------------*/
  znd=0; /* xnidx, ynidx; */
  double *xip, *yip;
  /* fill in the x2yIdx for the new type of indicies */
  maccnd=0;
  xnidx=mxGetNumberOfElements(prhs[xdimspecarg]);
  if( xnidx<xinfo.nd ) mexErrMsgTxt("tprod:Less X indicies than dimensions");
  ynidx=mxGetNumberOfElements(prhs[ydimspecarg]);
  if( ynidx<yinfo.nd ) mexErrMsgTxt("tprod:Less Y indicies than dimensions");

  x2yIdx=(int*)mxCalloc(xnidx+ynidx,sizeof(int));  	 
  xip=mxGetPr(prhs[xdimspecarg]); yip=mxGetPr(prhs[ydimspecarg]);

  /* find the max value of xip, this is num output dims */
  /* also check for non multiple instances of acc dims labels */
  znd=MAX((int)xip[0],0);
  for( i=1; i<xnidx; i++){
	 if ( xip[i] < .0 ) {
		for( j=0; j<i; j++)
		  if(xip[i]==xip[j]) mexErrMsgTxt("tprod: Duplicate x-dim label");
	 } else if ( xip[i] > .0 ) {
		znd=MAX(znd,(int)xip[i]); /* find the largest target dim */
	 } else if ( sz(xinfo,i)!=1 ) 
		mexErrMsgTxt("tprod: Ignored dims *MUST* have size 1");
  }
  /* same for yip */
  /* but also check that we always have matching x label to accumulate */
  znd=MAX(znd,(int)yip[0]); 
  for( i=1; i<ynidx; i++){
	 if ( yip[i] < .0 ) {
		for( j=0; j<i; j++)
		  if(yip[i]==yip[j]) mexErrMsgTxt("tprod: Duplicate y-dim label");
		for( j=0; j<xnidx; j++) if( yip[i]==xip[j] ) break;
		if( j==xnidx ) mexErrMsgTxt("tprod: Couldn't find matching x-idx for the y");
	 } else if ( yip[i] > .0 ) {
		znd=MAX(znd,(int)yip[i]); /* find the largest target dim */
	 } else if ( sz(yinfo,i)!=1 ) 
		mexErrMsgTxt("tprod: Ignored dims *MUST* have size 1");
  }

  /* compute the x->y mapping */
  for( i=0; i<xnidx; i++){
	 if ( xip[i] < .0 ) {
		/* search for the matching y */
		for( j=0; j<ynidx; j++) {
		  if ( yip[j]==xip[i] ) {
			 x2yIdx[i]=-(j+1);    x2yIdx[xnidx+j]=-(i+1);    maccnd++;
			 break;
		  }
		}
		if ( x2yIdx[i]==0 )
		  mexErrMsgTxt("tprod: Couldn't find a matching y-idx for the x");
		if( sz(xinfo,i) != sz(yinfo,j)) /* check sizes match too */
		  mexErrMsgTxt("tprod: Matched dims must have the same size!");
	 } else { /* just copy it through */
		x2yIdx[i]=(int)xip[i];
		/* search for the matching y, & check sizes match */
		for( j=0; j<ynidx && yip[j]!=xip[i]; j++);
		if ( j<ynidx ){ /* sequential dimension */
		  if ( sz(xinfo,i) != sz(yinfo,j) )
			 mexErrMsgTxt("tprod: Matched dims must have the same size!");
		  if ( sz(xinfo,i)!=1 ) { /* only *really* sequ if >1 size strides */
			 seqnd++; 
		  } 
		}
	 }
  }
  /* now set the y parts -- for the non-set non-accumulated dimensions */ 
  for( i=0; i<ynidx; i++){ 
	 if( yip[i] > .0 ) { x2yIdx[i+xnidx]=(int)yip[i]; }
  }
  /*   } */
  znd=MAX(znd,1); /* ensure at least 1 output dim */
  maccnd=MAX(maccnd,1); /* ensure at least 1 macc dim */
  
  /*-------------- ARGH! THIS IS HORRIBLE! -----------------------------*/
  /* swap X and Y if Y should come in the output before X */
  j = znd;
  for( i=0;   i < xnidx; i++ ) if ( x2yIdx[i]>0 ){ j = x2yIdx[i]; break; }/* +idx in X */
  for( i=xnidx+1; i<xnidx+ynidx; i++ ) if ( x2yIdx[i]>0 ) break; /* +idx in Y */
  if ( i<xnidx+ynidx && j>x2yIdx[i] ) { /* Y should come first, so swap X and Y */
	 xmaccinfo = xinfo; xinfo = yinfo; yinfo = xmaccinfo ; /* swap x and Y */
	 /* re-order x2yIdx for the new x/y labels */
	 int *nx2yIdx=(int*)mxCalloc(xnidx+ynidx,sizeof(int)); /* for the new order */
	 for( i=0; i<xnidx; i++) nx2yIdx[ynidx+i]=x2yIdx[i];    /* copy x info to end */
	 for( i=0; i<ynidx; i++) nx2yIdx[i]      =x2yIdx[xnidx+i];/* copy y info to start */
	 mxFree(x2yIdx); x2yIdx = nx2yIdx;
	 j=xnidx; xnidx=ynidx; ynidx=j; 
	 if ( indtype == 12 ) indtype = 21; 
	 else if( indtype == 21 ) indtype = 12; 
  }


  /* compute the mxInfo for the accumulated and rest sub-matrices */  
  xmaccinfo = mkemptymxInfo(maccnd);
  ymaccinfo = mkemptymxInfo(maccnd);
  /* N.B. xrestinfo.sz holds the	size of the combined x and y rest matrix */
  xrestinfo = mkemptymxInfo(znd); 
  yrestinfo = mkemptymxInfo(znd);

  initrestmaccmxInfo(znd, xinfo, yinfo, x2yIdx, xnidx, ynidx,
							&xrestinfo, &yrestinfo, &xmaccinfo, &ymaccinfo);

  /* N.B. otherwise miss chances to  use MATLAB */
  squeezemxInfoPair(&xrestinfo,&yrestinfo);
  squeezemxInfoPair(&xmaccinfo,&ymaccinfo);

  /* compute the size of the output matrix */
  /* this is the size of x followed by the size of y with accumulated dims
	  replaced with size 1 */
  zinfo=initzmxInfo(znd, xinfo, yinfo, x2yIdx, xnidx, ynidx); 
   
  if ( yinfo.numel==0 || xinfo.numel== 0 ) { /* deal with null inputs */
	 mexWarnMsgTxt("tprod: Empty matrix input!");
	 /* return an empty matrix */
	 plhs[0]=mxCreateNumericArray(zinfo.nd,zinfo.sz,mxDOUBLE_CLASS,
											(xinfo.ip==0&&yinfo.ip==0)?mxREAL:mxCOMPLEX);
	 	 
  } else if ( useMATLAB && /* allowed */
		 seqnd==0 &&                     /* no sequential dims */
		 xmaccinfo.nd <= 1 &&            /* at most 1 macc dim */
		 xrestinfo.nd <= 2 &&            /* at most 2 output dims */
		 (xrestinfo.nd<= 1 ||            /* 1 from X */
		  ((xrestinfo.stride[0]==0) || (stride(xrestinfo,1)==0))) && 
 		 (yrestinfo.nd<= 1 ||            /* 1 from Y */
		  ((yrestinfo.stride[0]==0) || (stride(yrestinfo,1)==0))) &&
				  (xmaccinfo.numel*((intype==0)?1:((intype>2)?4:2)) < MATLABMACCTHRESHOLDSIZE ) /* not tooo big! */
				  ){ 
	 /* Phew! we can use matlab! */
	 if ( xrestinfo.stride[0]>0 ) { /* x comes before y in output */
		plhs[0]=MATLAB_mm(zinfo,xinfo,yinfo,xrestinfo,yrestinfo,
								xmaccinfo,ymaccinfo);
	 } else { /* y comes before x in output, reverse order of inputs */
		plhs[0]=MATLAB_mm(zinfo,yinfo,xinfo,yrestinfo,xrestinfo,
								ymaccinfo,xmaccinfo);
	 }
  } else {

	 /* otherwise do it ourselves */
	 /* create the data for the z matrix and set its pointer */
	 plhs[0]=mxCreateNumericArray(zinfo.nd,zinfo.sz,zinfo.dtype,
											(xinfo.ip==0&&yinfo.ip==0)?mxREAL:mxCOMPLEX);
	 zinfo.rp = mxGetPr(plhs[0]);
	 zinfo.ip = mxGetPi(plhs[0]);

	 /* call tprod to do the real work */
	 /* do the (appropriately typed) operation */
	 switch ( indtype ) {

	 case 11: /* dd */	 
		retVal = ddtprod(zinfo,xrestinfo,yrestinfo,xmaccinfo,ymaccinfo,BLKSZ); 
		break;

	 case 12: /* ds */
		retVal = dstprod(zinfo,xrestinfo,yrestinfo,xmaccinfo,ymaccinfo,BLKSZ); 
		break;

	 case 21: /* sd */
		retVal = sdtprod(zinfo,xrestinfo,yrestinfo,xmaccinfo,ymaccinfo,BLKSZ); 
		break;

	 case 22: /* ss */
		retVal = sstprod(zinfo,xrestinfo,yrestinfo,xmaccinfo,ymaccinfo,BLKSZ); 
		break;
	 }
	 /* check for errors */
	 switch ( retVal ) {
	 case ZTYPEMISMATCH : 
		mexErrMsgTxt("tprod: Z is of unsupported type"); break;
	 case XTYPEMISMATCH :
		mexErrMsgTxt("tprod: X is of unsupported type"); break;
	 case YTYPEMISMATCH :
		mexErrMsgTxt("tprod: Y is of unsupported type"); break;
	 case INTYPEMISMATCH :
		mexErrMsgTxt("tprod: input real/complex mix unsupported"); break;
	 case OTHERERROR :
		mexErrMsgTxt("tprod: Something else went wrong, sorry!"); break;
	 default: ;
	 }
  }
  
  /* ensure the output has the size we want */
  if ( plhs[0] ) mxSetDimensions(plhs[0],zinfo.sz,zinfo.nd);
  
  /* free up all the memory we've allocated */
  /* N.B. not clear we need to do this from the matlab web-site should happen
	  automatically */
  delmxInfo(&xinfo);delmxInfo(&yinfo);delmxInfo(&zinfo);
  delmxInfo(&xmaccinfo);  delmxInfo(&ymaccinfo);
  delmxInfo(&xrestinfo);  delmxInfo(&yrestinfo);
  mxFree(x2yIdx);
}

/*-------------------------------------------------------------------------*/
/* the input problem could be reduced to a conventional 2D matrix product..
	so here we use the fast MATLAB version */
mxArray* MATLAB_mm(MxInfo zinfo, const MxInfo xinfo, const MxInfo yinfo,
						 const MxInfo xrestinfo, const MxInfo yrestinfo,
						 const MxInfo xmaccinfo, const MxInfo ymaccinfo){
  mxArray *Xmx, *Ymx, *Zmx, *args[2];
  /* call matlab to do the real work -- more reliable than dgemm */
  /* first create a new matrix with the right size to use in the matlab call*/
  /* create and empty array */
  Xmx = mxCreateDoubleMatrix(0,0,(xinfo.ip==0)?mxREAL:mxCOMPLEX);
  Ymx = mxCreateDoubleMatrix(0,0,(yinfo.ip==0)?mxREAL:mxCOMPLEX);
  /* and populate it with data */
  mxSetPr(Xmx,xinfo.rp); if ( xinfo.ip ) mxSetPi(Xmx,xinfo.ip);
  mxSetPr(Ymx,yinfo.rp); if ( yinfo.ip ) mxSetPi(Ymx,yinfo.ip);

  /* Set trap so errors return here so we can clean up correctly */
  mexSetTrapFlag(1);

  /* now do the calls to matlab to get the result */
  args[0] = Xmx; args[1]= Ymx;
  Zmx=0;  
  /* no accumulated dims -- just outer product, or just inner product */
  if ( xmaccinfo.nd == 1 && xmaccinfo.sz[0]==1 ) { 
	 /* set X as its vector version, implicitly transpose Y and then prod */
	 mxSetM(Xmx,MAX(sz(xrestinfo,0),1)); mxSetN(Xmx,1);
	 mxSetM(Ymx,1);                      mxSetN(Ymx,MAX(sz(yrestinfo,1),1)); 
	 mexCallMATLAB(1, &Zmx, 2, args, "*");

  } else if ( xmaccinfo.stride[0] == 1 && ymaccinfo.stride[0] == 1 ){ 
	 /* | x | == macc/rest * macc/rest*/
	 mxSetM(Xmx,xmaccinfo.sz[0]);      mxSetN(Xmx,xinfo.numel/xmaccinfo.sz[0]);
	 mxSetM(Ymx,ymaccinfo.sz[0]);      mxSetN(Ymx,yinfo.numel/ymaccinfo.sz[0]);
	 
	 /* transpose X and call matlab */
	 if( yinfo.numel+zinfo.numel>xinfo.numel*.75 ){/* cheaper to transpose X */
		mxArray *XmxT;
		mexCallMATLAB(1, &XmxT, 1, &Xmx, ".\'");/* N.B. this copies X!!! */
		args[0]=XmxT;
		mexCallMATLAB(1, &Zmx, 2, args, "*");
		mxDestroyArray(XmxT);

	 } else { /*cheaper to transpose y and z */
		mxArray *YmxT;
		mexCallMATLAB(1, &YmxT, 1, &Ymx, ".\'");/* N.B. this copies Y!!! */
		args[0]=YmxT; args[1]=Xmx;
		mexCallMATLAB(1, &Zmx, 2, args, "*");
		mxDestroyArray(YmxT);
		if( ymaccinfo.sz[0]!=yinfo.numel ) {/* transpose result if necessary */
		  mxArray *ZmxT;
		  mexCallMATLAB(1, &ZmxT, 1, &Zmx, ".\'");/* N.B. this copies Z!!! */
		  mxDestroyArray(Zmx);
		  Zmx=ZmxT;		
		}
	 }
  } else if ( xmaccinfo.stride[0] == 1 && ymaccinfo.stride[0] > 1 ){
	 /* | x _ == macc/rest * rest/macc */
	 mxSetM(Xmx,xmaccinfo.sz[0]);      
	 mxSetN(Xmx,xinfo.numel/xmaccinfo.sz[0]);
	 mxSetM(Ymx,ymaccinfo.stride[0]);  
	 mxSetN(Ymx,yinfo.numel/ymaccinfo.stride[0]);	 

	 if( yinfo.numel+xinfo.numel<zinfo.numel*.75 ){/* cheaper to transpose Z */
		/* reverse order of multiply, call matlab and transpose Z */	 
		mxArray *ZmxT;
		args[0]=Ymx; args[1]=Xmx;
		mexCallMATLAB(1, &Zmx, 2, args, "*");
		mexCallMATLAB(1, &ZmxT, 1, &Zmx, ".\'");/* N.B. this copies Z!!! */
		mxDestroyArray(Zmx);
		Zmx=ZmxT;

	 } else { /* cheaper to transpose twice */
		mxArray *XmxT, *YmxT;
		mexCallMATLAB(1, &XmxT, 1, &Xmx, ".\'");/* N.B. this copies X!!! */
		mexCallMATLAB(1, &YmxT, 1, &Ymx, ".\'");/* N.B. this copies Y!!! */
		args[0]=XmxT; args[1]=YmxT;
		mexCallMATLAB(1, &Zmx, 2, args, "*");
		mxDestroyArray(XmxT); mxDestroyArray(YmxT);
	 }
	 
  } else if ( xmaccinfo.stride[0] >  1 && ymaccinfo.stride[0] == 1 ){
	 /* _ x | == rest/macc * macc/rest */
	 mxSetM(Xmx,xmaccinfo.stride[0]); 
	 mxSetN(Xmx,xinfo.numel/xmaccinfo.stride[0]);
	 mxSetM(Ymx,ymaccinfo.sz[0]);     
	 mxSetN(Ymx,yinfo.numel/ymaccinfo.sz[0]); 
	 mexCallMATLAB(1, &Zmx, 2, args, "*");

  } else if ( xmaccinfo.stride[0] > 1 && ymaccinfo.stride[0] > 1 ){
	 /* _ x _ == rest/macc * rest/macc */ 
	 mxSetM(Xmx,xmaccinfo.stride[0]);  
	 mxSetN(Xmx,xinfo.numel/xmaccinfo.stride[0]);
	 mxSetM(Ymx,ymaccinfo.stride[0]);  
	 mxSetN(Ymx,yinfo.numel/ymaccinfo.stride[0]);	 

	 if( xinfo.numel+zinfo.numel>yinfo.numel*.75 ){/* cheaper to transpose Y */
		/* transpose Y and call matlab */
		mxArray *YmxT;
		mexCallMATLAB(1, &YmxT, 1, &Ymx, ".\'");/* N.B. this copies Y!!! */
		args[1]=YmxT;
		mexCallMATLAB(1, &Zmx, 2, args, "*");
		mxDestroyArray(YmxT);

	 } else { /* cheaper to transpose X and Z */
		mxArray *XmxT;
		mexCallMATLAB(1, &XmxT, 1, &Xmx, ".\'");/* N.B. this copies X!!! */
		args[0]=Ymx; args[1]=XmxT;
		mexCallMATLAB(1, &Zmx, 2, args, "*");
		mxDestroyArray(XmxT);
		if( xmaccinfo.sz[0]!=xinfo.numel ) {/* transpose result if necessary */
		  mxArray *ZmxT;
		  mexCallMATLAB(1, &ZmxT, 1, &Zmx, ".\'");/* N.B. this copies Z!!! */
		  mxDestroyArray(Zmx);
		  Zmx=ZmxT;		
		}
		
	 }

  } else {
	 mexErrMsgTxt("tprod: somethings gone horibbly wrong!");
  }

  /* Set trap so errors return here so we can clean up correctly */
  mexSetTrapFlag(0);

  /* set the tempory X and Y matrices back to empty without ref to data to
	  stop matlab "helpfully" double freeing them? */
  mxSetM(Xmx,0);mxSetN(Xmx,0);mxSetPr(Xmx,0);mxSetPi(Xmx,0);
  mxDestroyArray(Xmx);
  mxSetM(Ymx,0);mxSetN(Ymx,0);mxSetPr(Ymx,0);mxSetPi(Ymx,0);
  mxDestroyArray(Ymx);
  return Zmx;
}
