/*

tprod_mex.c

Generalised multiply accumulate operation. == tensor product with repeated
indicies

TODO:
* add a slicing ability to the input.  Format:
  [2 x d] start and end points for indexing in each col. 0 means all indices
    or 
  [3 x d] start,end,step points for indexing in each col, 0 means all indices
  
$Id$

*/

#include "mex.h"
#include "matrix.h"
#include "mxInfo.h"
#include "mxInfo_mex.h"
#include "tprod.h"

#ifndef MAX
#define MAX(A,B)  ((A) > (B) ? (A) : (B))
#define MIN(A,B)  ((A) < (B) ? (A) : (B))
#endif


/* matlab call to do the matrix product when possible */
mxArray* MATLAB_mm(MxInfo zinfo, const MxInfo xinfo, const MxInfo yinfo,
						 const MxInfo xrest, const MxInfo yrest,
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
  int i, maccnd=0, seqnd=0, BLKSZ=DEFAULTBLKSZ, err=OK; 
  const int xargi=0, ydimspecarg=3;
  int xdimspecarg=0, yargi=0; /* possibly other way round? */
  bool useMATLAB=true;
  int callType=0;
  MxInfo xinfo, yinfo, zinfo;
  MxInfo xmacc, ymacc, zrest, xrest, yrest;
  int znd, xnidx=0, ynidx=0;
  int *x2yIdx=0;
  if (nrhs < 4 || nrhs >6)	 { ERROR("tprod: Incorrect number of inputs."); err=OTHERERROR; }
  if (nlhs > 1) { ERROR("tprod: Too many output arguments."); err=OTHERERROR; }
  if ( mxGetNumberOfDimensions(prhs[ydimspecarg]) > 2 ){
	 ERROR("tprod: ydimspec must be a vector"); err=OTHERERROR; }
  /* parse the tprod options list */
  if( nrhs >= 5 && mxIsChar(prhs[4]) ) {
	 char *opnmStr=mxArrayToString(prhs[4]);
	 for (i=0 ; opnmStr[i] != 0; i++ ) {
		switch ( opnmStr[i] ) {
		case 'm': case 'M':   useMATLAB=false;  break;
		case 'n': case 'N':   callType=1;    break; /* new call type */
		case 'o': case 'O':   callType=-1;   break; /* old call type */
		default: 
		  WARNING("tprod: Unrecognised tprod option");
		}
	 }
	 mxFree(opnmStr); opnmStr=0;
  }

  if ( nrhs==6 && mxGetNumberOfElements(prhs[5])==1 ){
	 BLKSZ=(int)mxGetScalar(prhs[5]);
  }

  /* Get X */
  xinfo = mkmxInfoMxArray(prhs[xargi],0);
  /* remove trailing singlenton dimensions from x and y */
  for( i=xinfo.nd; i>1; i-- ) if ( xinfo.sz[i-1]!=1 ) break;  xinfo.nd=i;

  /*----------------------------------------------------------------------*/
  /* Identify the calling convention used */
  if ( callType == 0 ) {

	 /* Try and identify the calling convention used, i.e. 
		 X,Y,xdimspec,ydimspec, or X,xdimspec,Y,ydimspec (deprecated) */
	 if ( mxGetNumberOfDimensions(prhs[1])==2 &&		
			( (mxGetN(prhs[1])==1 && mxGetM(prhs[1])>=xinfo.nd) /* Xdimspec OK */
			  || (mxGetN(prhs[1])>=xinfo.nd && mxGetM(prhs[1])==1) ) ) {
		int ynd = 0; 
		const int *ysz = 0;
		yargi = 2; /* start by trying new call type */
		ynd=mxGetNumberOfDimensions(prhs[yargi]);
		ysz = mxGetDimensions(prhs[yargi]); /* size of poss Y*/
		for( i=ynd; i>1; i-- ) if ( ysz[i-1]!=1 ) break; ynd=i; 
		if( mxGetNumberOfElements(prhs[ydimspecarg]) >= ynd ){/*Ydimspec OK*/
		  callType = 1 ;  /* new call type */
		}
	 }
	 if( mxGetNumberOfDimensions(prhs[2])==2 &&
		  ((mxGetN(prhs[2])==1 && mxGetM(prhs[2])>=xinfo.nd)/* xdimspec OK */
			|| (mxGetN(prhs[2])>=xinfo.nd && mxGetM(prhs[2])==1) ) ) {
		/* Consitent so far, check the ydimspec is OK */
		int ynd = 0; 
		const int *ysz = 0;
		yargi  = 1; /* start by trying new call type */
		ynd=mxGetNumberOfDimensions(prhs[yargi]);
		ysz = mxGetDimensions(prhs[yargi]); /* size of poss Y*/
		for( i=ynd; i>1; i-- ) if ( ysz[i-1]!=1 ) break; ynd=i; 
		if ( mxGetNumberOfElements(prhs[ydimspecarg]) >= ynd ) {
		  
		  if ( callType == 0 ) { /* argument 3 is CONSISTENT, and 2 *WASN'T* */
			 callType = -1; 
			 
		  } else { /* argument 2 consistent ALSO */
			 /* of one of them matches exactly and the other doesn't */
			 int xnd = mxGetNumberOfDimensions(prhs[xargi]); /* num input X dims */
			 if ( xnd==2 && xinfo.nd == 1 ) xnd=1; /* col vec is special case*/
			 if ( xnd == mxGetNumberOfElements(prhs[1]) /* 1 matches *exactly* */
					&& xnd != mxGetNumberOfElements(prhs[2]) ) { /* 2 doesn't */
				callType = 1;
			 } else if( xnd == mxGetNumberOfElements(prhs[2])/* 2 *exact* match */
							&& xnd != mxGetNumberOfElements(prhs[1]) ){/* 1 doesn't */
				callType = -1;
			 } else { /* neither or both match exactly */
				callType = 1;
				WARNING("tprod: Could not unambigiously determine calling convention, tprod(X,xdimspec,Y,ydimspec) assumed. Use 'n' or 'o' to force new/old convention.");
			 }
		  }
		}
	 } 
  }
  switch ( callType ) {
  case 1:  xdimspecarg=1; yargi=2; break;/* new type: X,xdimspec,Y,xdimspec */
  case -1: xdimspecarg=2; yargi=1; break;/* old type: X,Y,xdimspec,xdimspec */
  default: ERROR("tprod: Couldnt identify calling convention."); err=OTHERERROR;
  }

  /* Now we know where the Y is we can get it too! */
  yinfo = mkmxInfoMxArray(prhs[yargi],0); 
  /* empty set as input for y means make it copy of x */
  if ( yinfo.numel==0 && yinfo.nd==2 && yinfo.sz[0]==0 && yinfo.sz[1]==0 ) { 
	 yinfo=copymxInfo(xinfo); 
  }
  /* remove trailing singlenton dimensions from x and y */
  for( i=yinfo.nd; i>1; i-- ) if ( yinfo.sz[i-1]!=1 ) break;  yinfo.nd=i;

  /* check the types are what we can use */
  if ( !(xinfo.dtype == DOUBLE_DTYPE || xinfo.dtype == SINGLE_DTYPE ) ){
	 ERROR("tprod: X type unsuppored: only full double/single"); err=UNSUPPORTEDINPUTS;
  }
  if ( mxIsSparse(prhs[xargi]) ){
	 ERROR("tprod: X is sparse, only full double/single supported"); err=UNSUPPORTEDINPUTS;
  }
  if ( !(yinfo.dtype == DOUBLE_DTYPE || yinfo.dtype == SINGLE_DTYPE ) ){
	 ERROR("tprod: Y type unsuppored: only double, single"); err=UNSUPPORTEDINPUTS;
  }
  if ( mxIsSparse(prhs[yargi]) ){
	 ERROR("tprod: X is sparse, only full double/single supported"); err=UNSUPPORTEDINPUTS;
  }
       
  /* fill in the x2yIdx for the new type of indicies */
  maccnd=0;
  xnidx=mxGetNumberOfElements(prhs[xdimspecarg]);
  ynidx=mxGetNumberOfElements(prhs[ydimspecarg]);
  err = compx2yIdx_dd(xinfo,xnidx,mxGetPr(prhs[xdimspecarg]), 
								 yinfo,ynidx,mxGetPr(prhs[ydimspecarg]),
								 &x2yIdx,&znd,&maccnd,&seqnd);
  if ( err != 0 ) { /* FREE and return */
	 delmxInfo(&xinfo);  delmxInfo(&yinfo);
	 if ( x2yIdx != 0 ) FREE(x2yIdx); 
	 return ; 
  }

  #ifdef LOGGING  
  FILE *fd = fopen("/tmp/tprod.log","a");
  if( fd==0 ){WARNING("Couldn't open log file /tmp/tprod.log\n"); fd=stderr;}
  fprintf(fd,"tprod( ");
  printMxInfoSummary(fd,xinfo);fprintf(fd," ,[");
  for(i=0; i<xnidx; i++) fprintf(fd,"%d ", (int)mxGetPr(prhs[xdimspecarg])[i]);
  fprintf(fd,"] , ");
  printMxInfoSummary(fd,yinfo);fprintf(fd," ,[");  
  for(i=0; i<ynidx; i++) fprintf(fd,"%d ", (int)mxGetPr(prhs[ydimspecarg])[i]);
  fprintf(fd,"] )\n"); 
  if ( fd!=stderr ) fclose(fd);
  #endif

  /* compute the mxInfo for the accumulated and rest sub-matrices */  
  xmacc = mkemptymxInfo(maccnd);
  ymacc = mkemptymxInfo(maccnd);
  /* N.B. xrest.sz holds the	size of the combined x and y rest matrix */
  xrest = mkemptymxInfo(znd); 
  yrest = mkemptymxInfo(znd);

  err = initrestmaccmxInfo(znd, xinfo, yinfo, 
										x2yIdx, xnidx, ynidx,
										&xrest, &yrest, &xmacc, &ymacc);
  if ( err != 0 ) { /* FREE and return */
	 delmxInfo(&xinfo);  delmxInfo(&yinfo);
	 delmxInfo(&xmacc);  delmxInfo(&ymacc);
	 delmxInfo(&xrest);  delmxInfo(&yrest);
	 FREE(x2yIdx);
	 return; 
  }

  /* compute the size of the output matrix */
  zinfo=initzmxInfo(znd, xinfo, yinfo, x2yIdx, xnidx, ynidx); 

  /* optimise/standardize the query so its the way tprod wants it */
  zrest = copymxInfo(zinfo);
  err =  optimisetprodQuery(&zrest, &xrest, &yrest, 
										 &xmacc, &ymacc);  
  if ( err != OK ) { /* free everything and return */
	 delmxInfo(&xinfo);  delmxInfo(&yinfo); delmxInfo(&zinfo);
	 delmxInfo(&xmacc);  delmxInfo(&ymacc);
	 delmxInfo(&xrest);  delmxInfo(&yrest); delmxInfo(&zrest);
	 FREE(x2yIdx);
	 return; 
  }
  
  if ( xrest.rp != xinfo.rp ) { /* swap xinfo/yinfo if needed */
	 MxInfo tmp = xinfo; xinfo=yinfo; yinfo=tmp;
  }

  /* Now do the actuall work to compute the result */  
  if ( yinfo.numel==0 || xinfo.numel== 0 ) { /* deal with null inputs */
	 WARNING("tprod: Empty matrix input!");
	 /* return an empty matrix */
	 plhs[0]=mxCreateNumericArray(zinfo.nd,zinfo.sz,mxDOUBLE_CLASS,
											(xinfo.ip==0&&yinfo.ip==0)?mxREAL:mxCOMPLEX);
	 	 
  } else if ( useMATLAB && /* allowed */
		 seqnd==0 &&                     /* no sequential dims */
		 xmacc.nd <= 1 &&            /* at most 1 macc dim */
		 xrest.nd <= 2 &&            /* at most 2 output dims */
		 (xrest.nd<= 1 ||            /* 1 from X */
		  ((xrest.stride[0]==0) || (stride(xrest,1)==0))) && 
 		 (yrest.nd<= 1 ||            /* 1 from Y */
		  ((yrest.stride[0]==0) || (stride(yrest,1)==0))) &&
				  (xmacc.numel*(2+(xinfo.ip==0)+(yinfo.ip==0)) < MATLABMACCTHRESHOLDSIZE ) /* not tooo big! */
				  ){ 
	 /* Phew! we can use matlab! */
	 if ( xrest.stride[0]>0 ) { /* x comes before y in output */
		plhs[0]=MATLAB_mm(zinfo,xinfo,yinfo,xrest,yrest,
								xmacc,ymacc);
	 } else { /* y comes before x in output, reverse order of inputs */
		plhs[0]=MATLAB_mm(zinfo,yinfo,xinfo,yrest,xrest,
								ymacc,xmacc);
	 }
  } else {

	 /* otherwise do it ourselves */
	 /* create the data for the z matrix and set its pointer */
	 plhs[0]=mxCreateNumericArray(zinfo.nd,zinfo.sz,zinfo.dtype,
											(xinfo.ip==0&&yinfo.ip==0)?mxREAL:mxCOMPLEX);
	 zinfo.rp = mxGetPr(plhs[0]); zrest.rp=zinfo.rp;
	 zinfo.ip = mxGetPi(plhs[0]); zrest.ip=zinfo.ip;

	 /* call tprod to do the real work */
	 /* do the (appropriately typed) operation */
	 if(        xrest.dtype==DOUBLE_DTYPE && yrest.dtype==DOUBLE_DTYPE ) {/*dd*/
		err= ddtprod(zrest,xrest,yrest,xmacc,ymacc,BLKSZ);

	 } else if( xrest.dtype==DOUBLE_DTYPE && yrest.dtype==SINGLE_DTYPE ) {/*ds*/
		err= dstprod(zrest,xrest,yrest,xmacc,ymacc,BLKSZ);

	 } else if( xrest.dtype==SINGLE_DTYPE && yrest.dtype==DOUBLE_DTYPE ) {/*sd*/
		err= sdtprod(zrest,xrest,yrest,xmacc,ymacc,BLKSZ);
 
	 } else if( xrest.dtype==SINGLE_DTYPE && yrest.dtype==SINGLE_DTYPE ){/*ss*/
		err= sstprod(zrest,xrest,yrest,xmacc,ymacc,BLKSZ);

	 } else {
		err= UNSUPPORTEDINPUTS;
		
	 }
	 /* check for errors */
	 switch ( err ) {
	 case ZTYPEMISMATCH : 
		ERROR("tprod: Z is of unsupported type"); break;
	 case XTYPEMISMATCH :
		ERROR("tprod: X is of unsupported type"); break;
	 case YTYPEMISMATCH :
		ERROR("tprod: Y is of unsupported type"); break;
	 case INTYPEMISMATCH :
		ERROR("tprod: input real/complex mix unsupported"); break;
	 case OTHERERROR :
		ERROR("tprod: Something else went wrong, sorry!"); break;
	 case UNSUPPORTEDINPUTS :
 		ERROR("tprod: Inputs of unsupported type: only double/single"); break;
	 default: ;
	 }
  }
  
  /* ensure the output has the size we want */
  if ( err==OK && plhs[0] ) mxSetDimensions(plhs[0],zinfo.sz,zinfo.nd);
  
  /* free up all the memory we've allocated */
  /* N.B. not clear we need to do this from the matlab web-site should happen
	  automatically */
  delmxInfo(&xinfo);delmxInfo(&yinfo);delmxInfo(&zinfo);
  delmxInfo(&xmacc);  delmxInfo(&ymacc);
  delmxInfo(&xrest);  delmxInfo(&yrest);
  FREE(x2yIdx);
}

/*-------------------------------------------------------------------------*/
/* the input problem could be reduced to a conventional 2D matrix product..
	so here we use the fast MATLAB version */
mxArray* MATLAB_mm(MxInfo zinfo, const MxInfo xinfo, const MxInfo yinfo,
						 const MxInfo xrest, const MxInfo yrest,
						 const MxInfo xmacc, const MxInfo ymacc){
  mxArray *Xmx, *Ymx, *Zmx, *args[2];
  /* call matlab to do the real work -- more reliable than dgemm */
  /* first create a new matrix with the right size to use in the matlab call*/
  /* create and empty array */
  Xmx= mxCreateNumericMatrix(0,0,xinfo.dtype,(xinfo.ip==0)?mxREAL:mxCOMPLEX);
  Ymx= mxCreateNumericMatrix(0,0,yinfo.dtype,(yinfo.ip==0)?mxREAL:mxCOMPLEX);
  /* and populate it with data */
  mxSetPr(Xmx,xinfo.rp); if ( xinfo.ip ) mxSetPi(Xmx,xinfo.ip);
  mxSetPr(Ymx,yinfo.rp); if ( yinfo.ip ) mxSetPi(Ymx,yinfo.ip);

  /* Set trap so errors return here so we can clean up correctly */
  mexSetTrapFlag(1);

  /* now do the calls to matlab to get the result */
  args[0] = Xmx; args[1]= Ymx;
  Zmx=0;  
  /* no accumulated dims -- just outer product, or just inner product */
  if ( xmacc.nd == 1 && xmacc.sz[0]==1 ) { 
	 /* set X as its vector version, implicitly transpose Y and then prod */
	 mxSetM(Xmx,MAX(sz(xrest,0),1)); mxSetN(Xmx,1);
	 mxSetM(Ymx,1);                  mxSetN(Ymx,MAX(sz(yrest,1),1)); 
	 mexCallMATLAB(1, &Zmx, 2, args, "*");

  } else if ( xmacc.stride[0] == 1 && ymacc.stride[0] == 1 ){ 
	 /* | x | == macc/rest * macc/rest*/
	 mxSetM(Xmx,xmacc.sz[0]);     mxSetN(Xmx,xinfo.numel/xmacc.sz[0]);
	 mxSetM(Ymx,ymacc.sz[0]);     mxSetN(Ymx,yinfo.numel/ymacc.sz[0]);
	 
	 /* transpose X and call matlab */
	 if( yinfo.numel+zinfo.numel>xinfo.numel*.75 ){/*cheaper to transpose X*/
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
		if( ymacc.sz[0]!=yinfo.numel ) {/* transpose result if necessary */
		  mxArray *ZmxT;
		  mexCallMATLAB(1, &ZmxT, 1, &Zmx, ".\'");/* N.B. this copies Z!!! */
		  mxDestroyArray(Zmx);
		  Zmx=ZmxT;		
		}
	 }
  } else if ( xmacc.stride[0] == 1 && ymacc.stride[0] > 1 ){
	 /* | x _ == macc/rest * rest/macc */
	 mxSetM(Xmx,xmacc.sz[0]);      
	 mxSetN(Xmx,xinfo.numel/xmacc.sz[0]);
	 mxSetM(Ymx,ymacc.stride[0]);  
	 mxSetN(Ymx,yinfo.numel/ymacc.stride[0]);	 

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
	 
  } else if ( xmacc.stride[0] >  1 && ymacc.stride[0] == 1 ){
	 /* _ x | == rest/macc * macc/rest */
	 mxSetM(Xmx,xmacc.stride[0]); 
	 mxSetN(Xmx,xinfo.numel/xmacc.stride[0]);
	 mxSetM(Ymx,ymacc.sz[0]);     
	 mxSetN(Ymx,yinfo.numel/ymacc.sz[0]); 
	 mexCallMATLAB(1, &Zmx, 2, args, "*");

  } else if ( xmacc.stride[0] > 1 && ymacc.stride[0] > 1 ){
	 /* _ x _ == rest/macc * rest/macc */ 
	 mxSetM(Xmx,xmacc.stride[0]);  
	 mxSetN(Xmx,xinfo.numel/xmacc.stride[0]);
	 mxSetM(Ymx,ymacc.stride[0]);  
	 mxSetN(Ymx,yinfo.numel/ymacc.stride[0]);	 

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
		if( xmacc.sz[0]!=xinfo.numel ) {/* transpose result if necessary */
		  mxArray *ZmxT;
		  mexCallMATLAB(1, &ZmxT, 1, &Zmx, ".\'");/* N.B. this copies Z!!! */
		  mxDestroyArray(Zmx);
		  Zmx=ZmxT;		
		}
		
	 }

  } else {
	 ERROR("tprod: somethings gone horibbly wrong!");
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
