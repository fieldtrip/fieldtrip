/* 
	File:
	
	Function to do tprod calls without the matlab dependencies.  Used for
	memory debugging or as an example for integration of tprod with other
	numeric libraries.

	Calling Convention:
	tprod_testcases xtypeSpec xsize xdataSrc ytypeSpec ysize ydataSrc xdimSpec ydimSpec

	Examples
	tprod_testcases dr 2x1 "1 2" dr 2x1 "1 2" "1 -2" "1 -2"


*/

#include <stdlib.h>
#include <stdio.h>
#include "mxInfo.h"
#include "tprod.h"
#include "mxUtils.h"

/* first provide necessary call-back functions */
void *CALLOC(size_t nmemb, size_t size){  return (void*)calloc(nmemb,size); }
void *MALLOC(size_t size) { return (void*)malloc(size); }
void FREE(void *ptr) {  free(ptr); }
void ERROR(const char *msg) { fprintf(stderr,"%s\n",msg); }
void WARNING(const char *msg) { fprintf(stderr,"%s\n",msg); }

static int BLKSZ=32;

int main(int argc, char **argv){
  char *xtypestr, *ytypestr, *xszstr, *yszstr, *xdatSrcStr, *ydatSrcStr, *xdimSpecStr, *ydimSpecStr, *blkszstr;
  MxInfoDTypes xdtype, ydtype;
  int xreal, yreal;
  int xsz[5], ysz[5];
  int xnd=5, xnumel, ynd=5, ynumel;
  int err=0;
  void *xrp, *xip, *yrp, *yip;
  MxInfo xinfo, yinfo;
  MxInfo zinfo;
  MxInfo xmacc, ymacc, zrest, xrest, yrest;
  int *x2yIdx=0, xnidx=5, xidx[5], ynidx=5, yidx[5];
  int znd=0,maccnd=0,seqnd=0;

  if ( argc < 8 ) { ERROR("Need at least 4 inputs."); exit(1);   }

  /* Get the different parts of the input arguments */
  xtypestr=argv[1]; xszstr=argv[2]; xdatSrcStr=argv[3]; xdimSpecStr=argv[4];
  ytypestr=argv[5]; yszstr=argv[6]; ydatSrcStr=argv[7]; ydimSpecStr=argv[8];
  blkszstr=argv[9];
  
  if ( argc > 9 ) sscanf(blkszstr,"%d",&BLKSZ);

  /* extract the requested size info from the inputs */
  xnd=parseInfo(xnd,xtypestr,xszstr,&xdtype,&xreal,xsz,&xnumel);
  ynd=parseInfo(ynd,ytypestr,yszstr,&ydtype,&yreal,ysz,&ynumel);

/*   fprintf(stdout,"X =");printInfo(stdout,xdtype,xreal,xsz);printf("\n");*/
/*   fprintf(stdout,"Y =");printInfo(stdout,ydtype,yreal,ysz);printf("\n");*/
  
  /* setup a mxInfo for each of these and allocate the required memory */
  xrp = memalloc(xnumel,xdtype);
  if ( !xreal ) xip = memalloc(xnumel,xdtype); else xip=0;
  yrp = memalloc(ynumel,ydtype);
  if ( !yreal ) yip = memalloc(ynumel,ydtype); else yip=0;

  /* Initialise the values */
  datFill(xdatSrcStr,xnumel,xrp,xip,xdtype);
  datFill(ydatSrcStr,ynumel,yrp,yip,ydtype);

  /* init a mxInfo equivalent */
  xinfo = mkmxInfo(xnd,xsz,xrp,xip,xdtype);
  yinfo = mkmxInfo(ynd,ysz,yrp,yip,ydtype);

  fprintf(stdout,"Xinfo =");printMxInfo(stdout,xinfo);printf("\n");
  fprintf(stdout,"Yinfo =");printMxInfo(stdout,yinfo);printf("\n");


  /* parse the spec */
  readMx(&xnidx, xdimSpecStr, xidx, INT32_DTYPE);
  if ( xnidx==1 ) { xnidx=2; xidx[1]=0; }
  readMx(&ynidx, ydimSpecStr, yidx, INT32_DTYPE);
  if ( ynidx==1 ) { ynidx=2; yidx[1]=0; }

  err=compx2yIdx(xinfo,xnidx,xidx, yinfo,ynidx,yidx, &x2yIdx,&znd,&maccnd,&seqnd);
  
  if ( !err ) {
	 /* compute the mxInfo for the accumulated and rest sub-matrices */  
	 xmacc = mkemptymxInfo(maccnd); 
	 ymacc = mkemptymxInfo(maccnd);
	 /* N.B. xrest.sz holds the	size of the combined x and y rest matrix */
	 xrest = mkemptymxInfo(znd); 
	 yrest = mkemptymxInfo(znd);

	 err = initrestmaccmxInfo(znd, xinfo, yinfo, 
									  x2yIdx, xnidx, ynidx,
									  &xrest, &yrest, &xmacc, &ymacc);
  }

  if ( ! err ) {
	 /* compute the size of the output matrix */
	 zinfo=initzmxInfo(znd, xinfo, yinfo, x2yIdx, xnidx, ynidx); 
	 
	 /* optimise/standardize the query so its the way tprod wants it */
	 zrest = copymxInfo(zinfo);
	 err = optimisetprodQuery(&zrest, &xrest, &yrest, &xmacc, &ymacc);  
  }

  if ( ! err ) {
	 if ( xrest.rp != xinfo.rp ) { /* swap xinfo/yinfo if optimiseQuery did */
		MxInfo tmp = xinfo; xinfo=yinfo; yinfo=tmp;
	 }

	 /* init the memory for the result */
	 zinfo.rp = memalloc(zinfo.numel,zinfo.dtype); zrest.rp=zinfo.rp;
	 if(xinfo.ip!=0 || yinfo.ip!=0) zinfo.ip=memalloc(zinfo.numel,zinfo.dtype);
	 zrest.ip=zinfo.ip;
  }
  /* do the actual tprod call */
  /* call tprod to do the real work */
  /* do the (appropriately typed) operation */
  if ( ! err ) {
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
  }
  fprintf(stdout,"Zinfo =");printMxInfo(stdout,zrest);printf("\n");
     

  /*------------------------------------------------------------------------*/
  /* error handling code */
  /* check for errors */
  switch ( err ) {
  case ZTYPEMISMATCH : 
	 ERROR("tprod: Z is of unsupported type"); break;
  case XTYPEMISMATCH :
	 ERROR("tprod: X is of unsupported type"); break;
  case YTYPEMISMATCH :
	 ERROR("tprod: YTYPEMISMATCH is of unsupported type"); break;
  case INTYPEMISMATCH :
	 ERROR("tprod: input real/complex mix unsupported"); break;
  case OTHERERROR :
	 ERROR("tprod: Something else went wrong, sorry!"); break;
  case UNSUPPORTEDINPUTS :
	 ERROR("tprod: Inputs of unsupported type: only double/single"); break;
  default: ;
  }
  
  /* free up all the memory we've allocated */
  /* N.B. not clear we need to do this from the matlab web-site should happen
	  automatically */
  delmxInfo(&xinfo);  delmxInfo(&yinfo); delmxInfo(&zinfo);
  delmxInfo(&xmacc);  delmxInfo(&ymacc);
  delmxInfo(&xrest);  delmxInfo(&yrest); delmxInfo(&zrest);
  FREE(x2yIdx);
  FREE(xinfo.rp); FREE(xinfo.ip);
  FREE(yinfo.rp); FREE(yinfo.ip); 
  FREE(zinfo.rp); FREE(zinfo.ip);   
  return err;
}
