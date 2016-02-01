/* 
	File:
	
	Function to do repop calls without the matlab dependencies.  Used for
	memory debugging or as an example for integration of repop with other
	numeric libraries.

	Calling Convention:
	repop_testcases xtypeSpec xsize xdataSrc opNm ytypeSpec ysize ydataSrc = ztypeSpec zsize zdataSrc OPTIONS

	Options:
	  i - inplace
     n - non-integer replication
     m - integer Multiple replication
	  [0-9]* - number of repitions to perform

	Examples
	repop_testcases dr 2x1 "1 2" "+" dr 2x1 "1 2" = dr 2x1 "2 4" 
	OR
	repop_testcases dr 2x1 "1 2" "+" dr 2x1 "1 2" 1000i

*/

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "mxInfo.h"
#include "repop.h"
#include "mxUtils.h"


/* first provide necessary call-back functions */
void *CALLOC(size_t nmemb, size_t size){  return (void*)calloc(nmemb,size); }
void *MALLOC(size_t size) { return (void*)malloc(size); }
void FREE(void *ptr) {  free(ptr); }
void ERROR(const char *msg) { fprintf(stderr,"%s\n",msg); }
void WARNING(const char *msg) { fprintf(stderr,"%s\n",msg); }

float EPS=1e-6;

void cmdLineProb(int argc,char **argv, char *opnmStr, MxInfo *xinfo, MxInfo *yinfo, MxInfo *zinfo, int *nReps, int *repNonUnitDim);
void fileProb(FILE *fd, char *opnmStr, MxInfo *xinfo, MxInfo *yinfo, MxInfo *zinfo);
MxInfo readMxInfoStr(char *typestr,char *szstr, char *datSrc);

double getVal(MxInfoDTypes dtype, void *p, int idx){
  double val=.0;
  switch ( dtype ) {
  case LOGICAL_DTYPE: val = ((char*)p)[idx]; break;
  case CHAR_DTYPE   : val = ((char*)p)[idx]; break;
  case DOUBLE_DTYPE : val = ((double*)p)[idx]; break;
  case SINGLE_DTYPE : val = ((float*)p)[idx]; break;
  }
  return val;
}


int main(int argc, char **argv){
  
  MxInfo xinfo, yinfo, truezinfo;
  char opnmStr[1024];
  int res=0, nReps=1, repNonUnitDim=1;
  clock_t eTime = 0l;
  
  eTime=clock();
  if ( argc > 2 ) { /* command line test problem specification */
	 cmdLineProb(argc,argv,opnmStr,&xinfo,&yinfo,&truezinfo,&nReps,&repNonUnitDim);
	 res = unitTest(opnmStr,xinfo,yinfo,truezinfo,nReps,repNonUnitDim);
  
	 /* free the RAM we've allocated */
	 delmxInfo(&xinfo); FREE(xinfo.rp); FREE(xinfo.ip);
	 delmxInfo(&yinfo); FREE(yinfo.rp); FREE(yinfo.ip); 
	 delmxInfo(&truezinfo); FREE(truezinfo.rp); FREE(truezinfo.ip); 

  } else { /* file contains the set of unit tests to run */
	 /* Get the different parts of the input arguments */
	 int subProb=1, spres=0;
	 FILE *fd;
	 char cp;
	 fd=fopen(argv[1],"r"); 
	 if ( fd == 0 ) { ERROR("Couldn't open the file"); exit(1); }

	 while ( !feof(fd) ) {
		fprintf(stdout,"Unit-Test #%d\n",subProb++);
		
		/* load the data */
		fileProb(fd,opnmStr,&xinfo,&yinfo,&truezinfo);
		/* run the test */
		spres = unitTest(opnmStr,xinfo,yinfo,truezinfo,nReps,repNonUnitDim);
		res = res + spres;
		
		/* free the RAM we've allocated */
		delmxInfo(&xinfo); FREE(xinfo.rp); FREE(xinfo.ip);
		delmxInfo(&yinfo); FREE(yinfo.rp); FREE(yinfo.ip); 
		delmxInfo(&truezinfo); FREE(truezinfo.rp); FREE(truezinfo.ip); 		
		fprintf(stdout,"---------\n\n");
		cp=getc(fd); while ( cp=='\n' || cp==' ' || cp=='\t' ) cp=getc(fd); 
		if ( !feof(fd) ) ungetc(cp,fd); /* consume white-space */
	 }
	 fclose(fd);
  }

  eTime=clock()-eTime;
  fprintf(stdout,"\n\nnReps:    \t%d\nTotalTime:\t%f\nAveTime:  \t%f\n", nReps,((double)eTime)/(double)CLOCKS_PER_SEC,((double)eTime)/((double)CLOCKS_PER_SEC*nReps));
 
  exit(res);
}


int unitTest(char *opnmStr, MxInfo xinfo, MxInfo yinfo, MxInfo truezinfo, int nReps, int repNonUnitDim){
  MxInfo zinfo;
  /* get the operator id */
  int opid = getOpid(opnmStr);  
  int retVal=0, repi=0; 
  if ( opid < 0 ) 
	 ERROR("repop: Unrecognised operator. Must be one of: + - * / \\ = < > ~= <= >=");

  fprintf(stdout,"Xinfo =");printMxInfoSummary(stdout,xinfo);printf("\n");
  fprintf(stdout,"%s\n",opnmStr);
  fprintf(stdout,"Yinfo =");printMxInfoSummary(stdout,yinfo);printf("\n");
  if ( truezinfo.numel ) {
	 fprintf(stdout,"=\n");  
	 fprintf(stdout,"trueZinfo =");printMxInfoSummary(stdout,truezinfo);printf("\n");
  }

  /* init the result info */
  zinfo=initzinfo(xinfo,yinfo,repNonUnitDim); /* get its size */
  if ( opid >= EQ && opid <= GE ) zinfo.dtype = LOGICAL_DTYPE;
  /* init the memory for the result */
  zinfo.rp = CALLOC(zinfo.numel,dsz_bytes(zinfo));
  if ( ( opid == POWER || xinfo.ip != 0 || yinfo.ip != 0 ) && ( zinfo.dtype == SINGLE_DTYPE || zinfo.dtype == DOUBLE_DTYPE ) ) {
	 zinfo.ip = CALLOC(zinfo.numel,dsz_bytes(zinfo));
  } else { 
	 zinfo.ip = 0;
  }
  
  /* do the actual repop call */
  for ( repi=0; repi<nReps; repi++){
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
		ERROR("repop: Inputs of unsupported type: only double/single");		
	 }
	 if ( opid == POWER ) removeZeroImag(&zinfo);
  }
  
  /* check for errors */
  switch ( retVal ) {
  case ZTYPEMISMATCH : 	 ERROR("repop: Z is of unsupported type"); break;
  case XTYPEMISMATCH :	 ERROR("repop: X is of unsupported type"); break;
  case YTYPEMISMATCH :	 ERROR("repop: Y is of unsupported type"); break;
  case INTYPEMISMATCH :  ERROR("repop: input real/complex mix unsupported"); break;
  case UNDEFOPERATOR :   ERROR("repop: Unsupported operator requested"); break;
  case OTHERERROR :      ERROR("repop: Something else went wrong, sorry!"); break;
  default: ;
  }  

  fprintf(stdout,"Zinfo =");printMxInfoSummary(stdout,zinfo);printf("\n");
  
  /* check if the result is numerically correct */
  if ( truezinfo.numel ) {
	 float numAcc=0, tmp=0; 
	 int i, numOK=1; 
	 if ( truezinfo.dtype != zinfo.dtype ) { 
		ERROR("repop: different types"); numOK=0; retVal=1;
	 }
	 if ( truezinfo.numel != zinfo.numel || truezinfo.nd != zinfo.nd ) { 
		ERROR("repop: different sizes"); numOK=0; retVal=1;
	 }
	 for( i=0; i<zinfo.nd; i++ ){ 
		if ( zinfo.sz[i]!=truezinfo.sz[i] ) { 
		  ERROR("repop: different sizes"); numOK=0; retVal=1; break; 
		}		
	 }
	 for( i=0; i<zinfo.nd; i++ ){ 
		if ( zinfo.stride[i]!=truezinfo.stride[i]) { 
		  ERROR("repop: different strides"); numOK=0; retVal=1; break; 
		} 
	 }
	 for( i=0; i<zinfo.numel; i++ ) { /* test numerical correctness */
		tmp = fabs(getVal(zinfo.dtype,zinfo.rp,i)-getVal(truezinfo.dtype,truezinfo.rp,i));
		numAcc = (numAcc>tmp)?numAcc:tmp; /* max */
		if (zinfo.ip) {
		  tmp = fabs(getVal(zinfo.dtype,zinfo.ip,i)-getVal(truezinfo.dtype,truezinfo.ip,i));
		  numAcc = (numAcc>tmp)?numAcc:tmp; /* max */
		}
		if ( numAcc > EPS ) { 
		  ERROR("repop: numerical accuracy"); fprintf(stdout," = %f",numAcc); numOK=0; retVal=1; break; 
		}
	 }
	 fprintf(stdout,"\nNumeric Test = "); if ( numOK ) fprintf(stdout,"Passed\n"); else fprintf(stdout,"FAILED!!!\n");
	 retVal = retVal & !numOK;
  }
  
  /* free up memory we've allocated */
  delmxInfo(&zinfo); FREE(zinfo.rp); FREE(zinfo.ip); 
  return retVal;
}


void cmdLineProb(int argc,char **argv, char *opnmStr, MxInfo *xinfo, MxInfo *yinfo, MxInfo *zinfo, int *nReps, int *repNonUnitDim){
  /* Get the different parts of the input arguments */
  char *xtypestr=0, *xszstr=0, *xdatSrc=0;
  char *ytypestr=0, *yszstr=0, *ydatSrc=0;
  char *ztypestr=0, *zszstr=0, *zdatSrc=0;
  char *optStr;
  int i, nch;
  int inPlace=0;
  if ( argc < 7 ) {	 ERROR("Need at least 7 inputs");	 exit(1);   }

  xtypestr= argv[1]; xszstr=argv[2]; xdatSrc=argv[3];
  for(i=0; argv[4][i]!=0; i++) opnmStr[i] = argv[4][i]; opnmStr[i]=0;  /* copy the name */
  ytypestr= argv[5]; yszstr=argv[6]; ydatSrc=argv[7]; 
  if ( argc>=11 && argv[8][0]=='=' ) {
	 ztypestr=argv[9]; zszstr=argv[10]; zdatSrc=argv[11];
	 optStr  = (argc>11)?argv[12]:0;
  } else {
	 optStr  = (argc>7)?argv[8]:0;
  }
      
  /* extract the requested size info from the inputs */
  *xinfo = readMxInfoStr(xtypestr,xszstr,xdatSrc);
  *yinfo = readMxInfoStr(ytypestr,yszstr,ydatSrc);
  if ( ztypestr ) { *zinfo = readMxInfoStr(ztypestr,zszstr,zdatSrc); 
  } else {          *zinfo = mkemptymxInfo(0); 
  }

  /*   fprintf(stdout,"X =");printInfo(stdout,xdtype,xreal,xsz);printf("\n"); */
  /*   fprintf(stdout,"Y =");printInfo(stdout,ydtype,yreal,ysz);printf("\n"); */

  /* get the options */
  if ( optStr ) {	 
	 for (i=0 ; optStr[i] != 0; i++ ) {
		switch ( optStr[i] ) {
		case 'i': case 'I':   inPlace=1;  break;
		case 'm': case 'M':   *repNonUnitDim=1;  break;
		case 'n': case 'N':   *repNonUnitDim=2;  break;
		case '0': case '1': case '2':	case '3': case '4': case '5': case '6': case '7': case '8': case '9': /* nReps */
		  nch=sscanf(optStr+i,"%d",nReps);
		  while ( optStr[i]>='0' && optStr[i]<='9' ) i++; i--; /* skip the numbers */
		  break;
		default: 
		  WARNING("repop: Unrecognised option");
		}
	 } 
  }
}


void fileProb(FILE *fd, char *opnmStr, MxInfo *xinfo, MxInfo *yinfo, MxInfo *zinfo){
  char tmpStr[1024];
  /* read X */
  *xinfo=loadMxInfo(fd);
  if ( feof(fd) ) return;
  
  /* read Operator */
  fscanf(fd,"%s",opnmStr);
  if ( feof(fd) ) return;
  
  /* read Y */
  *yinfo=loadMxInfo(fd);
  if ( feof(fd) ) return;

  /* read = */
  fscanf(fd,"%s",tmpStr);
  if ( feof(fd) ) return;
  if ( tmpStr[0]!='=' ) {
	 ungetc(tmpStr[0],fd);
	 zinfo->sz=0; /* set as empty zinfo */
  } else {
	 /* read z */
	 *zinfo=loadMxInfo(fd);
  }
}


MxInfo readMxInfoStr(char *typestr,char *szstr, char *datSrc){
  MxInfoDTypes dtype;
  int complexp;
  int nd=10, numel;
  int sz[10+1];
  void *rp, *ip;

  nd=parseInfo(nd,typestr,szstr,&dtype,&complexp,sz,&numel);
  rp = CALLOC(numel,sizeof_dtype(dtype));
  if ( complexp ) ip = CALLOC(numel,sizeof_dtype(dtype)); else ip=0;
  datFill(datSrc,numel,rp,ip,dtype);   /* Initialise the values */
  return mkmxInfo(nd,sz,rp,ip,dtype); /* init a mxInfo equivalent */
}
