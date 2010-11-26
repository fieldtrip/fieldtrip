/* File: mxUtils.C
	
	Header for set of utility functions to setup mxInfo arrays for testing 
*/
#include "mxUtils.h"

int parseInfo(int nd, char *typestr, char *szstr, int *type, int *complexp, 
					int *sz, int *numel){
  /* read a mxInfo format description string */
 char *cp=szstr;
 int i, nch;
 if ( typestr[0]=='s' ) *type=SINGLE_DTYPE; 
 else if ( typestr[0]=='d' ) *type=DOUBLE_DTYPE; 
 else if ( typestr[0]=='l' || typestr[0]=='b' ) *type=LOGICAL_DTYPE;
 else ERROR("Unsupported type string");
 if ( typestr[1]=='c' ) *complexp=1; else if ( typestr[1]=='r' ) *complexp=0; 
 else ERROR("Unsupported complex type"); 
 *numel=1;
 if ( *cp=='[' ) cp++;
 for (i=0; *cp!=0 && i < nd ; i++) { 
	nch = sscanf(cp,"%d",sz+i); 	
	if ( nch==0 || *cp==']' ) break;
	*numel *= sz[i];
	for ( ; (*cp>='0' && *cp<='9') || *cp=='.' || *cp=='-'; cp++ );
	if ( *cp=='x' ) cp++;
 } 
 if ( *cp==']' ) cp++;
 for ( nd--  ; i < nd ; nd--) { sz[nd]=1; } sz[nd]=1; /* count back to num dims used */
 return nd;
}

char* readMx(int *numel, char *datstr, void *data, int dtype){
  /* read a matrix of data values from a string */
 char *cp=datstr;
 int i, nch=0, tmp;
 if ( *cp=='[' ) cp++;
 for (i=0; *cp!=0 && i < *numel ; i++) { 
	switch (dtype){
	case INT32_DTYPE:  nch = sscanf(cp,"%d", ((int*)data)+i); 	    break;
	case DOUBLE_DTYPE: nch = sscanf(cp,"%lf",((double*)data)+i); 	 break;
	case SINGLE_DTYPE: nch = sscanf(cp,"%f", ((float*)data)+i); 	 break;
	case LOGICAL_DTYPE:nch = sscanf(cp,"%d", &tmp); ((char*)data)[i]=(char)tmp; break;
	default: ERROR("Unsupported data type");
	}
	if ( nch==0 || *cp==']' ) break;
	for ( ; (*cp>='0' && *cp<='9') || *cp=='.' || *cp=='-'; cp++ );
	if ( *cp==',' || *cp==' ' ) cp++;
 }
 if ( *cp==']' ) cp++;
 *numel=i;
 return cp; /* return the un-consumed characters */
}

FILE* freadMx(int *numel, FILE *fd, void *data, int dtype){
  /* read a matrix of data values from a file */
 char cp;
 int i, nch=0, tmp;
 cp=getc(fd); while ( cp=='\n' || cp==' ' || cp=='\t' ) cp=getc(fd); 
 if ( cp!='[' ) ungetc(cp,fd);
 for (i=0; i < *numel && !feof(fd) ; i++) { 
	switch (dtype){
	case INT32_DTYPE:  nch = fscanf(fd,"%d", ((int*)data)+i); 	    break;
	case DOUBLE_DTYPE: nch = fscanf(fd,"%lf",((double*)data)+i); 	 break;
	case SINGLE_DTYPE: nch = fscanf(fd,"%f", ((float*)data)+i); 	 break;
	case LOGICAL_DTYPE:nch = fscanf(fd,"%d", &tmp); ((char*)data)[i]=(char)tmp; break;
	default: ERROR("Unsupported data type");
	}
	if ( nch==0 ) break;
	cp = fgetc(fd);	
	if ( cp==',' || cp==' ') ; else ungetc(cp,fd);
 }
 cp=getc(fd); while ( cp=='\n' || cp==' ' || cp=='\t' ) cp=getc(fd); 
 if ( cp!=']' ) ungetc(cp,fd);
 *numel=i;
 return fd;
}

void randFill(int numel,double *dat, int dtype){
  int i;
  if ( dtype==DOUBLE_DTYPE ) {
	 for( i=0; i < numel; i++ ) dat[i]=((double)rand())/RAND_MAX;
  } else {
	 for( i=0; i < numel; i++ ) ((float*)dat)[i]=((float)rand())/RAND_MAX;
  }
}


int datFill(char *datSrc, int numel,void *rp, void *ip, int dtype){
  const char *RAND="rand"; 
  int i=0; for(; datSrc[i]!=0 && RAND[i]!=0 && datSrc[i]==RAND[i]; i++); /* strcmp */
  if ( datSrc[i]==0 ) { /*strcmp(datSrc,"rand")==0 ) {*/
	 randFill(numel,rp,dtype);
	 if ( ip ) randFill(numel,ip,dtype);	
  } else {
	 FILE *fd = fopen(datSrc,"r");
	 if ( !fd ) {
		datSrc = readMx(&numel,datSrc,rp,dtype);
		if ( ip ) datSrc=readMx(&numel,datSrc,ip,dtype);			
	 } else {
		fd=freadMx(&numel,fd,rp,dtype);
		if ( ip ) fd=freadMx(&numel,fd,ip,dtype);	
		fclose(fd);
	 }
  }
 return numel;
}


MxInfo loadMxInfo(FILE *fd){
  /* load mxInfo from a file */
  char typeStr[1024], szStr[1024]; /* string to read stuff into */
  MxInfo mx;
  int complex;
  int di;
  
  /* read and setup the meta-info */
  fscanf(fd,"%s",typeStr);
  fscanf(fd,"%s",szStr);
  mx=mkemptymxInfo(10);
  mx.nd=parseInfo(mx.nd, typeStr, szStr, &mx.dtype, &complex, mx.sz, &mx.numel);

  /* setup the strides */
  mx.stride[0]=1; for( di=1; di<mx.nd; di++ ) mx.stride[di]=mx.stride[di-1]*mx.sz[di-1]; mx.stride[mx.nd]=mx.numel;

  /* allocate ram to hold the data */
  mx.rp = CALLOC(mx.numel,dsz_bytes(mx));
  if ( complex==1 )  mx.ip = CALLOC(mx.numel,dsz_bytes(mx)); else mx.ip=0;

  /* read the data */
  freadMx(&mx.numel,fd,mx.rp,mx.dtype);
  if ( mx.ip ) 
	 freadMx(&mx.numel,fd,mx.ip,mx.dtype);

  return mx;
}
