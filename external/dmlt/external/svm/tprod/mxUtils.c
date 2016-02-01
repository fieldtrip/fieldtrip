/* File: mxUtils.C
	
	Header for set of utility functions to setup mxInfo arrays for testing 
*/
#include "mxUtils.h"

void printInfo(FILE *fd,int xtype, int xreal, int* xsz){
  int i;
  fprintf(fd,"[");
  for ( i=0; i < 5; i++ ) fprintf(fd,"%d x ",xsz[i]);
  fprintf(fd," ]");
  if ( xtype==SINGLE_DTYPE ) fprintf(fd," (single"); else fprintf(fd," (double");
  if ( xreal ) fprintf(fd,")"); else fprintf(fd," complex)");
}

int parseInfo(int nd, char *typestr, char *szstr, int *type, int *real, 
					int *sz, int *numel){
 char *cp=szstr;
 int i, nch;
 if ( typestr[0]=='s' ) *type=SINGLE_DTYPE; else *type=DOUBLE_DTYPE;
 if ( typestr[1]=='c' ) *real=0; else *real=1;
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
 nd=i;
 for (   ; i < 5 ; i++) { sz[i]=1; }
 return nd;
}

char* readMx(int *numel, char *datstr, void *data, int dtype){
  /* read a matrix of data values from a string */
 char *cp=datstr;
 int i, nch;
 if ( *cp=='[' ) cp++;
 for (i=0; *cp!=0 && i < *numel ; i++) { 
	switch (dtype){
	case INT32_DTYPE:  nch = sscanf(cp,"%d", ((int*)data)+i); 	    break;
	case DOUBLE_DTYPE: nch = sscanf(cp,"%lf",((double*)data)+i); 	 break;
	case SINGLE_DTYPE: nch = sscanf(cp,"%f", ((float*)data)+i); 	 break;
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
 int i, nch;
 cp=getc(fd); if ( cp!='[' ) ungetc(cp,fd);
 for (i=0; i < *numel ; i++) { 
	switch (dtype){
	case INT32_DTYPE:  nch = fscanf(fd,"%d", ((int*)data)+i); 	    break;
	case DOUBLE_DTYPE: nch = fscanf(fd,"%lf",((double*)data)+i); 	 break;
	case SINGLE_DTYPE: nch = fscanf(fd,"%f", ((float*)data)+i); 	 break;
	default: ERROR("Unsupported data type");
	}
	if ( nch==0 ) break;
	cp = fgetc(fd);	
	if ( cp==',' || cp==' ') ; else ungetc(cp,fd);
 }
 cp=getc(fd); if ( cp!=']' ) ungetc(cp,fd);
 *numel=i;
 return fd;
}

void *memalloc(int numel, int dtype){
  int datasize;
  switch(dtype){
  case DOUBLE_DTYPE : datasize=sizeof(double); break;
  case SINGLE_DTYPE : datasize=sizeof(float);  break;
  case INT32_DTYPE : datasize=sizeof(int);     break;
  case LOGICAL_DTYPE : datasize=sizeof(char);  break;
  otherwise: ERROR("Unsupported datatype");
  }
  return CALLOC(datasize,numel);
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
if ( strcmp(datSrc,"rand")==0 ) {
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
