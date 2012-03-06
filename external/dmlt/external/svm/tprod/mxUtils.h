#ifndef mxUtilsH
#define mxUtilsH
/* File: mxUtils.H
	
	Header for set of utility functions to setup mxInfo arrays for testing 
*/
#include "mxInfo.h"

void printInfo(FILE *fd,int xtype, int xreal, int* xsz);
int parseInfo(int nd, char *typestr, char *szstr, int *type, int *real, 
				  int *sz, int *numel);
char* readMx(int *numel, char *datstr, void *data, int dtype);
FILE* freadMx(int *numel, FILE *fd, void *data, int dtype);
int datFill(char *datSrc, int numel,void *rp, void *ip, int dtype);
void *memalloc(int numel, int dtype);

#endif
