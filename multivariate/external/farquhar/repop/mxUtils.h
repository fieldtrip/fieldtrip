#ifndef mxUtilsH
#define mxUtilsH
/* File: mxUtils.H
	
	Header for set of utility functions to setup mxInfo arrays for testing 
*/
#include "mxInfo.h"

int parseInfo(int nd, char *typestr, char *szstr, int *type, int *complexp, 
				  int *sz, int *numel);
char* readMx(int *numel, char *datstr, void *data, int dtype);
FILE* freadMx(int *numel, FILE *fd, void *data, int dtype);
int datFill(char *datSrc, int numel,void *rp, void *ip, int dtype);
MxInfo loadMxInfo(FILE *fd);

#endif
