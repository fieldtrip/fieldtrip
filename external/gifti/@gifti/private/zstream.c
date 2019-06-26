/*
 * $Id: zstream.c 6417 2015-04-21 16:03:44Z guillaume $
 * Guillaume Flandin
 */

/* mex -O CFLAGS='$CFLAGS -std=c99' -largeArrayDims zstream.c */

/* setenv CFLAGS "`mkoctfile -p CFLAGS` -std=c99" */
/* mkoctfile --mex zstream.c */

/* miniz: http://code.google.com/p/miniz/ */
#define MINIZ_NO_STDIO
#define MINIZ_NO_ARCHIVE_APIS
#define MINIZ_NO_TIME
#define MINIZ_NO_ZLIB_APIS
#include "miniz.c"

#include "mex.h"

/* --- GATEWAY FUNCTION --- */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

char *action;
unsigned char *IN = NULL, *OUT = NULL;
size_t INlen, OUTlen;

/* Check for proper number of arguments */
if (nrhs < 2)
    mexErrMsgTxt("Not enough input arguments.");
else if (nrhs > 2)
    mexErrMsgTxt("Too many input arguments.");
else if (nlhs > 1)
    mexErrMsgTxt("Too many output arguments.");

/* The input ACTION must be a string */
if (!mxIsChar(prhs[0]))
    mexErrMsgTxt("Input ACTION must be a string.");
action = mxArrayToString(prhs[0]);

/* The input IN must be a real uint8 array */
if (!mxIsUint8(prhs[1]) || mxIsComplex(prhs[1]))
    mexErrMsgTxt("Input IN must be a real uint8 array.");

INlen = mxGetNumberOfElements(prhs[1]);
IN = mxGetData(prhs[1]);

if (!strcmp(action,"D")) {

    /* Decompress data */
    OUT = tinfl_decompress_mem_to_heap(IN, INlen, &OUTlen, TINFL_FLAG_PARSE_ZLIB_HEADER);

    if (OUT == NULL)
        mexErrMsgTxt("Error when decompressing data.");
}
else if (!strcmp(action,"C")) {
    /* Compress data */
    OUT = tdefl_compress_mem_to_heap(IN, INlen, &OUTlen, TDEFL_WRITE_ZLIB_HEADER);
    
    if (OUT == NULL)
        mexErrMsgTxt("Error when compressing data.");
}
else {
    mexErrMsgTxt("Unknown ACTION type.");
}

/*  */
plhs[0] = mxCreateNumericMatrix(OUTlen,1,mxUINT8_CLASS,mxREAL);
if (plhs[0] == NULL)
    mexErrMsgTxt("Error when creating output variable.");

memcpy(mxGetData(plhs[0]), OUT, OUTlen);

mxFree(action);
mz_free(OUT);

}
