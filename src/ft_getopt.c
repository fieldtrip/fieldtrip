/*
 * MEX file implementation of keyval: this function is called very frequently
 * and therefore worthwile to speed up with a compiled version.
 * 
 * Use as
 *   [val] = keyval(s, key, default)
 * where s is a structure or a cell-array. The default argument can be
 * empty, in which case [] is returned.
 * 
 * Copyright (C) 2011, Robert Oostenveld
 *
 */

#include "mex.h"
#include "matrix.h"
#include "compiler.h"

#if defined (COMPILER_LCC)
#include <string.h>
#define strcasecmp strcmpi
#elif defined (COMPILER_MSVC)
#include <string.h>
#define strcasecmp stricmp
#else
#include <strings.h>
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
		int num, i;
		char *key = NULL, *str = NULL;
		mxArray *field = NULL, *defaultval = NULL;
		int index;

		if (nrhs<2 || nrhs>3)
				mexErrMsgTxt("incorrect number of input arguments");

		if (!mxIsChar(prhs[1]))
				mexErrMsgTxt("the key should be specified as a string");

		key = mxArrayToString(prhs[1]);
		num = mxGetNumberOfElements(prhs[0]);

		/* the default output will be dealt with later */
		plhs[0] = NULL;

		if (mxIsStruct(prhs[0])) {
				if (num!=1)
						mexErrMsgTxt("the first input should be a single structure");

				field = mxGetField(prhs[0], 0, key);
				if (field)
						plhs[0] = mxDuplicateArray(field);
		}
		else if (mxIsCell(prhs[0])) {
				if ((num % 2)!=0)
						mexErrMsgTxt("the first input should contain key-value pairs");

				for (i=0; i<num; i+=2) {
						if ((str = mxArrayToString(mxGetCell(prhs[0], i)))==NULL) {
								mexErrMsgTxt("the first input should contain key-value pairs");
						}
						else {
								if (strcasecmp(str, key)==0) {
										field = mxGetCell(prhs[0], i+1);
										if (field)
												plhs[0] = mxDuplicateArray(field);
										break;
								}
						}
				}
		}
		else if (mxIsEmpty(prhs[0])) {
				/* do nothing, the default value will be assigned below */
		}
		else {
				mexErrMsgTxt("the first input argument should be a cell-array or structure");
		}

		if (plhs[0]!=NULL && mxIsEmpty(plhs[0])) {
				/* use the default value instead of the empty input that was specified:
				   this applies for example if you do functionname('key', []), where
				   the empty is meant to indicate that the user does not know or care
				   what the value is */
				mxDestroyArray(plhs[0]);
				plhs[0] = NULL;
		}

		if (plhs[0]==NULL) {
				/* the output value has not yet been assigned */
				if (nrhs==3)
						/* return the default value */
						plhs[0] = mxDuplicateArray(prhs[2]);
				else
						/* return an empty array */
						plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
		}
}

