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
#include <strings.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
		int num, i;
		char *key = NULL, *str = NULL;
		mxArray *field = NULL, *defaultval = NULL;
		int index;

		if (nrhs<2 || nrhs>3)
				mexErrMsgTxt("incorrect number of input arguments");

		if (!mxIsChar(prhs[0]))
				mexErrMsgTxt("the key should be specified as a string");

		key = mxArrayToString(prhs[0]);
		num = mxGetNumberOfElements(prhs[1]);

		/* the default output will be dealt with later */
		plhs[0] = NULL;

		if (mxIsStruct(prhs[1])) {
				if (num!=1)
						mexErrMsgTxt("the second input should be a single structure");

				field = mxGetField(prhs[1], 0, key);
				if (field)
						plhs[0] = mxDuplicateArray(field);
		}
		else if (mxIsCell(prhs[1])) {
				if ((num % 2)!=0)
						mexErrMsgTxt("the second input should contain key-value pairs");

				for (i=0; i<num; i+=2) {
						if ((str = mxArrayToString(mxGetCell(prhs[1], i)))==NULL) {
								mexErrMsgTxt("the second input should contain key-value pairs");
						}
						else {
								if (strcasecmp(str, key)==0) {
										field = mxGetCell(prhs[1], i+1);
										if (field)
												plhs[0] = mxDuplicateArray(field);
										break;
								}
						}
				}
		}
		else {
				mexErrMsgTxt("the first input argument should be a cell-array or structure");
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

