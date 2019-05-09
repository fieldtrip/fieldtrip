/*
 * MEX file implementation of ft_getopt: this function is called very frequently
 * and therefore worthwile to speed up with a compiled version.
 * 
 * Use as
 *   val = ft_getopt(s, key, default)
 * where s is a structure or a cell-array.
 *
 * It will return the value of the option, or an empty array if the option was
 * not present.
 *
 * The optional fourth argument allows you to specify whether
 * or not an empty value in the configuration structure/cell-array should be
 * interpreted as meaningful. If emptymeaningful = 1, then an empty
 * configuration option will be returned if present. If emptymeaningful = 0,
 * then the specified default will be returned if an empty value is
 * encountered. The default value for emptymeaningful = 0.
 * 
 * Copyright (C) 2011, Robert Oostenveld
 *
 */

#include "mex.h"
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
		int num, i, emptymeaningful;
		char *key = NULL, *str = NULL;
		mxArray *field = NULL, *defaultval = NULL;
		int index;

		if (nrhs<2 || nrhs>4)
				mexErrMsgTxt("incorrect number of input arguments");

		if (!mxIsChar(prhs[1]))
				mexErrMsgTxt("the key should be specified as a string");
                
        if (nrhs == 4 && !(mxIsLogical(prhs[3]) || mxIsNumeric(prhs[3]))) {
            mexErrMsgTxt("if specified, input argument emptymeaningful should be a logical or numeric value");
        }
        
        if (nrhs < 4) {
            emptymeaningful = 0;
        } else {
            emptymeaningful = (int)mxGetScalar(prhs[3]);
        }

		key = mxArrayToString(prhs[1]);
		num = mxGetNumberOfElements(prhs[0]);

		/* the default output will be dealt with later */
		plhs[0] = NULL;

		if (mxIsClass(prhs[0], "config")) {
				/* the config object has to be converted to a struct object */
                /* this fixes bug 885 */
				mexPutVariable("caller", "bcks4i37yr3_cwb", prhs[0]);
				mexEvalString("bcks4i37yr3_cwb = struct(bcks4i37yr3_cwb);");
				prhs[0] = mexGetVariable("caller", "bcks4i37yr3_cwb");
				mexEvalString("clear bcks4i37yr3_cwb;");
		}

		if (mxIsStruct(prhs[0])) {
				/* it will also end up here if the input is an object, in which case this code fails */
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

		if (plhs[0]!=NULL && mxIsEmpty(plhs[0]) && !emptymeaningful) {
				/* use the default value instead of the empty input that was specified:
				   this applies for example if you do functionname('key', []), where
				   the empty is meant to indicate that the user does not know or care
				   what the value is */
				mxDestroyArray(plhs[0]);
				plhs[0] = NULL;
		}

		if (plhs[0]==NULL) {
				/* the output value has not yet been assigned */
				if (nrhs>=3)
						/* return the default value */
						plhs[0] = mxDuplicateArray(prhs[2]);
				else
						/* return an empty array */
						plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
		}
}

