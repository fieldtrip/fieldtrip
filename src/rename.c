/*
 * Copyright (C) 2011, Robert Oostenveld
 *
 * This implements an atomic rename by calling the UNIX/POSIX operating system command.
 *
 */

#include <stdio.h>

#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
		char *old, *new;
		int status;

		if (nrhs!=2)
				mexErrMsgTxt("incorrect number of input arguments");

		if (!mxIsChar(prhs[0]))
				mexErrMsgTxt("invalid input argument #1");

		if (!mxIsChar(prhs[1]))
				mexErrMsgTxt("invalid input argument #2");

		old = mxArrayToString(prhs[0]);
		new = mxArrayToString(prhs[1]);
		status = rename(old, new);

		plhs[0] = mxCreateDoubleScalar(status);

		return;
}
