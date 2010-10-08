#include <mex.h>
#include <stdio.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	/*
	int fid;
	const double *pn;

	if (nrhs!=1 || !mxIsDouble(prhs[0]) || mxGetM(prhs[0])*mxGetN(prhs[0])!=1) {
		mexErrMsgTxt("usage: flush(fid);");
	}

	pn  = mxGetPr(prhs[0]);
	fid = (int) pn[0];
	*/
	fflush(NULL);
}