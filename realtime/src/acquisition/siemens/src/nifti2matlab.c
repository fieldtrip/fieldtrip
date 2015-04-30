/*
 * Copyright (C) 2010, Stefan Klanke
 * Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
 * Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 */

#include <mex.h>
#include <matrix.h>
#include <nifti1.h>

/* incomplete, but more to follow */
const char *nifti1_field_names[] = {
    "dim",
    "pixdim",
	"datatype"
};

static int nifti1_num_fields = 3;

mxArray *createFromNifti1(nifti_1_header *NH) {
	mxArray *S, *A;
	double *a;
	int i,dims;
	
	S = mxCreateStructMatrix(1, 1, nifti1_num_fields, nifti1_field_names);
	
	dims = NH->dim[0];
	
	/* dim */
	A = mxCreateDoubleMatrix(dims, 1, mxREAL);
	a = mxGetPr(A);
	for (i=0;i<dims;i++) {
		a[i] = NH->dim[i+1];
	}
	mxSetFieldByNumber(S, 0, 0, A);
	
	/* pixdim */
	A = mxCreateDoubleMatrix(dims, 1, mxREAL);
	a = mxGetPr(A);
	for (i=0;i<dims;i++) {
		a[i] = NH->pixdim[i+1];
	}
	mxSetFieldByNumber(S, 0, 1, A);
	
	/* datatype */
	A = mxCreateDoubleScalar(NH->datatype);
	mxSetFieldByNumber(S, 0, 2, A);
	
	return S;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	nifti_1_header *NH;
	int size;
	
	if (nrhs!=1 || !mxIsUint8(prhs[0])) mexErrMsgTxt("This function needs exactly one (uint8) argument.");
	
	size = mxGetNumberOfElements(prhs[0]);
	
	if (size < sizeof(nifti_1_header)) {
		plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
		return;
	}
	
	NH = (nifti_1_header *) mxGetData(prhs[0]);
	if (strcmp(NH->magic, "ni1")==0 || strcmp(NH->magic,"n+1")==0) {
		plhs[0] = createFromNifti1(NH);
		return;
	} 
	plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
}