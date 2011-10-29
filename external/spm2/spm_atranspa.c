#ifndef lint
static char sccsid[]="@(#)spm_atranspa.c	2.1 (c) John Ashburner 99/01/15";
#endif

#include "mex.h"

/* C = A'*A */
void atranspa(m,n,A,C)
int m,n;
double A[/* m,n */], C[/* n,n */];
{
	int i, j1,j2;
	double *p1, *p2, c;

	/* Generate half of symmetric matrix C */
	for (j1=0;j1<n;j1++)
	{
		p1 = &(A[j1*m]);
		for (j2=0;j2<=j1;j2++)
		{
			p2 = &(A[j2*m]);
			c = 0.0;

			/* Work down columns in inner loop
			   to reduce paging */
			for(i=0; i<m; i++)
				c += p1[i]*p2[i];

			C[j1*n+j2] = c;
		}
	}

	/* Generate other half */
	for(j1=0; j1<n; j1++)
		for (j2=0;j2<j1;j2++)
			C[j2*n+j1] = C[j1*n+j2];
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	unsigned int n;
	unsigned int m;
	double *C;
	double *A;

	if (nrhs == 0) mexErrMsgTxt("Incorrect usage.");
	if (nrhs != 1) mexErrMsgTxt("Only 1 input argument required.");
	if (nlhs > 1) mexErrMsgTxt("Only 1 output argument required.");

	if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]))
		mexErrMsgTxt("spm_atranspa: A must be numeric, real, full and double");
	A = mxGetPr(prhs[0]);
	m = mxGetM(prhs[0]);
	n = mxGetN(prhs[0]);

	plhs[0] = mxCreateDoubleMatrix(n,n, mxREAL);
	C = mxGetPr(plhs[0]);

	atranspa(m,n,A,C);

}

