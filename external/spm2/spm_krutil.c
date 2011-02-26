#ifndef lint
static char sccsid[]="@(#)spm_krutil.c	2.2 John Ashburner 03/01/21";
#endif

/********************************************************************************/
/* beta = kron(b2,b1)'*img(:)
 * m1  - rows in b1
 * m2  - rows in b2
 * n1  - columns in b1
 * n2  - columns in b2
 * img - m2*m2 vector
 * b1  - basis functions - x
 * b2  - basis functions - y
 * beta - resulting vector
*/
void kronutil1(int n1, int n2, int m1, int m2,
	double img[], double b1[], double b2[], double beta[])
{
	int j1,j2, i1,i2;
	double beta1[32];

	/* Zero beta */
	for(j1=0; j1<n1*n2; j1++) beta[j1]=0.0;

	for(i2=0; i2<m2; i2++)
	{
		/* generate small beta */
		for(j1=0; j1<n1; j1++) beta1[j1]=0.0;
		for(i1=0; i1<m1; i1++)
		{
			double wt = img[i1+m1*i2];
			for(j1=0; j1<n1; j1++)
				beta1[j1] += b1[i1+j1*m1]*wt;
		}

		/* kronecker tensor product to increment beta */
		for(j2=0; j2<n2; j2++)
		{
			double wt = b2[i2+j2*m2];
			double *ptrb = beta+(n1*j2);
			for(j1=0; j1<n1; j1++)
				ptrb[j1] += wt*beta1[j1];
		}
	}
}

/********************************************************************************/
/* alpha = kron(b2,b1)'*diag(img(:))*kron(b2,b1)
 * m1  - rows in b1
 * m2  - rows in b2
 * n1  - columns in b1
 * n2  - columns in b2
 * img - m1*m2 vector
 * b1  - basis functions - x
 * b2  - basis functions - y
 * alpha - resulting matrix
*/
void kronutil2(int n1, int n2, int m1, int m2,
	double img[], double b1[], double b2[], double alpha[])
{
	int j11,j12, j21,j22, i1, i2;
	double alpha1[1024];

	/* Zero alpha */
	for(j21=0; j21<n1*n2; j21++)
		for(j11=0; j11<=j21; j11++)
			alpha[j11+j21*n1*n2]=0.0;

	for(i2=0; i2<m2; i2++)
	{
		/* zero small alpha */
		for(j21=0; j21<n1; j21++)
			for(j11=0; j11<=j21; j11++)
				alpha1[j11+j21*n1]=0.0;

		/* generate upper half of small alpha */
		for(i1=0; i1<m1; i1++)
		{
			double wt = img[i1+m1*i2];
			for(j21=0; j21<n1; j21++)
			{
				double wt2 = wt*b1[i1+j21*m1];
				for(j11=0; j11<=j21; j11++)
					alpha1[j11+j21*n1] += wt2*b1[i1+j11*m1];
			}
		}

		/* kronecker tensor product to increment upper half of large alpha */
		for(j22=0; j22<n2; j22++)
		{
			double wt = b2[i2+j22*m2];
			for(j12=0; j12<=j22; j12++)
			{
				double *ptra = alpha+(n1*j12+n1*n1*n2*j22);
				double wt2 = wt*b2[i2+j12*m2];

				for(j21=0; j21<n1; j21++)
					for(j11=0; j11<=j21; j11++)
						ptra[j11+j21*n1*n2] += wt2*alpha1[j11+j21*n1];
			}
		}
	}
	/* fill in lower triangle of symmetric submatrices of alpha */
	for(j22=0; j22<n2; j22++)
		for(j12=0; j12<j22; j12++)
		{
			double *ptra = alpha+(n1*j12+n1*n1*n2*j22);
			for(j21=0; j21<n1; j21++)
				for(j11=0; j11<j21; j11++)
					ptra[j21+j11*n1*n2] = ptra[j11+j21*n1*n2];
		}

	/* fill in lower triangle of alpha */
	for(j22=0; j22<n2*n1; j22++)
		for(j12=0; j12<j22; j12++)
			alpha[j22+j12*n1*n2]=alpha[j12+j22*n1*n2];
}

/********************************************************************************/
/* alpha = kron(B1y,B1x)'*diag(img(:))*kron(B2y,B2x)
*/
void kronutil3(int n1x, int n1y, int n2x, int n2y, int m1, int m2,
	double img[], double b1x[], double b1y[], double b2x[], double b2y[], double alpha[])
{
	int j1x,j1y, j2x,j2y, i1, i2;
	double alpha1[1024];

	/* Zero alpha */
	for(j2x=0; j2x<n2y*n2x; j2x++)
		for(j1x=0; j1x<n1y*n1x; j1x++)
			alpha[j1x+j2x*n1y*n1x]=0.0;

	for(i2=0; i2<m2; i2++)
	{
		/* zero small alpha */
		for(j2x=0; j2x<n2x; j2x++)
			for(j1x=0; j1x<n1x; j1x++)
				alpha1[j1x+j2x*n1x]=0.0;

		/* generate upper half of small alpha */
		for(i1=0; i1<m1; i1++)
		{
			double wt = img[i1+m1*i2];
			for(j2x=0; j2x<n2x; j2x++)
			{
				double wt2 = wt*b2x[i1+j2x*m1];
				for(j1x=0; j1x<n1x; j1x++)
					alpha1[j1x+j2x*n1x] += wt2*b1x[i1+j1x*m1];
			}
		}

		/* kronecker tensor product to increment upper half of large alpha */
		for(j2y=0; j2y<n2y; j2y++)
		{
			double wt2 = b2y[i2+j2y*m2];
			for(j1y=0; j1y<n1y; j1y++)
			{
				double *ptra = alpha+(n1x*j1y + n1y*n1x*n2x*j2y);
				double wt1   = wt2*b1y[i2+j1y*m2];

				for(j2x=0; j2x<n2x; j2x++)
					for(j1x=0; j1x<n1x; j1x++)
						ptra[j1x+j2x*n1y*n1x] += wt1*alpha1[j1x+j2x*n1x];
			}
		}
	}
}

/********************************************************************************/
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	unsigned int n1x, n1y, n2x, n2y, m1, m2;

	double *alpha, *beta, *img, *b1x, *b1y, *b2x, *b2y;

	if (nrhs == 0) mexErrMsgTxt("Incorrect usage");

	if (nrhs != 4 & nrhs != 5) mexErrMsgTxt("4 or 5 input arguments required");
	if (nlhs > 1) mexErrMsgTxt("only 1 output argument required");

	if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]))
		mexErrMsgTxt("img must be numeric, real, full and double");
	img = mxGetPr(prhs[0]);
	m1 = mxGetM(prhs[0]);
	m2 = mxGetN(prhs[0]);

	if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1]))
		mexErrMsgTxt("b1x must be numeric, real, full and double");
	b1x = mxGetPr(prhs[1]);
	if (mxGetM(prhs[1]) != m1)
		mexErrMsgTxt("b1x has incompatible m dimension");
	n1x = mxGetN(prhs[1]);

	if (!mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2]) || mxIsSparse(prhs[2]) || !mxIsDouble(prhs[2]))
		mexErrMsgTxt("b1y must be numeric, real, full and double");
	b1y = mxGetPr(prhs[2]);
	if (mxGetM(prhs[2]) != m2)
		mexErrMsgTxt("b1y has incompatible m dimension");
	n1y = mxGetN(prhs[2]);

	if (nrhs == 4)
	{
		if (!mxIsNumeric(prhs[3]) || mxIsComplex(prhs[3]) || mxIsSparse(prhs[3]) || !mxIsDouble(prhs[3]))
			mexErrMsgTxt("flag must be numeric, real, full and double");
		if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 1)
			mexErrMsgTxt("flag must contain only one element");

		if (!*mxGetPr(prhs[3]))
		{
			plhs[0] = mxCreateDoubleMatrix(n1x*n1y,1, mxREAL);
			beta = mxGetPr(plhs[0]);
			kronutil1(n1x,n1y,m1,m2,img,b1x,b1y, beta);
		}
		else
		{
			plhs[0] = mxCreateDoubleMatrix(n1x*n1y,n1x*n1y, mxREAL);
			alpha = mxGetPr(plhs[0]);
			kronutil2(n1x,n1y,m1,m2,img,b1x,b1y, alpha);
		}
	}
	else
	{
		if (!mxIsNumeric(prhs[3]) || mxIsComplex(prhs[3]) || mxIsSparse(prhs[3]) || !mxIsDouble(prhs[3]))
			mexErrMsgTxt("b2x must be numeric, real, full and double");
		b2x = mxGetPr(prhs[3]);
		if (mxGetM(prhs[3]) != m1)
			mexErrMsgTxt("b2x has incompatible m dimension");
		n2x = mxGetN(prhs[3]);

		if (!mxIsNumeric(prhs[4]) || mxIsComplex(prhs[4]) || mxIsSparse(prhs[4]) || !mxIsDouble(prhs[4]))
			mexErrMsgTxt("b1y must be numeric, real, full and double");
		b2y = mxGetPr(prhs[4]);
		if (mxGetM(prhs[4]) != m2)
			mexErrMsgTxt("b2y has incompatible m dimension");
		n2y = mxGetN(prhs[4]);

		plhs[0] = mxCreateDoubleMatrix(n1x*n1y,n2x*n2y, mxREAL);
		alpha = mxGetPr(plhs[0]);
		kronutil3(n1x,n1y,n2x,n2y,m1,m2,img,b1x,b1y,b2x,b2y, alpha);
	}
}
