#include <mex.h>
#include <MultiChannelFilter.h>


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	const double *A, *B, *X;
	double *Y;
	int downsample;
	int dimX, N, Nout;
	int order;
	int skipSamples = 0;

	if (nrhs!=4) mexErrMsgTxt("call myFilter(B,A,downsample,X)");

	order = mxGetNumberOfElements(prhs[0]);
	if (order != mxGetNumberOfElements(prhs[1])) mexErrMsgTxt("B and A must have same length");
	order--;

	B = mxGetPr(prhs[0]);
	A = mxGetPr(prhs[1]);

	downsample = (int) mxGetScalar(prhs[2]);
	if (downsample < 1) downsample = 1;

	dimX = mxGetM(prhs[3]);
	N = mxGetN(prhs[3]);
	X = mxGetPr(prhs[3]);

	MultiChannelFilter<double> *filter = new MultiChannelFilter<double>(dimX, order);
	filter->setCoefficients(B, A);

	Nout = (N - skipSamples + downsample - 1)/downsample;

	plhs[0] =  mxCreateDoubleMatrix(dimX, Nout, mxREAL);
	Y = mxGetPr(plhs[0]);

	for (int i=0;i<N;i++) {
		if (skipSamples == 0) {
			filter->process(Y, X);
			Y += dimX;
		} else {
			filter->process(X);
		}
		X+=dimX;
		if (--skipSamples < 0) skipSamples = downsample-1;
	}

	delete filter;
}
