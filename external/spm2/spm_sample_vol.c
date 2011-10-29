#ifndef lint
static char sccsid[]="@(#)spm_sample_vol.c	2.1 (c) John Ashburner 99/01/15";
#endif

#include "spm_mapping.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
/*
void mexFunction(nlhs, plhs, nrhs, prhs)
int nlhs, nrhs;
mxArray *plhs[];
const mxArray *prhs[];
*/
{
	MAPTYPE *map, *get_maps();
	int m,n, k, hold;
	double background=0.0;

	if (nrhs != 5 || nlhs > 4)
		mexErrMsgTxt("Inappropriate usage.");

	map=get_maps(prhs[0], &n);
	if (n!=1)
	{
		free_maps(map, n);
		mexErrMsgTxt("Bad image handle dimensions.");
	}

	for(k=1; k<=3; k++)
		if (!mxIsNumeric(prhs[k]) || mxIsComplex(prhs[k]) ||
			mxIsSparse(prhs[k]) || !mxIsDouble(prhs[k]))
		{
			free_maps(map, 1);
			mexErrMsgTxt("Coordinates must be numeric, real, full and double.");
		}

	m = mxGetM(prhs[1]);
	n = mxGetN(prhs[1]);
	if (mxGetM(prhs[2]) != m || mxGetN(prhs[2]) != n ||
		mxGetM(prhs[3]) != m || mxGetN(prhs[3]) != n)
		mexErrMsgTxt("Coordinates must have compatible dimensions.");

	if (!mxIsNumeric(prhs[4]) || mxIsComplex(prhs[4]) ||
		mxIsSparse(prhs[4]) || !mxIsDouble(prhs[4]) ||
		(mxGetM(prhs[4])*mxGetN(prhs[4]) != 1 && mxGetM(prhs[4])*mxGetN(prhs[4]) != 2))
	{
		free_maps(map, 1);
		mexErrMsgTxt("Bad hold & background argument.");
	}

	hold = (int)(*(mxGetPr(prhs[4])));

	if (abs(hold) > 127)
	{
		free_maps(map, 1);
		mexErrMsgTxt("Bad hold value.");
	}

	if (mxGetM(prhs[4])*mxGetN(prhs[4]) > 1)
		background = mxGetPr(prhs[4])[1];

	if (nlhs<=1)
	{
		plhs[0] = mxCreateDoubleMatrix(m,n,mxREAL);

		resample(m*n, map, mxGetPr(plhs[0]),
			mxGetPr(prhs[1]),mxGetPr(prhs[2]),mxGetPr(prhs[3]),
			hold, background);
	}
	else
	{
		if (hold==0)
		{
			free_maps(map, 1);
			mexErrMsgTxt("This wont work for nearest neighbour resampling.");
		}
		plhs[0] = mxCreateDoubleMatrix(m,n,mxREAL);
		plhs[1] = mxCreateDoubleMatrix(m,n,mxREAL);
		plhs[2] = mxCreateDoubleMatrix(m,n,mxREAL);
		plhs[3] = mxCreateDoubleMatrix(m,n,mxREAL);

		resample_d(m*n, map, mxGetPr(plhs[0]),mxGetPr(plhs[1]),mxGetPr(plhs[2]),mxGetPr(plhs[3]),
			mxGetPr(prhs[1]),mxGetPr(prhs[2]),mxGetPr(prhs[3]),
			hold, background);
	}
	free_maps(map, 1);
}
