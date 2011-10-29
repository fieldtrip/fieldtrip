/*
 * $Id$
 * John Ashburner
 */

#include "mex.h"
#include "spm_mapping.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	MAPTYPE *map, *get_maps();
	int m,n, k, hold, status;
	double *mat, *ptr, *img, background=0.0;

	if (nrhs != 4 || nlhs > 1)
	{
		mexErrMsgTxt("Incorrect usage.");
	}

	map = get_maps(prhs[0], &n);
	if (n!=1)
	{
		free_maps(map, n);
		mexErrMsgTxt("Bad image handle dimensions.");
	}

	for(k=1; k<=3; k++)
	{
		if (!mxIsNumeric(prhs[k]) || mxIsComplex(prhs[k]) ||
			mxIsSparse(prhs[k]) || !mxIsDouble(prhs[k]))
		{
			free_maps(map, 1);
			mexErrMsgTxt("Arguments must be numeric, real, full and double.");
		}
	}

	/* get transformation matrix */
	if (mxGetM(prhs[1]) != 4 && mxGetN(prhs[1]) != 4)
	{
		free_maps(map, 1);
		mexErrMsgTxt("Transformation matrix must be 4 x 4.");
	}
	mat = mxGetPr(prhs[1]);

	/* get output dimensions */
	if (mxGetM(prhs[2]) * mxGetN(prhs[2]) != 2)
	{
		free_maps(map, 1);
		mexErrMsgTxt("Output dimensions must have two elements.");
	}
	ptr = mxGetPr(prhs[2]);
	m = abs((int)ptr[0]);
	n = abs((int)ptr[1]);
	plhs[0] = mxCreateDoubleMatrix(m,n,mxREAL);
	img = mxGetPr(plhs[0]);

	if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 1 && mxGetM(prhs[3])*mxGetN(prhs[3]) != 2)
	{
		free_maps(map, 1);
		mexErrMsgTxt("Hold & background argument must have one or two element(s).");
	}
	hold = (int)(*(mxGetPr(prhs[3])));
	if (abs(hold) > 127)
	{
		free_maps(map, 1);
		mexErrMsgTxt("Bad hold value.");
	}

	if (mxGetM(prhs[3])*mxGetN(prhs[3]) > 1)
		background = mxGetPr(prhs[3])[1];

	status = slice(mat, img, m, n, map, hold, background);
	free_maps(map, 1);
	if (status)
	{
		mexErrMsgTxt("Slicing failed.");
	}
}
