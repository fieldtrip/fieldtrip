/*
 * $Id: spm_global.c 247 2005-10-04 17:20:34Z guillaume $
 */

/* returns the global mean for a memory mapped volume image
  FORMAT [G] = spm_global(V)
  V   - memory mapped volume
  G   - mean global activity
 ____________________________________________________________________________
 
  spm_global returns the mean counts integrated over all the  
  slices from the volume
 
  The mean is estimated after discounting voxels outside the object
  using a criteria of greater than > (global mean)/8
*/

#include <math.h>
#include "mex.h"
#include "spm_mapping.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int i, j, m, n;
	double s1=0.0, s2=0.0;
	double *dat;
	MAPTYPE *map, *get_maps();
	static double M[] = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};

	if (nrhs != 1 || nlhs > 1)
	{
		mexErrMsgTxt("Incorrect usage.");
	}

	map = get_maps(prhs[0], &j);

	if (j != 1)
	{
		free_maps(map, j);
		mexErrMsgTxt("Incorrect usage.");
	}
	n = (map->dim[0])*(map->dim[1]);
	dat = (double *)mxCalloc(n, sizeof(double));

	s1 = 0.0;
	m = 0;
	for (i=0; i<map->dim[2]; i++)
	{
		M[14] = i+1;
		slice(M, dat, map->dim[0],map->dim[1], map, 0,0);
		for(j=0;j<n; j++)
			if (mxIsFinite(dat[j]))
			{
				s1 += dat[j];
				m ++;
			}
	}
	s1/=(8.0*m);

	s2=0.0;
	m =0;
	for (i=0; i<map->dim[2]; i++)
	{
		M[14] = i+1;
		slice(M, dat, map->dim[0],map->dim[1], map, 0,0);
		for(j=0;j<n; j++)
			if (mxIsFinite(dat[j]) && dat[j]>s1)
			{
				m++;
				s2+=dat[j];
			}
	}
	s2/=m;

	plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
	mxGetPr(plhs[0])[0]=s2;
	free_maps(map, 1);
}
