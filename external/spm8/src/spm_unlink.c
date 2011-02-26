/*
 * $Id: spm_unlink.c 938 2007-10-12 19:09:31Z john $
 * John Ashburner
 */

/* Do a silent deletion of files on disk */

#include <stdio.h>
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int i;
	if (nlhs != 0) mexErrMsgTxt("Too many output arguments.");

	for(i=0; i<nrhs; i++)
	{
		const mxArray *matptr = prhs[i];
		if (!mxIsNumeric(matptr))
		{
			char *str = NULL;
			int k, stlen;
			
			stlen = mxGetN(matptr);
			str = (char *)mxCalloc(stlen+1, sizeof(char));
			mxGetString(matptr,str,stlen+1);

			/* delete white space */
			for(k=0; k<stlen; k++)
				if (str[k] == ' ')
				{
					str[k] = '\0';
					break;
				}
			remove(str); /* not bothered about return status */
			mxFree(str);
		}
		else
			mexErrMsgTxt("Filename should be a string.");

	}
}
