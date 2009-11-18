#include "mex.h"

mxArray *mxCreateNumericArrayE(mwSize ndim, const mwSize *dims, 
			       mxClassID classid, mxComplexity ComplexFlag);
mxArray *mxCreateNumericMatrixE(mwSize m, mwSize n, mxClassID classid, 
				mxComplexity ComplexFlag);
mxArray *mxCreateDoubleMatrixE(mwSize m, mwSize n, 
			       mxComplexity ComplexFlag);
