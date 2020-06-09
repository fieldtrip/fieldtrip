#include <mex.h>
#include <stdlib.h>
#include <string.h>
#include "compiler.h"

#if defined (COMPILER_MSVC)
#include <math.h>
#define isnan _isnan
#define INFINITY (HUGE_VAL+HUGE_VAL)
#define NAN (INFINITY - INFINITY)
#elif defined(COMPILER_LCC)
#include <math.h>
#define INFINITY (DBL_MAX+DBL_MAX)
#define NAN (INFINITY - INFINITY)
#else
#include <math.h>
#endif


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    /* declare variables */
    const mwSize *dims;
    mwSize *dimsout;
    mwIndex indx;
    int i, numdims, dim;
    int numelin, numelout, x0, x1, y1;
    mxClassID classid;
    
    /* used with double precision input */
    double *inputr_p,  *inputi_p,  *output1r_p,  *output1i_p;
    
    /* used with single precision input */
    float  *inputr_ps, *inputi_ps, *output1r_ps, *output1i_ps;

    /* figure out the classid */
    if (nrhs>0) {
      classid = mxGetClassID(prhs[0]);
    } else {
      classid = mxUNKNOWN_CLASS;
    }

    /* check inputs */
    if (nrhs > 2) {
        mexErrMsgTxt("Too many input arguments, maximum of 2 supported.");
    } else if (nrhs < 1) {
        mexErrMsgTxt("Too few input arguments, at least 1 required.");
    } else if (nrhs == 2) {
        if (mxGetM(prhs[1])!=1 || mxGetN(prhs[1])!=1) {
            mexErrMsgTxt ("Invalid dimension for input argument 2, scalar required.");
        }
        if (mxGetScalar(prhs[1])<=0) {
            mexErrMsgTxt ("Invalid value for input argument 2, positive integer required.");
        }
    }

    if (mxIsEmpty(prhs[0])) {
        plhs[0] = mxCreateDoubleScalar(NAN);
        return;
    } else if (!mxIsNumeric(prhs[0]) && !mxIsLogical(prhs[0]) && !mxIsChar(prhs[0])) {
        mexErrMsgTxt ("Input argument 1 should be either numeric or logical.");
    }

    /* figure out dimension info and number of elements */
    dims    = mxGetDimensions(prhs[0]);
    numdims = mxGetNumberOfDimensions(prhs[0]);
    numelin = mxGetNumberOfElements(prhs[0]);

    if (nrhs==2) {
        dim = mxGetScalar(prhs[1]) - 1;
    } else if (nrhs==1) {
        /* figure out the averaging dimension when only 1 input argument is given */
        dim = 0;
        for (i=0; i<numdims; i++) {
            if (dims[i]>1) {
                dim = i;
                break;
            }
        }
    }

    /* helper variable needed to kick out the last dimension, if this is the averaging dimension */
    x0 = 0;
    if (numdims==dim+1) {
        x0 = -1;
    }

    /* create the vector which contains the dimensionality of the output */
    dimsout = mxMalloc((numdims+x0) * sizeof(mwSize));
    for (i=0; i<numdims+x0; i++) {
        dimsout[i] = dims[i];
    }

    /* make the dimension over which the averaging is done singleton in the output */
    if (numdims>dim+1) {
        dimsout[dim] = 1;
    }

    /* compute the number of output elements */
    if (numdims >= dim+1) {
        numelout = numelin / dims[dim];
    } else {
        numelout = numelin;
    }

    /* compute helper variables x1 and y1 needed for the indexing */
    if (dim+1 > numdims) {
        /* this essentially means that no averaging is done */
        x1 = numelin;
        y1 = numelin;
    } else {
        x1 = 1;
        for (i=0; i<numdims+x0; i++) {
            if (i==dim) {
                break;
            }

            x1 = x1 * dims[i];
        }
        y1 = x1 * dims[dim];
    }

    if (classid == mxDOUBLE_CLASS) {
        /* associate inputs */
        inputr_p = mxGetData(prhs[0]);
        inputi_p = mxGetImagData(prhs[0]);

        /* assign the outputs */
        if (inputi_p == NULL) {
            plhs[0]    = mxCreateNumericArray((numdims+x0), dimsout, classid, mxREAL);
            output1r_p = mxGetData(plhs[0]);
        } else {
            plhs[0]    = mxCreateNumericArray((numdims+x0), dimsout, classid, mxCOMPLEX);
            output1r_p = mxGetData(plhs[0]);
            output1i_p = mxGetImagData(plhs[0]);
        }

        if (inputi_p == NULL) {
            /* compute running sum */
            for (i=0; i<numelin; i++) {
                if (!isnan(inputr_p[i])) {
                    indx             = i%x1 + (i/y1) * x1;
                    output1r_p[indx] = output1r_p[indx] + inputr_p[i];
                } else if (dim+1 > numdims) {
                    /* we are summing across a singleton dimension, nothing */
                    /* should happen except all nans become 0 (as per Mathworks'  implementation) */
                    output1r_p[i] = 0.0;
                }

            }
        } else {
            /* handle the complex valued case separately */
            
            /* compute running sum */
            for (i=0; i<numelin; i++) {
                if (!isnan(inputr_p[i]) && !isnan(inputi_p[i])) {
                    indx             = i%x1 + (i/y1) * x1;
                    output1r_p[indx] = output1r_p[indx] + inputr_p[i];
                    output1i_p[indx] = output1i_p[indx] + inputi_p[i];
                } else if (dim+1>numdims) {
                    /* we are summing across a singleton dimension, nothing */
                    /* should happen except all nans become 0 (as per Mathworks'  implementation) */
                    output1r_p[i] = 0.0;
                    output1i_p[i] = 0.0;
                }
            }
        }

        /* free memory */
        mxFree(dimsout);

        return;
    }

    else if (classid == mxSINGLE_CLASS) {
        /* associate inputs */
        inputr_ps = mxGetData(prhs[0]);
        inputi_ps = mxGetImagData(prhs[0]);

        /* assign the outputs */
        if (inputi_ps == NULL) {
            plhs[0]     = mxCreateNumericArray((numdims+x0), dimsout, classid, mxREAL);
            output1r_ps = mxGetData(plhs[0]);
        } else {
            plhs[0]     = mxCreateNumericArray((numdims+x0), dimsout, classid, mxCOMPLEX);
            output1r_ps = mxGetData(plhs[0]);
            output1i_ps = mxGetImagData(plhs[0]);
        }

        if (inputi_ps == NULL) {
            /* compute running sum */
            for (i=0; i<numelin; i++) {
                if (!isnan(inputr_ps[i])) {
                    indx              = i%x1 + (i/y1) * x1;
                    output1r_ps[indx] = output1r_ps[indx] + inputr_ps[i];
                } else if (dim+1>numdims) {
                    /* we are summing across a singleton dimension, nothing */
                    /* should happen except all nans become 0 (as per Mathworks'  implementation) */
                    output1r_p[i] = 0.0;
                }
            }

        } else {
            /* handle the complex valued case separately */
            
            /* compute running sum */
            for (i=0; i<numelin; i++) {
                if (!isnan(inputr_ps[i]) && !isnan(inputi_ps[i])) {
                    indx              = i%x1 + (i/y1) * x1;
                    output1r_ps[indx] = output1r_ps[indx] + inputr_ps[i];
                    output1i_ps[indx] = output1i_ps[indx] + inputi_ps[i];
                } else if (dim+1>numdims) {
                    /* we are summing across a singleton dimension, nothing */
                    /* should happen except all nans become 0 (as per Mathworks'  implementation) */
                    output1r_ps[i] = 0.0;
                    output1i_ps[i] = 0.0;
                }
            }
        }

        /* free memory */
        mxFree(dimsout);

        return;
    } else {
        /* we now at this point the input data is either numeric, char, or logical, but not double or single precision */
        /* since only double or single can be NaN, simply call matlab's sum() function to do the work, we can safely ignore nans */
        mexCallMATLAB(nlhs, plhs, nrhs, (mxArray **)prhs, "sum");
        return;
    }
}

