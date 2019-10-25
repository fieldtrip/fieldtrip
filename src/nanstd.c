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
    double *inputr_p,  *inputi_p,  *output1r_p,  *cnt,  *ssqr,  *ssqi,  *sumi,  biasterm;
    
    /* used with single precision input */
    float  *inputr_ps, *inputi_ps, *output1r_ps, *cnts, *ssqrs, *ssqis, *sumis, biasterms;

    /* figure out the classid */
    if (nrhs>0) {
      classid = mxGetClassID(prhs[0]);
    } else {
      classid = mxUNKNOWN_CLASS;
    }

    /* check inputs */
    if (nrhs > 3) {
        mexErrMsgTxt("Too many input arguments.");
    } else if (nrhs == 3) {
        if (~mxIsEmpty(prhs[2]) && (mxGetM(prhs[2])!=1 || mxGetN(prhs[2])!=1)) {
            mexErrMsgTxt("Invalid dimension for input argument 2, must be scalar. Weights vector is not supported in this implementation of nanstd.");
        }
        if (mxGetScalar(prhs[2]) <= 0) {
            mexErrMsgTxt("Invalid value for input argument 2.");
        }
    } else if (nrhs < 1) {
        mexErrMsgTxt("Too few input arguments, at least 1 required.");
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

    /* extract dimension for standard deviation */
    if (nrhs==3) {
        dim = mxGetScalar(prhs[2]) - 1;
    } else {
        /* figure out the averaging dimension when only 1 input argument is given */
        dim = 0;
        for (i=0; i<numdims; i++) {
            if (dims[i]>1) {
                dim = i;
                break;
            }
        }
    }

    /* determine the normalisation term (either N or N-1), this depends on the numeric precision */
    if (classid==mxDOUBLE_CLASS) {
        if (nrhs==1) {
            biasterm = 1.0;
        } else if (mxIsEmpty(prhs[1])) {
            biasterm = 1.0;
        } else if (mxGetScalar(prhs[1])==1) {
            biasterm = 0.0;
        } else if (mxGetScalar(prhs[1])==0) {
            biasterm = 1.0;
        } else {
            mexErrMsgTxt("Invalid value for second input argument: this should either be [], 0, or 1.");
        }
    } else if (classid==mxSINGLE_CLASS) {
        if (nrhs==1) {
            biasterms = 1.0;
        } else if (mxIsEmpty(prhs[1])) {
            biasterms = 1.0;
        } else if (mxGetScalar(prhs[1])==1) {
            biasterms = 0.0;
        } else if (mxGetScalar(prhs[1])==0) {
            biasterms = 1.0;
        } else {
            mexErrMsgTxt("Invalid value for second input argument: this should either be [], 0, or 1.");
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
    if (numdims>=dim+1) {
        numelout = numelin / dims[dim];
    } else {
        numelout = numelin;
    }

    /* compute helper variables x1 and y1 needed for the indexing */
    if (dim+1>numdims)
        /* this essentially means that no averaging is done */
    {
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

    if (classid==mxDOUBLE_CLASS) {
        /* allocate memory for denominator */
        plhs[1] = mxCreateNumericArray((numdims+x0), dimsout, classid, mxREAL);
        cnt     = mxGetData(plhs[1]);

        /* associate inputs */
        inputr_p = mxGetData(prhs[0]);
        inputi_p = mxGetImagData(prhs[0]);

        /* assign the outputs */
        if (inputi_p == NULL) {
            plhs[0]    = mxCreateNumericArray((numdims+x0), dimsout, classid, mxREAL);
            output1r_p = mxGetData(plhs[0]);
            plhs[2]    = mxCreateNumericArray((numdims+x0), dimsout, classid, mxREAL);
            ssqr       = mxGetData(plhs[2]);
        } else {
            plhs[0]    = mxCreateNumericArray((numdims+x0), dimsout, classid, mxREAL);
            output1r_p = mxGetData(plhs[0]);
            plhs[2]    = mxCreateNumericArray((numdims+x0), dimsout, classid, mxCOMPLEX);
            ssqr       = mxGetData(plhs[2]);
            ssqi       = mxGetImagData(plhs[2]);
            plhs[3]    = mxCreateNumericArray((numdims+x0), dimsout, classid, mxREAL);
            sumi       = mxGetData(plhs[3]);
        }

        if (inputi_p == NULL) { /* input data is real-valued */
            
            /* compute variance using 'online' algorithm */
            /* Knuth (1998). The Art of Computer Programming, volume 2: Seminumerical Algorithms, 3rd edn., p. 232. */
            double delta;
            for (i=0; i<numelin; i++) {
                if (!isnan(inputr_p[i])) {
                    indx             = i%x1 + (i/y1) * x1;
                    cnt[indx]        = cnt[indx]        + 1.0;
                    
                    delta = inputr_p[i] - output1r_p[indx];
                    output1r_p[indx] = output1r_p[indx] + delta/cnt[indx];
                    ssqr[indx]       = ssqr[indx] + delta*(inputr_p[i] - output1r_p[indx]);
                    
                } else if (dim+1>numdims) {
                    output1r_p[i] = inputr_p[i];
                    ssqr[i]       = inputr_p[i]*inputr_p[i];
                    cnt[i]        = 1.0;
                }

            }
            
            /* compute the variance */
            for (i=0; i<numelout; i++) {
                output1r_p[i] = sqrt(ssqr[i]/(cnt[i] - biasterm));
            }

            
        } else
            /* handle the complex valued case separately */
        {
        
            /* compute variance using 'online' algorithm */
            /* Knuth (1998). The Art of Computer Programming, volume 2: Seminumerical Algorithms, 3rd edn., p. 232. */
            double delta;
            for (i=0; i<numelin; i++) {
                if (!isnan(inputr_p[i]) && !isnan(inputi_p[i])) {
                    indx             = i%x1 + (i/y1) * x1;
                    cnt[indx]        = cnt[indx]        + 1.0;
                    
                    delta = inputr_p[i] - output1r_p[indx];
                    output1r_p[indx] = output1r_p[indx] + delta/cnt[indx];
                    ssqr[indx] = ssqr[indx] + delta*(inputr_p[i] - output1r_p[indx]);
                    
                    delta = inputi_p[i] - sumi[indx];
                    sumi[indx] = sumi[indx] + delta/cnt[indx];
                    ssqi[indx] = ssqi[indx] + delta*(inputi_p[i] - sumi[indx]);
                    
                } else if (dim+1>numdims) {
                    output1r_p[i] = inputr_p[i];
                    sumi[i] = inputi_p[i];
                    ssqr[i]       = inputr_p[i]*inputr_p[i];
                    ssqi[i]       = inputi_p[i]*inputi_p[i];
                    cnt[i]        = 1.0;
                }

            }
            
            /* compute the variance */
            for (i=0; i<numelout; i++) {
                output1r_p[i] = sqrt((ssqr[i] + ssqi[i])/(cnt[i] - biasterm));
            }
        }

        /* free memory */
        mxFree(dimsout);

        return;
    }

    else if (classid==mxSINGLE_CLASS) {
        /* allocate memory for denominator */
        plhs[1] = mxCreateNumericArray((numdims+x0), dimsout, classid, mxREAL);
        cnts    = mxGetData(plhs[1]);

        /* associate inputs */
        inputr_ps = mxGetData(prhs[0]);
        inputi_ps = mxGetImagData(prhs[0]);

        /* assign the outputs */
        if (inputi_ps == NULL) {
            plhs[0]     = mxCreateNumericArray((numdims+x0), dimsout, classid, mxREAL);
            output1r_ps = mxGetData(plhs[0]);
            plhs[2]     = mxCreateNumericArray((numdims+x0), dimsout, classid, mxREAL);
            ssqrs       = mxGetData(plhs[2]);
        } else {
            plhs[0]     = mxCreateNumericArray((numdims+x0), dimsout, classid, mxREAL);
            output1r_ps = mxGetData(plhs[0]);
            plhs[2]     = mxCreateNumericArray((numdims+x0), dimsout, classid, mxCOMPLEX);
            ssqrs       = mxGetData(plhs[2]);
            ssqis       = mxGetImagData(plhs[2]);
            plhs[3]     = mxCreateNumericArray((numdims+x0), dimsout, classid, mxREAL);
            sumis       = mxGetData(plhs[3]);
        }

        if (inputi_ps == NULL) { /* real-valued input data */
        
            /* compute variance using 'online' algorithm */
            /* Knuth (1998). The Art of Computer Programming, volume 2: Seminumerical Algorithms, 3rd edn., p. 232. */
            double delta;
            for (i=0; i<numelin; i++) {
                if (!isnan(inputr_ps[i])) {
                    indx             = i%x1 + (i/y1) * x1;
                    cnts[indx]        = cnts[indx]        + 1.0;
                    
                    delta = inputr_ps[i] - output1r_ps[indx];
                    output1r_ps[indx] = output1r_ps[indx] + delta/cnts[indx];
                    ssqrs[indx]       = ssqrs[indx] + delta*(inputr_ps[i] - output1r_ps[indx]);
                    
                } else if (dim+1>numdims) {
                    output1r_ps[i] = inputr_ps[i];
                    ssqrs[i]       = inputr_ps[i]*inputr_ps[i];
                    cnts[i]        = 1.0;
                }

            }
            
            /* compute the variance */
            for (i=0; i<numelout; i++) {
                output1r_ps[i] = sqrt(ssqrs[i]/(cnts[i] - biasterms));
            }
            
        } else
            /* handle the complex valued case separately */
        {
        
            /* compute variance using 'online' algorithm */
            /* Knuth (1998). The Art of Computer Programming, volume 2: Seminumerical Algorithms, 3rd edn., p. 232. */
            double delta;
            for (i=0; i<numelin; i++) {
                if (!isnan(inputr_ps[i]) && !isnan(inputi_ps[i])) {
                    indx             = i%x1 + (i/y1) * x1;
                    cnts[indx]        = cnts[indx]        + 1.0;
                    
                    delta = inputr_ps[i] - output1r_ps[indx];
                    output1r_ps[indx] = output1r_ps[indx] + delta/cnts[indx];
                    ssqrs[indx] = ssqrs[indx] + delta*(inputr_ps[i] - output1r_ps[indx]);
                    
                    delta = inputi_ps[i] - sumis[indx];
                    sumis[indx] = sumis[indx] + delta/cnts[indx];
                    ssqis[indx] = ssqis[indx] + delta*(inputi_ps[i] - sumis[indx]);
                    
                } else if (dim+1>numdims) {
                    output1r_ps[i] = inputr_ps[i];
                    sumis[i] = inputi_ps[i];
                    ssqrs[i]       = inputr_ps[i]*inputr_ps[i];
                    ssqis[i]       = inputi_ps[i]*inputi_ps[i];
                    cnts[i]        = 1.0;
                }

            }
            
            /* compute the variance */
            for (i=0; i<numelout; i++) {
                output1r_ps[i] = sqrt((ssqrs[i] + ssqis[i])/(cnts[i] - biasterms));
            }
        }

        /* free memory */
        mxFree(dimsout);

        return;
    }

    else {
        /* we now at this point the input data is either numeric, char, or logical, but not double or single precision */
        /* since only double or single can be NaN, simply call matlab's std() function to do the work, we can safely ignore nans */
        mexCallMATLAB(nlhs, plhs, nrhs, (mxArray **)prhs, "std");
        return;
    }
}
