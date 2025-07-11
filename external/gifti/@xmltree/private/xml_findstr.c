#include "mex.h"

/*
 * $Id: xml_findstr.c 6480 2015-06-13 01:08:30Z guillaume $
 * Guillaume Flandin <guillaume@artefact.tk>
 */

/*
    Differences with built-in findstr:
        - allows to search only the n first occurences of a pattern
        - allows to search only in a substring (given an index of the beginning)
   
    MATLAB hack:
        - doesn't use mxGetString to prevent a copy of the string.
        - assumes MATLAB stores strings as unsigned short (Unicode 16 bits)
          matrix.h: typedef uint16_T mxChar;
          (that's the case for MATLAB 5.*, 6.* and 7.* but MATLAB 4.* stores
           strings as double and GNU Octave as char, see src/mxarray.h)
*/

/* Comment the following line to use standard mxGetString (slower) */
#if !defined (HAVE_OCTAVE)
#define __HACK_MXCHAR__
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    unsigned int i, j, stext, spattern, nbmatch = 0, ind = 1, occur = 0, nboccur = 0;
#ifdef __HACK_MXCHAR__
    unsigned short int *text = NULL, *pattern = NULL;
#else
    char *text = NULL, *pattern = NULL;
#endif
    unsigned int *k = NULL;
    mxArray *out = NULL;
    
    /* Check for proper number of arguments. */
    if ((nrhs == 0) || (nrhs == 1))
        mexErrMsgTxt("Not enough input arguments.");
    else if (nrhs > 4)
        mexErrMsgTxt("Too many input arguments.");
    else if (nlhs > 1)
        mexErrMsgTxt("Too many output arguments.");
    
    /* The input TEXT must be a string */
    if (!mxIsChar(prhs[0]))
        mexErrMsgTxt("Inputs must be character arrays.");
    stext = mxGetM(prhs[0]) * mxGetN(prhs[0]);
#ifdef __HACK_MXCHAR__
    text = mxGetData(prhs[0]);
#else
    text = mxCalloc(stext+1, sizeof(char));
    mxGetString(prhs[0], text, stext+1);
#endif
        
    /* The input PATTERN must be a string */
    if (!mxIsChar(prhs[1]))
        mexErrMsgTxt("Inputs must be character arrays.");
    spattern = mxGetM(prhs[1]) * mxGetN(prhs[1]);
#ifdef __HACK_MXCHAR__
    pattern = mxGetData(prhs[1]);
#else
    pattern = mxCalloc(spattern+1, sizeof(char));
    mxGetString(prhs[1], pattern, spattern+1);
#endif

    /* The input INDEX must be an integer */
    if (nrhs > 2) {
        if ((!mxIsNumeric(prhs[2]) || (mxGetM(prhs[2]) * mxGetN(prhs[2]) !=  1)))
            mexErrMsgTxt("Index input must be an integer.");
        ind = (unsigned int)mxGetScalar(prhs[2]);
        if (ind < 1)
            mexErrMsgTxt("Index must be greater than 1.");
    }
    
    /* The input OCCUR must be an integer */
    if (nrhs == 4) {
        if ((!mxIsNumeric(prhs[3]) || (mxGetM(prhs[3]) * mxGetN(prhs[3]) !=  1)))
            mexErrMsgTxt("Index input must be an integer.");
        nboccur = (unsigned int)mxGetScalar(prhs[3]);
    }
    
    /* Find pattern in text */
    for (i=ind-1;i<stext;i++) {
        for (j=0;j<spattern && i+j<stext;j++) {
            if (pattern[j] == text[i+j]) {
                if (j == spattern-1) {
                    nbmatch += 1;
                    k = mxRealloc(k,nbmatch*sizeof(unsigned int));
                    k[nbmatch-1] = i+1;
                    if (++occur == nboccur) i = stext;
                }
            }
            else break;
        }
    }
    
    /* Allocate output */
    out = mxCreateDoubleMatrix((nbmatch) ? 1:0, nbmatch, mxREAL);
    
    /* Copy index array into output */
    for (i=0;i<nbmatch;i++)
        mxGetPr(out)[i] = (double)k[i];
    
    /* Assign pointer to output */
    plhs[0] = out;
    
    /* Free memory */
    if (k) mxFree(k);
}
