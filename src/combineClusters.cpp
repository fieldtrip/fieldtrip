#include "mex.h"

#include <algorithm>
#include <vector>

/* The computational routine */
void combineClusters_impl(unsigned int *labelmat, unsigned int *total, mwSize spatdimlength, mwSize timefreqlength, mxLogical *neighbours, unsigned int *out) {
    
    /* increase the total by one because indices in labelmat are 1-based and our array is zero-based */
    (*total)++;
    
    /* fill array with 1:total */
    unsigned int *replaceby;
    replaceby = (unsigned int *) malloc( (*total) * sizeof(unsigned int));
    int n;
    for (n = 0; n < *total; n++) {
        replaceby[n] = n;
    }

    mwSize i;
    mwSize j;
    mwSize k;

    /* iterate over channels */
    for (i = 0; i < spatdimlength; i++) {

        /* iterate over possible neighbours for this channel */
        for (j = 0; j < spatdimlength; j++) {
        
            if ( *(neighbours + i*spatdimlength + j) ) {
                /* channel is a neighbour */

                for (k = 0; k < timefreqlength; k++) {
                    unsigned int a = *(labelmat + k*spatdimlength + i);
                    unsigned int b = *(labelmat + k*spatdimlength + j);
                    if (a > 0 && b > 0) {
                        if (replaceby[a] == replaceby[b]) {
                            continue;
                        } else if (replaceby[a] < replaceby[b]) {
                            for (n = 0; n < *total; n++) {
                                if (replaceby[n] == replaceby[b]) {
                                    replaceby[n] = replaceby[a];
                                }
                            }
                        } else if (replaceby[b] < replaceby[a]) {
                            for (n = 0; n < *total; n++) {
                                if (replaceby[n] == replaceby[a]) {
                                    replaceby[n] = replaceby[b];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    /* copy and sort replaceby, retain only unique elements */
    std::vector<int> replacebySorted (replaceby, replaceby+*total);
    std::sort(replacebySorted.begin(), replacebySorted.end());
    std::vector<int>::iterator it;
    it = std::unique(replacebySorted.begin(), replacebySorted.end());
    replacebySorted.resize(std::distance(replacebySorted.begin(), it));
    
    /* generate sequential cluster numbers */
    unsigned int *clusternums;
    unsigned int sortedSize = replacebySorted.size();

    clusternums = (unsigned int *) malloc(sortedSize * sizeof(unsigned int));
    for (n = 0; n < sortedSize; n++) {
        clusternums[n] = n; /* the first element will be 0 (i.e. no cluster present), as is the case in replacebySorted */
    }
    
    for (i = 0; i < spatdimlength; i++) {
        for (j = 0; j < timefreqlength; j++) {
            unsigned int val = *(labelmat + j*spatdimlength + i);
            if (val > 0) {
            
                /* look for the cluster number */
                for (n = 0; n < sortedSize; n++) {
                    if (replacebySorted[n] == replaceby[val]) {
                        *(out + j*spatdimlength + i) = clusternums[n];
                        break;
                    }
                }
            }
        }
    }
    
    /* clean up temp variables */
    free((void *) replaceby);
    free((void *) clusternums);
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    unsigned int *labelmat;
    unsigned int *total;
    mwSize spatdimlength;
    mwSize timefreqlength;
    mxLogical *neighbours;
    unsigned int *out;
    unsigned int *replaceby;

    /* check for proper number of arguments */
    if(nrhs != 3) {
        mexErrMsgIdAndTxt("FieldTrip:findcluster:nrhs", "three inputs required: labelmat, neighbours, total");
    }
    if(nlhs != 1) {
        mexErrMsgIdAndTxt("FieldTrip:findcluster:nlhs", "one output required");
    }
    
    /* make sure the first input argument is real uint32 */
    if (!mxIsUint32(prhs[0]) || mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt("FieldTrip:findcluster:notUint32", "first input must be a matrix of uint32");
    }

    if (!mxIsLogical(prhs[1])) {
        mexErrMsgIdAndTxt("FieldTrip:findcluster:notLogical", "second input must be logical matrix");
    }
    
    if (!mxIsUint32(prhs[2]) || mxGetNumberOfElements(prhs[2]) > 1) {
        mexErrMsgIdAndTxt("FieldTrip:findcluster:notUint32", "third input must be a scalar of uint32");
    }
    
    /* get dimensions of labelmat */
    spatdimlength = (mwSize)mxGetM(prhs[0]);
    timefreqlength = (mwSize)mxGetN(prhs[0]);
    
    /* perform dimension check */
    if (spatdimlength != mxGetM(prhs[1]) || spatdimlength != mxGetN(prhs[1])) {
        mexErrMsgIdAndTxt("FieldTrip:findcluster:neighboursNotOK","second input must be square matrix with one row and column for each channel");
    }
    
    /* get the other inputs */
    labelmat = (unsigned int *)mxGetData(prhs[0]);
    neighbours = (mxLogical *)mxGetData(prhs[1]);
    total = (unsigned int *)mxGetData(prhs[2]);

    /* create the output matrix */
    plhs[0] = mxCreateNumericMatrix(spatdimlength, timefreqlength, mxUINT32_CLASS, mxREAL);
    out = (unsigned int *)mxGetData(plhs[0]);

    /* call the computational routine */
    combineClusters_impl(labelmat, total, spatdimlength, timefreqlength, neighbours, out);
}
