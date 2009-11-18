#include "mex.h"

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]) {
    
    double *C;
    double *Ca;
    double scalar12;
    double scalar2;
    
    double *output;
    
    int i, j;
    int nrows, ncols;
    double val;
    
    
    /* Get matrix x */
    C = mxGetPr(prhs[0]);
    ncols = mxGetN(prhs[0]);
    nrows = mxGetM(prhs[0]);
    
    /* Get matrix y */
    Ca = mxGetPr(prhs[1]);
    scalar12 = mxGetScalar(prhs[2]);
    scalar2 = mxGetScalar(prhs[3]);
    
    
    /* create the output vector */
    plhs[0] = mxCreateDoubleMatrix(nrows, ncols, mxREAL);
    
    /* get a pointer to the real data in the output matrix */
    output = mxGetPr(plhs[0]);
    
    for(i=0;i<nrows;i++) {
        for(j=i;j<ncols;j++) {
            val = C[(i*ncols)+j] + Ca[i]*Ca[j]*scalar2 - Ca[i]*Ca[j]*scalar12;
            output[(i*ncols)+j] = val;
            output[i+(ncols*j)] = val;
        }
    }
}
