// $Id: gc_aux_mex.cpp 2 2009-06-16 19:24:10Z gramfort $
// $LastChangedBy: gramfort $
// $LastChangedDate: 2009-06-16 15:24:10 -0400 (Mar, 16 jui 2009) $
// $Revision: 2 $

#include "mex.h"
#include "maxflow/maxflow.h"

#include <limits>

void mexFunction( int  nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if  ( (nrhs  != 2) ) {
        mexErrMsgTxt("Fast Graph Cut: Wrong number of input arguments");
    }

    size_t K = mxGetM(prhs[0]);
    size_t T = mxGetN(prhs[0]);

    // mexPrintf("Nb trials : %d\n",K);
    // mexPrintf("Nb time samples : %d\n",T);

    if (K*T == 0)
    {
        mexErrMsgTxt("Trials cannot be empty !");
    }
    if (mxGetM(prhs[1])*mxGetN(prhs[1]) == 0)
    {
        mexErrMsgTxt("Lambda cannot be empty !");
    }

    double *trials = (double*)mxGetPr(prhs[0]);
    double *palpha = (double*)mxGetPr(prhs[1]);
    double alpha = palpha[0];

    // mexPrintf("Lambda : %f\n",alpha);

    // ============================
    // = start graph construction =
    // ============================
    typedef Graph<double,double,double> GraphType;
    GraphType *g = new GraphType(K*T,(K-1)*T + (T-1)*K);

    for(size_t i = 0; i < K*T; ++i) {
        g -> add_node();
    }

    for(size_t i = 0; i < K*T; ++i) {
        if(i < K) {
            g -> add_tweights( i,   std::numeric_limits<double>::max() / 2.0, 0.0 );
        } else if (i >= K*(T-1)) {
            g -> add_tweights( i,   0.0, std::numeric_limits<double>::max() / 2.0 );
        } else {
            g -> add_tweights( i,   0.0, 0.0 );
        }
    }

    for(size_t k = 0; k < K; ++k) {
        for(size_t t = 0; t < T-1; ++t) {
            g -> add_edge( t*K+k, (t+1)*K+k, trials[t*K+k], std::numeric_limits<double>::max() / 2.0 );
        }
    }

    for(size_t k = 0; k < K-1; ++k) {
        for(size_t t = 0; t < T; ++t) {
            g -> add_edge( t*K+k, t*K+k+1, std::numeric_limits<double>::max() / 2.0 , alpha);
        }
    }

    plhs[0] = mxCreateDoubleMatrix(K,1,mxREAL);
    double* lags = mxGetPr(plhs[0]);

    plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
    double* flow = mxGetPr(plhs[1]);

    *flow = g -> maxflow();

    for(size_t k = 0; k < K; ++k) {
        for(size_t t = 0; t < T-1; ++t) {
            if (g->what_segment(t*K+k) == GraphType::SOURCE && g->what_segment((t+1)*K+k) != GraphType::SOURCE) {
                lags[k] = (double)t + 1.0; // add 1 to match with matlab index
                // lags[k] = (double)t;
            }
        }
    }

    if(nlhs == 3) {
        plhs[2] = mxCreateDoubleMatrix(K,T,mxREAL);
        double* labels = mxGetPr(plhs[2]);
        for(size_t k = 0; k < K; ++k) {
            for(size_t t = 0; t < T; ++t) {
                if(g->what_segment(t*K+k) == GraphType::SOURCE) {
                    labels[t*K+k] = 1;
                } else {
                    labels[t*K+k] = 0;
                }
            }
        }
    }

    delete g;

    return;
}
