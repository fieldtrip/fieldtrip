/* TRCOV     Evaluate covariance matrix.
 *
 *         Description
 *         C = GPEXPTRCOV(GPCF, TX) takes in Gaussian process GP and
 *         matrix TX that contains input vectors to GP. Returns
 *         covariance matrix C. Every element ij of C contains covariance
 *         between inputs i and j in TX.
 *
 *
 * Last modified: 2011-09-07 10:22:23 EEST
 *
 */

/* Copyright (C) 1998-2001 Aki Vehtari
 * Copyright (C) 2008      Jarno Vanhatalo
 *
 * This software is distributed under the GNU General Public
 * License (version 3 or later); please refer to the file
 * License.txt, included with the software, for details.
 *
 */

#include<string.h>
#include <stdlib.h>
#include <math.h>
#include "mex.h"
#include "matrix.h"
#define max(a, b) (((a) > (b)) ? (a) : (b))
#define PI (3.141592653589793)
void cumsum2(mwIndex *p, mwIndex *c, mwIndex n);

void mexFunction(const int nlhs, mxArray *plhs[],
        const int nrhs, const mxArray *prhs[]) {
  
  if (nlhs!=1 && nlhs!=3)
    mexErrMsgTxt( "Wrong number of output arguments.");
  
  if (nrhs!=2)
    mexErrMsgTxt( "Wrong number of input arguments." );
  
  {
    double *x, lms, ms, *l, rr, *rr2, *C, d, eps, c, *Ct, c1, c2, c3, D, D2;
    double percent_sparse, alpha, *period, decay, *cc, *s_sexp, pp, *s_sexp2;
    const mwSize *dims;
    const mxArray *selectedVariables;
    char *type;
    mwIndex i, j, k, ind, *I, *J, *It, *Jt, *Jc, *w2, *w;
    mwSize m, n, nnz, lr, pr;
    mxArray *field;
    
    dims = mxGetDimensions(prhs[1]);
    x = mxGetPr(prhs[1]);
    m = dims[0];
    
    if((selectedVariables=mxGetField(*prhs, 0, "selectedVariables"))!=NULL) {
      dims = mxGetDimensions(selectedVariables);
      if (dims[0]!=1 && dims[1]!=1)
        mexErrMsgTxt( "gpcf.selectedVariables must be a vector." );
      n = max(dims[0], dims[1]);
      cc = mxGetPr(selectedVariables);
    }else{
      n = dims[1];
      cc = NULL;
    }
    
    if((field=mxGetField(*prhs, 0, "magnSigma2"))==NULL)
      mexErrMsgTxt("Could not get gpcf.magnSigma2");
    dims = mxGetDimensions(field);
    if (dims[0]!=1 || dims[1]!=1)
      mexErrMsgTxt( "gpcf.magnSigma2 must be a scalar." );
    lms = log(mxGetScalar(field));
    ms = exp(lms);
    
    if((field=mxGetField(*prhs, 0, "lengthScale"))==NULL)
      mexErrMsgTxt("Could not get gpcf.lengthScale");
    dims = mxGetDimensions(field);
    if (dims[0]!=1 && dims[1]!=1)
      mexErrMsgTxt( "gpcf.lengthScale must be a scalar or a vector." );
    lr=max(dims[0], dims[1]);
    if (lr!=1 && lr!=n            )
      mexErrMsgTxt( "gpcf.lengthScale must be a scalar or its lenght must be the same as the number of columns in X." );
    l = mxGetPr(field);
    
    if((field=mxGetField(*prhs, 0, "type"))==NULL)
      mexErrMsgTxt("Could not get gpcf.type");
    if (mxIsChar(field) !=1)
      mexErrMsgTxt( "gpcf.type must be a string." );
    type = mxArrayToString(field);
    
    /*
     * squared exponential covariance
     */
    if( strcmp( type, "gpcf_sexp" ) == 0 ) {
      plhs[0]=mxCreateDoubleMatrix(m, m, mxREAL);
      C = mxGetPr(plhs[0]);
      eps=mxGetEps();
      
      for (i=0;i<n;i++) {
        rr=(lr>1)?(l[i]*l[i]):(l[0]*l[0]);   
        for (j=0;j<m;j++) {
          for (k=0;k<j;k++) {
            if(cc!=NULL) {
              d=x[j+m*((mwIndex)cc[i]-1)]-x[k+m*((mwIndex)cc[i]-1)];
            }else{
              d=x[j]-x[k];
            }
            C[j*m+k]+=d*d/rr/2.0;
          }
        }
        if(cc==NULL)
          x+=m;
      }
      for (j=0;j<m;j++) {
        for (k=0;k<j;k++) {
          d=exp(lms-C[j*m+k]);
          d=(d>eps) ? d : 0;
          C[j*m+k]=d;
          C[j+k*m]=d;
        }
        C[j*(m+1)]=ms;
      }
    }
    /*
     * exponential covariance
     */
    else if(strcmp( type, "gpcf_exp" ) == 0 ){
      plhs[0]=mxCreateDoubleMatrix(m, m, mxREAL);
      C = mxGetPr(plhs[0]);
      eps=mxGetEps();
      
      for (i=0;i<n;i++) {
        rr=(lr>1)?(l[i]*l[i]):(l[0]*l[0]);
        for (j=0;j<m;j++) {
          for (k=0;k<j;k++) {
            if(cc!=NULL) {
              d=x[j+m*((mwIndex)cc[i]-1)]-x[k+m*((mwIndex)cc[i]-1)];
            }else{
              d=x[j]-x[k];
            }
            C[j*m+k]+=d*d/rr;
          }
        }
        if (cc==NULL)
          x+=m;
      }
      
      for (j=0;j<m;j++) {
        for (k=0;k<j;k++) {
          d=exp(lms-sqrt(C[j*m+k]));
          d=(d>eps) ? d : 0;
          C[j*m+k]=d;
          C[j+k*m]=d;
        }
        C[j*(m+1)]=ms;
      }
    }
    /*
     * matern nu = 3/2 covariance
     */
    else if(strcmp( type, "gpcf_matern32" ) == 0 ){
      plhs[0]=mxCreateDoubleMatrix(m, m, mxREAL);
      C = mxGetPr(plhs[0]);
      eps=mxGetEps();
      
      for (i=0;i<n;i++) {
        rr=(lr>1)?(l[i]*l[i]):(l[0]*l[0]);
        for (j=0;j<m;j++) {
          for (k=0;k<j;k++) {
            if(cc!=NULL) {
              d=x[j+m*((mwIndex)cc[i]-1)]-x[k+m*((mwIndex)cc[i]-1)];
            }else{
              d=x[j]-x[k];
            }
            C[j*m+k]+=d*d/rr;
          }
        }
        if(cc==NULL)
          x+=m;
      }
      for (j=0;j<m;j++) {
        for (k=0;k<j;k++) {
          c = sqrt(3.0*C[j*m+k]);
          d=(1+c)*exp(lms-c);
          d=(d>eps) ? d : 0;
          C[j*m+k]=d;
          C[j+k*m]=d;
        }
        C[j*(m+1)]=ms;
      }
    }
    /*
     * matern nu = 5/2 covariance
     */
    else if(strcmp( type, "gpcf_matern52" ) == 0 ){
      plhs[0]=mxCreateDoubleMatrix(m, m, mxREAL);
      C = mxGetPr(plhs[0]);
      eps=mxGetEps();
      
      for (i=0;i<n;i++) {
        rr=(lr>1)?(l[i]*l[i]):(l[0]*l[0]);
        for (j=0;j<m;j++) {
          for (k=0;k<j;k++) {
            if(cc!=NULL) {
              d=x[j+m*((mwIndex)cc[i]-1)]-x[k+m*((mwIndex)cc[i]-1)];
            }else{
              d=x[j]-x[k];
            }
            C[j*m+k]+=d*d/rr;
          }
        }
        if(cc==NULL)
          x+=m;
      }
      for (j=0;j<m;j++) {
        for (k=0;k<j;k++) {
          c = sqrt(5.0*C[j*m+k]);
          d=(1+c+5.0*C[j*m+k]/3.0)*exp(lms-c);
          d=(d>eps) ? d : 0;
          C[j*m+k]=d;
          C[j+k*m]=d;
        }
        C[j*(m+1)]=ms;
      }
    }
    /*
     * piece wise polynomial 0 covariance
     */
    else if(strcmp( type, "gpcf_ppcs0" ) == 0 ){
      if((field=mxGetField(*prhs, 0, "l"))==NULL)
        mexErrMsgTxt("Could not get gpcf.l");
      dims = mxGetDimensions(field);
      if (dims[0]!=1 || dims[1]!=1)
        mexErrMsgTxt( "gpcf.l must be a scalar." );
      D = mxGetScalar(field);
      
      percent_sparse = 0.05;
      nnz = (mwSize) max(1, floor(percent_sparse*(double)m*(double)m));
      It = mxCalloc(nnz, sizeof(mwIndex));
      Jt = mxCalloc(nnz, sizeof(mwIndex));
      Ct = mxCalloc(nnz, sizeof(double));
      ind = 0;
      
      /* Set the length-scales in vector of length of number of inputs */
      rr2 = mxCalloc(n, sizeof(double));
      for (i=0;i<n;i++) {
        rr2[i]=(lr>1)?(l[i]*l[i]):(l[0]*l[0]);
      }
      
      /* Evaluate the distances that are less than one,
       * and evaluate the covariance at them.
       * This is strictly upper triangular matrix */
      D2 = D;
      for (j=0;j<m;j++) {
        for (k=0;k!=j;k++) {
          c = 0.0;
          for (i=0;i<n;i++) {
            if(cc!=NULL) {
              d=x[j+m*((mwIndex)cc[i]-1)]-x[k+m*((mwIndex)cc[i]-1)];
            }else{
              d=x[j+i*m]-x[k+i*m];
            }
            c+=d*d/rr2[i];
          }
          if (c<1.0){   /* store the covariance */
            if (ind==nnz){ /* allocate more memory */
              nnz=(mwSize)2*nnz;
              It = mxRealloc(It, nnz*sizeof(mwIndex));
              Jt = mxRealloc(Jt, nnz*sizeof(mwIndex));
              Ct = mxRealloc(Ct, nnz*sizeof(double));
            }
            d = ms*pow(1.0-sqrt(c), D2);
            It[ind] = k;
            Jt[ind] = j;
            Ct[ind] = d;
            ind++;
          }
        }
      }
      
      /* resize the vectors */
      It = mxRealloc(It, (mwSize)ind*sizeof(mwIndex));
      Jt = mxRealloc(Jt, (mwSize)ind*sizeof(mwIndex));
      Ct = mxRealloc(Ct, (mwSize)ind*sizeof(double));
      
      /* evaluate the row and column counts */
      w = mxCalloc(m, sizeof(mwIndex));              /* workspace */
      w2 = mxCalloc(m, sizeof(mwIndex));             /* workspace */
      for (k=0;k<ind;k++) w[It[k]]++;               /* row counts of the upper triangular */
      for (k=0;k<ind;k++) w2[Jt[k]]++;              /* column counts of the upper triangular */
      for (k=0;k<m;k++) w[k] += w2[k] + (mwIndex)1; /* column counts of the covariance matrix */
      Jc = mxCalloc(m+(mwSize)1, sizeof(mwIndex));
      cumsum2(Jc, w2, m);                           /* column starting points of the upper triangle */
      
      /* Create sparse matrix. Note! The matrix can contain only real numbers  */
      nnz = (mwSize)2*(mwSize)ind+m;
      plhs[0] = mxCreateSparse(m, m, nnz, mxREAL);
      I = mxGetIr(plhs[0]);
      J = mxGetJc(plhs[0]);
      C = mxGetPr(plhs[0]);
      
      /* Set the elements in the sparse matrix */
      cumsum2(J, w, m);                      /* column starting points */
      for (j = 0 ; j < m ; j++){             /* fill the upper triangular */
        for (k = Jc[j] ; k != Jc[j+1] ; k++){
          I[i = w[j]++] = It[k] ;
          if (C) C[i] = Ct[k] ;
        }
      }
      for (j = 0 ; j < m ; j++){             /* fill the diagonal */
        I[i = w[j]++] = j ;
        if (C) C[i] = ms;
      }
      for (j = 0 ; j < m ; j++){             /* fill the lower triangular */
        for (k = Jc[j] ; k != Jc[j+1] ; k++){
          I[i = w[ It[k]]++] = j ;
          if (C) C[i] = Ct[k] ;
        }
      }
      
      mxFree(It);
      mxFree(Jt);
      mxFree(Jc);
      mxFree(Ct);
      mxFree(rr2);
      mxFree(w);
      mxFree(w2);
      
    }
    /*
     * piece wise polynomial 1 covariance
     */
    else if(strcmp( type, "gpcf_ppcs1" ) == 0 ){
      if((field=mxGetField(*prhs, 0, "l"))==NULL)
        mexErrMsgTxt("Could not get gpcf.l");
      dims = mxGetDimensions(field);
      if (dims[0]!=1 || dims[1]!=1)
        mexErrMsgTxt( "gpcf.l must be a scalar." );
      D = mxGetScalar(field);
      
      percent_sparse = 0.05;
      nnz = (mwSize) max(1, floor(percent_sparse*(double)m*(double)m));
      It = mxCalloc(nnz, sizeof(mwIndex));
      Jt = mxCalloc(nnz, sizeof(mwIndex));
      Ct = mxCalloc(nnz, sizeof(double));
      ind = 0;
      
      /* Set the length-scales in vector of length of number of inputs */
      rr2 = mxCalloc(n, sizeof(double));
      for (i=0;i<n;i++) {
        rr2[i]=(lr>1)?(l[i]*l[i]):(l[0]*l[0]);
      }
      
      /* Evaluate the distances that are less than one,
       * and evaluate the covariance at them.
       * This is strictly upper triangular matrix */
      c1 = D + 1;
      D2 = D+1.0;
      for (j=0;j<m;j++) {
        for (k=0;k!=j;k++) {
          c = 0.0;
          for (i=0;i<n;i++) {
            if(cc!=NULL) {
              d=x[j+m*((mwIndex)cc[i]-1)]-x[k+m*((mwIndex)cc[i]-1)];
            }else{
              d=x[j+i*m]-x[k+i*m];
            }
            c+=d*d/rr2[i];
          }
          if (c<1.0){   /* store the covariance */
            if (ind==nnz){ /* allocate more memory */
              nnz=(mwSize)2*nnz;
              It = mxRealloc(It, nnz*sizeof(mwIndex));
              Jt = mxRealloc(Jt, nnz*sizeof(mwIndex));
              Ct = mxRealloc(Ct, nnz*sizeof(double));
            }
            d = c1*sqrt(c) + 1.0;
            d = ms*pow(1.0-sqrt(c), D2)*d;
            It[ind] = k;
            Jt[ind] = j;
            Ct[ind] = d;
            ind++;
          }
        }
      }
      
      /* resize the vectors */
      It = mxRealloc(It, (mwSize)ind*sizeof(mwIndex));
      Jt = mxRealloc(Jt, (mwSize)ind*sizeof(mwIndex));
      Ct = mxRealloc(Ct, (mwSize)ind*sizeof(double));
      
      /* evaluate the row and column counts */
      w = mxCalloc(m, sizeof(mwIndex));              /* workspace */
      w2 = mxCalloc(m, sizeof(mwIndex));             /* workspace */
      for (k=0;k<ind;k++) w[It[k]]++;               /* row counts of the upper triangular */
      for (k=0;k<ind;k++) w2[Jt[k]]++;              /* column counts of the upper triangular */
      for (k=0;k<m;k++) w[k] += w2[k] + (mwIndex)1; /* column counts of the covariance matrix */
      Jc = mxCalloc(m+(mwSize)1, sizeof(mwIndex));
      cumsum2(Jc, w2, m);                           /* column starting points of the upper triangle */
      
      /* Create sparse matrix. Note! The matrix can contain only real numbers  */
      nnz = (mwSize)2*(mwSize)ind+m;
      plhs[0] = mxCreateSparse(m, m, nnz, mxREAL);
      I = mxGetIr(plhs[0]);
      J = mxGetJc(plhs[0]);
      C = mxGetPr(plhs[0]);
      
      /* Set the elements in the sparse matrix */
      cumsum2(J, w, m);                      /* column starting points */
      for (j = 0 ; j < m ; j++){             /* fill the upper triangular */
        for (k = Jc[j] ; k != Jc[j+1] ; k++){
          I[i = w[j]++] = It[k] ;
          if (C) C[i] = Ct[k] ;
        }
      }
      for (j = 0 ; j < m ; j++){             /* fill the diagonal */
        I[i = w[j]++] = j ;
        if (C) C[i] = ms;
      }
      for (j = 0 ; j < m ; j++){             /* fill the lower triangular */
        for (k = Jc[j] ; k != Jc[j+1] ; k++){
          I[i = w[ It[k]]++] = j ;
          if (C) C[i] = Ct[k] ;
        }
      }
      
      mxFree(It);
      mxFree(Jt);
      mxFree(Jc);
      mxFree(Ct);
      mxFree(rr2);
      mxFree(w);
      mxFree(w2);
      
    }
    /*
     * piece wise polynomial 2 covariance
     */
    else if(strcmp( type, "gpcf_ppcs2" ) == 0 ){
      if((field=mxGetField(*prhs, 0, "l"))==NULL)
        mexErrMsgTxt("Could not get gpcf.l");
      dims = mxGetDimensions(field);
      if (dims[0]!=1 || dims[1]!=1)
        mexErrMsgTxt( "gpcf.l must be a scalar." );
      D = mxGetScalar(field);
      
      percent_sparse = 0.05;
      nnz = (mwSize) max(1, floor(percent_sparse*(double)m*(double)m));
      It = mxCalloc(nnz, sizeof(mwIndex));
      Jt = mxCalloc(nnz, sizeof(mwIndex));
      Ct = mxCalloc(nnz, sizeof(double));
      ind = 0;
      
      /* Set the length-scales in vector of length of number of inputs */
      rr2 = mxCalloc(n, sizeof(double));
      for (i=0;i<n;i++) {
        rr2[i]=(lr>1)?(l[i]*l[i]):(l[0]*l[0]);
      }
      
      /* Evaluate the distances that are less than one,
       * and evaluate the covariance at them.
       * This is strictly upper triangular matrix */
      c1 = D*D + 4.0*D + 3.0;
      c2 = 3.0*D + 6.0;
      D2 = D+2.0;
      for (j=0;j<m;j++) {
        for (k=0;k!=j;k++) {
          c = 0.0;
          for (i=0;i<n;i++) {
            if(cc!=NULL) {
              d=x[j+m*((mwIndex)cc[i]-1)]-x[k+m*((mwIndex)cc[i]-1)];
            }else{
              d=x[j+i*m]-x[k+i*m];
            }
            c+=d*d/rr2[i];
          }
          if (c<1.0){   /* store the covariance */
            if (ind==nnz){ /* allocate more memory */
              nnz=(mwSize)2*nnz;
              It = mxRealloc(It, nnz*sizeof(mwIndex));
              Jt = mxRealloc(Jt, nnz*sizeof(mwIndex));
              Ct = mxRealloc(Ct, nnz*sizeof(double));
            }
            d = c1*c + c2*sqrt(c) + 3.0;
            d = ms*pow(1.0-sqrt(c), D2)*d/3.0;
            It[ind] = k;
            Jt[ind] = j;
            Ct[ind] = d;
            ind++;
          }
        }
      }
      
      /* resize the vectors */
      It = mxRealloc(It, (mwSize)ind*sizeof(mwIndex));
      Jt = mxRealloc(Jt, (mwSize)ind*sizeof(mwIndex));
      Ct = mxRealloc(Ct, (mwSize)ind*sizeof(double));
      
      /* evaluate the row and column counts */
      w = mxCalloc(m, sizeof(mwIndex));              /* workspace */
      w2 = mxCalloc(m, sizeof(mwIndex));             /* workspace */
      for (k=0;k<ind;k++) w[It[k]]++;               /* row counts of the upper triangular */
      for (k=0;k<ind;k++) w2[Jt[k]]++;              /* column counts of the upper triangular */
      for (k=0;k<m;k++) w[k] += w2[k] + (mwIndex)1; /* column counts of the covariance matrix */
      Jc = mxCalloc(m+(mwSize)1, sizeof(mwIndex));
      cumsum2(Jc, w2, m);                           /* column starting points of the upper triangle */
      
      /* Create sparse matrix. Note! The matrix can contain only real numbers  */
      nnz = (mwSize)2*(mwSize)ind+m;
      plhs[0] = mxCreateSparse(m, m, nnz, mxREAL);
      I = mxGetIr(plhs[0]);
      J = mxGetJc(plhs[0]);
      C = mxGetPr(plhs[0]);
      
      /* Set the elements in the sparse matrix */
      cumsum2(J, w, m);                      /* column starting points */
      for (j = 0 ; j < m ; j++){             /* fill the upper triangular */
        for (k = Jc[j] ; k != Jc[j+1] ; k++){
          I[i = w[j]++] = It[k] ;
          if (C) C[i] = Ct[k] ;
        }
      }
      for (j = 0 ; j < m ; j++){             /* fill the diagonal */
        I[i = w[j]++] = j ;
        if (C) C[i] = ms;
      }
      for (j = 0 ; j < m ; j++){             /* fill the lower triangular */
        for (k = Jc[j] ; k != Jc[j+1] ; k++){
          I[i = w[ It[k]]++] = j ;
          if (C) C[i] = Ct[k] ;
        }
      }
      
      mxFree(It);
      mxFree(Jt);
      mxFree(Jc);
      mxFree(Ct);
      mxFree(rr2);
      mxFree(w);
      mxFree(w2);
      
    }
    /*
     * piece wise polynomial 3 covariance
     */
    else if(strcmp( type, "gpcf_ppcs3" ) == 0 ){
      if((field=mxGetField(*prhs, 0, "l"))==NULL)
        mexErrMsgTxt("Could not get gpcf.l");
      dims = mxGetDimensions(field);
      if (dims[0]!=1 || dims[1]!=1)
        mexErrMsgTxt( "gpcf.l must be a scalar." );
      D = mxGetScalar(field);
      
      percent_sparse = 0.05;
      nnz = (mwSize) max(1, floor(percent_sparse*(double)m*(double)m));
      It = mxCalloc(nnz, sizeof(mwIndex));
      Jt = mxCalloc(nnz, sizeof(mwIndex));
      Ct = mxCalloc(nnz, sizeof(double));
      ind = 0;
      
      /* Set the length-scales in vector of length of number of inputs */
      rr2 = mxCalloc(n, sizeof(double));
      for (i=0;i<n;i++) {
        rr2[i]=(lr>1)?(l[i]*l[i]):(l[0]*l[0]);
      }
      
      /* Evaluate the distances that are less than one,
       * and evaluate the covariance at them.
       * This is strictly upper triangular matrix */
      c1 = D*D*D + 9.0*D*D + 23.0*D + 15.0;
      c2 = 6.0*D*D + 36.0*D + 45.0;
      c3 = 15.0*D + 45.0;
      D2 = D+3.0;
      for (j=0;j<m;j++) {
        for (k=0;k!=j;k++) {
          c = 0.0;
          for (i=0;i<n;i++) {
            if(cc!=NULL) {
              d=x[j+m*((mwIndex)cc[i]-1)]-x[k+m*((mwIndex)cc[i]-1)];
            }else{
              d=x[j+i*m]-x[k+i*m];
            }
            c+=d*d/rr2[i];
          }
          if (c<1.0){   /* store the covariance */
            if (ind==nnz){ /* allocate more memory */
              nnz=(mwSize)2*nnz;
              It = mxRealloc(It, nnz*sizeof(mwIndex));
              Jt = mxRealloc(Jt, nnz*sizeof(mwIndex));
              Ct = mxRealloc(Ct, nnz*sizeof(double));
            }
            d = c1*c*sqrt(c) + c2*c + c3*sqrt(c) + 15.0;
            d = ms*pow(1.0-sqrt(c), D2)*d/15.0;
            It[ind] = k;
            Jt[ind] = j;
            Ct[ind] = d;
            ind++;
          }
        }
      }
      
      /* resize the vectors */
      It = mxRealloc(It, (mwSize)ind*sizeof(mwIndex));
      Jt = mxRealloc(Jt, (mwSize)ind*sizeof(mwIndex));
      Ct = mxRealloc(Ct, (mwSize)ind*sizeof(double));
      
      /* evaluate the row and column counts */
      w = mxCalloc(m, sizeof(mwIndex));              /* workspace */
      w2 = mxCalloc(m, sizeof(mwIndex));             /* workspace */
      for (k=0;k<ind;k++) w[It[k]]++;               /* row counts of the upper triangular */
      for (k=0;k<ind;k++) w2[Jt[k]]++;              /* column counts of the upper triangular */
      for (k=0;k<m;k++) w[k] += w2[k] + (mwIndex)1; /* column counts of the covariance matrix */
      Jc = mxCalloc(m+(mwSize)1, sizeof(mwIndex));
      cumsum2(Jc, w2, m);                           /* column starting points of the upper triangle */
      
      /* Create sparse matrix. Note! The matrix can contain only real numbers  */
      nnz = (mwSize)2*(mwSize)ind+m;
      plhs[0] = mxCreateSparse(m, m, nnz, mxREAL);
      I = mxGetIr(plhs[0]);
      J = mxGetJc(plhs[0]);
      C = mxGetPr(plhs[0]);
      
      /* Set the elements in the sparse matrix */
      cumsum2(J, w, m);                      /* column starting points */
      for (j = 0 ; j < m ; j++){             /* fill the upper triangular */
        for (k = Jc[j] ; k != Jc[j+1] ; k++){
          I[i = w[j]++] = It[k] ;
          if (C) C[i] = Ct[k] ;
        }
      }
      for (j = 0 ; j < m ; j++){             /* fill the diagonal */
        I[i = w[j]++] = j ;
        if (C) C[i] = ms;
      }
      for (j = 0 ; j < m ; j++){             /* fill the lower triangular */
        for (k = Jc[j] ; k != Jc[j+1] ; k++){
          I[i = w[ It[k]]++] = j ;
          if (C) C[i] = Ct[k] ;
        }
      }
      
      mxFree(It);
      mxFree(Jt);
      mxFree(Jc);
      mxFree(Ct);
      mxFree(rr2);
      mxFree(w);
      mxFree(w2);
      
    }
    /*
     * Rational quadratic covariance
     */
    else if( strcmp( type, "gpcf_rq" ) == 0 ) {
      if((field=mxGetField(*prhs, 0, "alpha"))==NULL)
        mexErrMsgTxt("Could not get gpcf.alpha");
      alpha = mxGetScalar(field);
      dims = mxGetDimensions(field);
      plhs[0]=mxCreateDoubleMatrix(m, m, mxREAL);
      C = mxGetPr(plhs[0]);
      eps=mxGetEps();
      for (i=0;i<n;i++) {
        rr=(lr>1)?(l[i]*l[i]):(l[0]*l[0]);    
        for (j=0;j<m;j++) {
          for (k=0;k<j;k++) {
            if(cc!=NULL) {
              d=x[j+m*((mwIndex)cc[i]-1)]-x[k+m*((mwIndex)cc[i]-1)];
            }else{
            d=x[j]-x[k];
            }
            C[j*m+k]+=d*d/(2.0*alpha*rr);
          }
        }
        if(cc==NULL)
          x+=m;
      }
      for (j=0;j<m;j++) {
        for (k=0;k<j;k++) {
          d=ms*pow((C[j*m+k]+1),-alpha);
          d=(d>eps) ? d : 0;
          C[j*m+k]=d;
          C[j+k*m]=d;
        }
        C[j*(m+1)]=ms;
      }
    }
    /*
     * Periodic covariance
     */
    else if( strcmp( type, "gpcf_periodic" ) == 0 ) {
      if((field=mxGetField(*prhs, 0, "period"))==NULL)
        mexErrMsgTxt("Could not get gpcf.period");
      period = mxGetPr(field);
      dims = mxGetDimensions(field);
      pr = max(dims[0], dims[1]);
      if((field=mxGetField(*prhs, 0, "decay"))==NULL)
        mexErrMsgTxt("Could not get gpcf.decay");
      decay = mxGetScalar(field);
      if(decay==1) {
        if((field=mxGetField(*prhs, 0, "lengthScale_sexp"))==NULL)
          mexErrMsgTxt("Could not get gpcf.lengthScale_sexp");
        dims = mxGetDimensions(field);
        s_sexp = mxGetPr(field);
        s_sexp2 = mxCalloc(max(dims[0], dims[1]), sizeof(s_sexp));
        if(max(dims[0], dims[1])>1) {
          for (i=0;i<max(dims[0], dims[1]);i++) {
            s_sexp2[i] = 1.0/pow(s_sexp[i],2);
          }
        } else {
          s_sexp2[0] = (1.0/s_sexp[0]);
          s_sexp2[0] = pow(s_sexp2[0],2);
        }
      }
      plhs[0]=mxCreateDoubleMatrix(m, m, mxREAL);
      C = mxGetPr(plhs[0]);
      eps=mxGetEps();
      for (i=0;i<n;i++, x+=m) {
        rr=(lr>1)?(l[i]*l[i]):(l[0]*l[0]);
        pp=(pr>1)?(period[i]):(period[0]);
        for (j=0;j<m;j++) {
          for (k=0;k<j;k++) {
            d=sin(PI*(x[j]-x[k])/pp);
            d=2*d*d/rr;
            if(decay==1 && max(dims[0], dims[1])>1) {
              d+=s_sexp2[j]/2*((x[j]-x[k])*(x[j]-x[k]));
            } else if(decay==1 && max(dims[0], dims[1]) == 1) {
              d+=s_sexp2[0]*pow((x[j]-x[k]),2)/2;
            }
            C[j*m+k]+=d;
          }
        } 
      }
      if(decay==1)
        mxFree(s_sexp2);
      for (j=0;j<m;j++) {
        for (k=0;k<j;k++) {
          d=exp(lms-C[j*m+k]);
          d=(d>eps) ? d : 0;
          C[j*m+k]=d;
          C[j+k*m]=d;
        }
        C[j*(m+1)]=ms;
      }
    }
    else{
      mexErrMsgTxt( "Undefined type of covariance function." );
    }
  }
    return;
}

void cumsum2(mwIndex *p, mwIndex *c, mwIndex n) {
  mwIndex i;
  mwIndex nz = 0;
  if(!p || !c) return;
  for (i=0;i<n;i++){
    p[i]=nz;
    nz+=c[i];
    c[i]=p[i];
  }
  p[n]=nz;
  
}
