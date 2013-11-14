#include "mex.h"

/*
 * mex implementation of Rue's method.
 *
 * Marcel van Gerven, 2009
 */

static void
finv(mxArray *S, const mxArray *L, double *cirS, double *cjcS)
{
  double  *prL, *prS;
  mwIndex *irL, *jcL;
  mwIndex *irS, *jcS;
  mwIndex *diagL;
  mwSize  n, nzL, nzS, k, row, col, i, idxS, idxL, ngtcol, rmin;
  
  double *rowL;
   
  /* number of rows and columns in L */
  n = mxGetN(L);
  
  /* Get the starting positions of all data arrays in L */ 
  prL = mxGetPr(L);
  irL = mxGetIr(L);
  jcL = mxGetJc(L);

  /* get the number of non-zeros in L and S */
  nzL = jcL[n];
  nzS = mxGetNzmax(S);
 
  /* initialize the covariance matrix */
	mxSetPr(S, mxCalloc(nzS, sizeof(double)));  
  prS = mxGetPr(S);
 
  irS = mxCalloc(nzS, sizeof(mwIndex));  
  for (k=0;k<nzS;k++) {
    irS[k] = (mwIndex)cirS[k];
  }    
  mxSetIr(S, irS);  

  jcS = mxCalloc(n+1, sizeof(mwIndex));  
  for (k=0;k<=n;k++) {
    jcS[k] = (mwIndex)(cjcS[k]);
  }    
  mxSetJc(S, jcS);
     
  /* find indices belonging to diagonal elements */
  
  diagL = mxCalloc(n, sizeof(mwIndex));
  
  for (col=0; col<n; col++) {
        
    for (k=jcL[col]; k<jcL[col+1]; k++) {
      
      if (irL[k]==col) {
        diagL[col] = k;
        break;
      }
      
    }       
  }
  
  /* iterate over all elements using sparse column and row loops */

  /* create space for the L and S rows */
  rowL = mxCalloc(n, sizeof(double));
  
  /* iterate over columns */
  for (col=n-1; col>=0; col--) {
   
      /* iterate over non-zero rows in this column */
      for (k=jcL[col+1]-1; k>=jcL[col]; k--) { 
    
        /* k is the linear index into the nonzero elements wrt L*/

        /* zero-based row */
        row = irL[k];
    
        /* find corresponding index of S; try deterministically  */
        for (idxS=jcS[col]; idxS<jcS[col+1]; idxS++) {
       
           if (irS[idxS]==row) {          
            break;
          }
        }                
        
        /* fill the rows */        
        rmin = n;
        for (i=jcL[col+1]-1; i>=jcL[col]; i--) {

          if (irL[i]<=col)
            break; /* jump out */
     
          if (irL[i]<rmin)
            rmin = irL[i]; /* minimum row in L */
          
          rowL[irL[i]] = prL[i];          
        }
        
        if (rmin<n) {

         for (i=jcS[row+1]-1; i>=jcS[row]-1; i--) {

          if (irS[i]<=col || irS[i]<rmin)
            break; /* jump out */     
            
          prS[idxS] -= rowL[irS[i]] * prS[i];
         }
         
        }
       
        /* set rows to zero */
        for (i=jcL[col+1]-1; i>=jcL[col]; i--) {

          if (irL[i]<=col)
            break; /* jump out */
     
          rowL[irL[i]] = 0.0;
        }
        
        /* divide by L(c,c) */
        prS[idxS] /= prL[diagL[col]];
              
        if (row==col) {
          
          prS[idxS] += 1.0/(prL[diagL[col]]*prL[diagL[col]]);
  
        }  
        else if (prS[idxS] != 0) {
                      
           /* find corresponding index of S symmetric; try deterministically */
           for (i=jcS[row]; i<jcS[row+1]; i++) {
             
             if (irS[i]==col) {
               
               prS[i] = prS[idxS];
               
               break;
             }
           }
           
         }
         
         
                 
        
      }    
  }  
  
  /* clean up the mess */
  mxFree(diagL);
  mxFree(rowL);
 
}


/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    
  mwSize nzS, n;
  nzS = mxGetM(prhs[1]);

  /* create matrix */
  n = mxGetN(prhs[0]);
  plhs[0] = mxCreateSparse(n,n, nzS, mxREAL);
  
  /* input parameters are L, irS and jcS */
  finv(plhs[0],prhs[0],mxGetPr(prhs[1]),mxGetPr(prhs[2]));
  
  
}
