#include "mex.h"

/*
 * mex implementation of Rue's method.
 *
 * compile with mex -largeArrayDims fastinvre64.c
 *
 * Marcel van Gerven, 2009
 *
 * 01-06-2012: Memory leak fixed by Tomi Peltola
 *
 */

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
     
  mwSize  n64, nzL64, nzS64;
  mwIndex *irL64, *jcL64, *irS64, *jcS64, *diagL64;  
  mwIndex row64, idxS64, i64, k64, col64, rmin64;   
  double  *prL64, *prS64, *rowL64;

  const mxArray *L, *CS;
  mxArray *S; 
  
  double *cirS, *cjcS;
  
  n64 = mxGetN(prhs[0]);
  nzS64 = mxGetM(prhs[1]);

  L = prhs[0];

  /* create matrix */
  S = mxCreateSparse(n64,n64, nzS64, mxREAL);
  
  /* initialize the covariance matrix */
  prS64 = mxGetPr(S);
    
  cirS = mxGetPr(prhs[1]);
  irS64 = mxGetIr(S);
  for (k64=0;k64<nzS64;k64++) {
    irS64[k64] = (mwIndex)cirS[k64];
  }    

   
  cjcS = mxGetPr(prhs[2]);
  jcS64 = mxGetJc(S);
  for (k64=0;k64<=n64;k64++) {
    jcS64[k64] = (mwIndex)(cjcS[k64]);
  }    

  /* number of rows and columns in L */
  n64 = mxGetN(L);
    
  /* Get the starting positions of all data arrays in L */ 
  prL64 = mxGetPr(L);
  irL64 = mxGetIr(L);
  jcL64 = mxGetJc(L);
 
  /* Get the starting positions of all data arrays in S */ 
  prS64 = mxGetPr(S);
  irS64 = mxGetIr(S);
  jcS64 = mxGetJc(S);    
  
  /* number of non-zero elements in S */
  mxSetNzmax(S,nzS64);
 
  /* set entries to zero */
  for (i64=0; i64<nzS64; i64++) {
   prS64[i64] = 0; 
 }
 
  /* find indices belonging to diagonal elements */
  
  diagL64 = mxCalloc(n64, sizeof(mwIndex));
  
  for (col64=0; col64<n64; col64++) {
        
    for (k64=jcL64[col64]; k64<jcL64[col64+1]; k64++) {

      if (irL64[k64]==col64) {
        diagL64[col64] = k64; 
        break;
      }
      
    }       
  }
  
  /* iterate over all elements using sparse column and row loops */

  /* create space for the L and S rows */
  rowL64 = mxCalloc(n64, sizeof(double));   
     
  for (col64=n64-1; (int)(col64)>=0; col64--) {
 
      for (k64=jcL64[col64+1]-1; (int)(k64-jcL64[col64])>=0; k64--) { 
    
        row64 = irL64[k64];
        
        for (idxS64=jcS64[col64]; (int)(idxS64-jcS64[col64+1])<0; idxS64++) {
         
          if ((int)(irS64[idxS64]-row64)==0) {          
            break;
          }
        }         

        rmin64 = n64;
        for (i64=jcL64[col64+1]-1; (int)(i64-jcL64[col64])>=0; i64--) {

          
          if ((int)(irL64[i64]-col64)<=0)
            break;
     
          if ((int)(irL64[i64]-rmin64)<0)
            rmin64 = irL64[i64];
          
          rowL64[irL64[i64]] = prL64[i64];          
          
        }                

        
        if ((int)(rmin64-n64)<0) {

         for (i64=jcS64[row64+1]-1; (int)(i64-jcS64[row64]-1)>=0; i64--) {

          if ((int)(irS64[i64]-col64)<=0 || (int)(irS64[i64]-rmin64)<0)
            break; 
            
          prS64[idxS64] -= rowL64[irS64[i64]] * prS64[i64];
         }
         
        }
               
        for (i64=jcL64[col64+1]-1; (int)(i64-jcL64[col64])>=0; i64--) {

          if ((int)(irL64[i64]-col64)<=0)
            break; 
     
          rowL64[irL64[i64]] = 0.0;
        }
        
        prS64[idxS64] /= prL64[diagL64[col64]];
              
        if ((int)(row64-col64)==0) {
          
          prS64[idxS64] += 1.0/(prL64[diagL64[col64]]*prL64[diagL64[col64]]);
  
        }  
        else  {
                      
           for (i64=jcS64[row64]; (int)(i64-jcS64[row64+1])<0; i64++) {
             
             if ((int)(irS64[i64]-col64)==0) {
               
               prS64[i64] = prS64[idxS64];
               
               break;
             }
           }
           
         }
       
      }  
  }     

  /* clean up the mess */
    
  mxFree(diagL64);
  mxFree(rowL64);    
  
  plhs[0] = S;  
  
  
}
