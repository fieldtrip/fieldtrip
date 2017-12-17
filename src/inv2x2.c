#include <math.h>
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /*declare variables*/
  const mwSize *dims;
  mwIndex indx;
  int i, numdims;
  int numelin;
  mxClassID classid;
  double *input1r, *input1i, *output1r, *output1i;
  double a,  b,  c,  d,  denom;
  double ai, bi, ci, di, denomi, denomabs;
  
  /*figure out the classid*/
  classid = mxGetClassID(prhs[0]);
     
  /*check inputs*/
  if (nrhs>1)
    mexErrMsgTxt("Too many input arguments");
  
  /*associate inputs*/
  input1r = mxGetData(prhs[0]);
  input1i = mxGetImagData(prhs[0]);
  
  /*figure out dimension info and number of elements*/
  dims    = mxGetDimensions(prhs[0]);
  numdims = mxGetNumberOfDimensions(prhs[0]);
  numelin = mxGetNumberOfElements(prhs[0]);
  /*associate output*/
  if (input1i == NULL)
  {
    plhs[0]  = mxCreateNumericArray(numdims, dims, classid, mxREAL);
    output1r = mxGetData(plhs[0]);
  }
  else
  {
    plhs[0]  = mxCreateNumericArray(numdims, dims, classid, mxCOMPLEX);
    output1r = mxGetData(plhs[0]);
    output1i = mxGetImagData(plhs[0]);
  }
  
  /* do the computation*/
  if (input1i == NULL)
  {  
    for (i=0; i<numelin/4; i++)
    {
      a = input1r[i*4  ];
      b = input1r[i*4+1];
      c = input1r[i*4+2];
      d = input1r[i*4+3];
      denom = a*d - b*c;
    
      output1r[i*4  ] =  d/denom;
      output1r[i*4+1] = -b/denom;
      output1r[i*4+2] = -c/denom;
      output1r[i*4+3] =  a/denom;
    }
    return;
  }
  else
  {  
    for (i=0; i<numelin/4; i++)
    {
      /*matrix 1*/
      a  = input1r[i*4  ]; b  = input1r[i*4+1]; c  = input1r[i*4+2]; d  = input1r[i*4+3];
      ai = input1i[i*4  ]; bi = input1i[i*4+1]; ci = input1i[i*4+2]; di = input1i[i*4+3];
      
      /*get the determinant*/
      denom    = (a*d-ai*di) - (b*c-bi*ci);
      denomi   = (ai*d+a*di) - (bi*c+b*ci);
      denomabs = denom*denom + denomi*denomi;
      
      /*fill in the real part of the output matrix*/
      output1r[i*4  ] =  (d*denom+di*denomi)/denomabs;
      output1r[i*4+1] = -(b*denom+bi*denomi)/denomabs;
      output1r[i*4+2] = -(c*denom+ci*denomi)/denomabs;
      output1r[i*4+3] =  (a*denom+ai*denomi)/denomabs;
      
      /*fill in the imaginary part of the output matrix*/
      output1i[i*4  ] = -(d*denomi-di*denom)/denomabs;
      output1i[i*4+1] =  (b*denomi-bi*denom)/denomabs;
      output1i[i*4+2] =  (c*denomi-ci*denom)/denomabs;
      output1i[i*4+3] = -(a*denomi-ai*denom)/denomabs;
    }
    return;
  }
  
  /* do the computation*/
  for (i=0; i<numelin/4; i++)
  {
  }
  return;
}
