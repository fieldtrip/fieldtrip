#include <math.h>
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /*declare variables*/
  const mwSize *dims;
  mwSize  *dimsout;
  mwIndex indx;
  int i, numdims;
  int numelin;
  mxClassID classid;
  double *input1r, *input1i, *output1r, *output1i;
  double a,  b,  c,  d,  e,  f,  g,  h,  j;
  double ai, bi, ci, di, ei, fi, gi, hi ,ji;
  
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
  
  dimsout    = mxMalloc(numdims * sizeof(mwSize));
  for (i=0; i<numdims; i++)
  {
    dimsout[i] = dims[i];
  }
  dimsout[0] = 1;
  dimsout[1] = 1;
  
  /*associate output*/
  if (input1i == NULL)
  {
    plhs[0]  = mxCreateNumericArray(numdims, dimsout, classid, mxREAL);
    output1r = mxGetData(plhs[0]);
  }
  else
  {
    plhs[0]  = mxCreateNumericArray(numdims, dimsout, classid, mxCOMPLEX);
    output1r = mxGetData(plhs[0]);
    output1i = mxGetImagData(plhs[0]);
  }
  
  /* do the computation*/
  if (input1i == NULL)
  {  
    for (i=0; i<numelin/9; i++)
    {
      a = input1r[i*9  ];
      b = input1r[i*9+1];
      c = input1r[i*9+2];
      d = input1r[i*9+3];
      e = input1r[i*9+4];
      f = input1r[i*9+5];
      g = input1r[i*9+6];
      h = input1r[i*9+7];
      j = input1r[i*9+8];
      
      output1r[i] = a*e*j - a*h*f - d*b*j + d*h*c +g*b*f - g*e*c;
    }
    return;
  }
  else
  {  
    for (i=0; i<numelin/9; i++)
    {
      /*matrix 1*/
      a  = input1r[i*9  ]; b  = input1r[i*9+1]; c  = input1r[i*9+2]; d  = input1r[i*9+3];
      e  = input1r[i*9+4]; f  = input1r[i*9+5]; g  = input1r[i*9+6]; h  = input1r[i*9+7];
      j  = input1r[i*9+8];
      
      ai = input1i[i*9  ]; bi = input1i[i*9+1]; ci = input1i[i*9+2]; di = input1i[i*9+3];
      ei = input1i[i*9+4]; fi = input1i[i*9+5]; gi = input1i[i*9+6]; hi = input1i[i*9+7];
      ji = input1i[i*9+8];
      
      /*fill in the real part of the output matrix*/
      output1r[i] = (a*e*j-ai*ei*j-a*ei*ji-ai*e*ji) - (a*h*f-ai*hi*f-a*hi*fi-ai*h*fi) - (d*b*j-di*bi*j-d*bi*ji-di*b*ji) + (d*h*c-di*hi*c-d*hi*ci-di*h*ci) + (g*b*f-gi*bi*f-g*bi*fi-gi*b*fi) - (g*e*c-gi*ei*c-g*ei*ci-gi*e*ci);
      
      /*fill in the imaginary part of the output matrix*/
      output1i[i] = (a*ei*j+ai*e*j+a*e*ji-ai*ei*ji) - (a*hi*f+ai*h*f+a*h*fi-ai*hi*fi) - (d*bi*j+di*b*j+d*b*ji-di*bi*ji) + (d*hi*c+di*h*c+d*h*ci-di*hi*ci) + (g*bi*f+gi*b*f+g*b*fi-gi*bi*fi) - (g*ei*c+gi*e*c+g*e*ci-gi*ei*ci);
    }
    return;
  }
}
