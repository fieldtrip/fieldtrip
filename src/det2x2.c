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
  double a,  b,  c,  d;
  double ai, bi, ci, di;
  
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
    for (i=0; i<numelin/4; i++)
    {
      a = input1r[i*4  ];
      b = input1r[i*4+1];
      c = input1r[i*4+2];
      d = input1r[i*4+3];
      output1r[i] = a*d - b*c;
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
      
      /*fill in the real part of the output matrix*/
      output1r[i] = (a*d-ai*di) - (b*c-bi*ci);
      
      /*fill in the imaginary part of the output matrix*/
      output1i[i] = (ai*d+a*di) - (bi*c+b*ci);
    }
    return;
  }
}
