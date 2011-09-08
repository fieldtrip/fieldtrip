#include <math.h>
#include <matrix.h>
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /*declare variables*/
  const mwSize *dims;
  mwIndex indx;
  int i, numdims;
  int numelin;
  mxClassID classid;
  double *input1r, *input1i, *input2r, *input2i, *output1r, *output1i;
  double a,  b,  c,  d,  e,  f,  g,  h;
  double ai, bi, ci, di, ei, fi, gi, hi;
  
  /*figure out the classid*/
  classid = mxGetClassID(prhs[0]);
     
  /*check inputs*/
  if (nrhs!=2)
    mexErrMsgTxt("Wrong number of input arguments");
  
  /*associate inputs*/
  input1r = mxGetData(prhs[0]);
  input1i = mxGetImagData(prhs[0]);
  
  input2r = mxGetData(prhs[1]);
  input2i = mxGetImagData(prhs[1]);
  
  /*figure out dimension info and number of elements*/
  dims    = mxGetDimensions(prhs[0]);
  numdims = mxGetNumberOfDimensions(prhs[0]);
  numelin = mxGetNumberOfElements(prhs[0]);
    
  /*associate output*/
  if (input1i == NULL && input2i == NULL)
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
  if (input1i == NULL && input2i == NULL)
  {  
    for (i=0; i<numelin/4; i++)
    {
      a = input1r[i*4  ];
      b = input1r[i*4+1];
      c = input1r[i*4+2];
      d = input1r[i*4+3];
    
      e = input2r[i*4  ];
      f = input2r[i*4+1];
      g = input2r[i*4+2];
      h = input2r[i*4+3];
    
      output1r[i*4  ] = a*e + c*f;
      output1r[i*4+1] = b*e + d*f;
      output1r[i*4+2] = a*g + c*h;
      output1r[i*4+3] = b*g + d*h;
    }
    return;
  }
  else if (input1i == NULL)
  {  
    for (i=0; i<numelin/4; i++)
    {
      /*matrix 1*/
      a  = input1r[i*4  ]; b  = input1r[i*4+1]; c  = input1r[i*4+2]; d  = input1r[i*4+3];
      
      /*matrix 2*/
      e  = input2r[i*4  ]; f  = input2r[i*4+1]; g  = input2r[i*4+2]; h  = input2r[i*4+3];
      ei = input2i[i*4  ]; fi = input2i[i*4+1]; gi = input2i[i*4+2]; hi = input2i[i*4+3];
    
      /*fill in the real part of the output matrix*/
      output1r[i*4  ] = a*e + c*f;
      output1r[i*4+1] = b*e + d*f;
      output1r[i*4+2] = a*g + c*h;
      output1r[i*4+3] = b*g + d*h;
      
      /*fill in the imaginary part of the output matrix*/
      output1i[i*4  ] = a*ei + c*fi;
      output1i[i*4+1] = b*ei + d*fi;
      output1i[i*4+2] = a*gi + c*hi;
      output1i[i*4+3] = b*gi + d*hi;
    }
    return;
  }
  else if (input2i == NULL)
  {  
    for (i=0; i<numelin/4; i++)
    {
      /*matrix 1*/
      a  = input1r[i*4  ]; b  = input1r[i*4+1]; c  = input1r[i*4+2]; d  = input1r[i*4+3];
      ai = input1i[i*4  ]; bi = input1i[i*4+1]; ci = input1i[i*4+2]; di = input1i[i*4+3];
      
      /*matrix 2*/
      e  = input2r[i*4  ]; f  = input2r[i*4+1]; g  = input2r[i*4+2]; h  = input2r[i*4+3];
      
      /*fill in the real part of the output matrix*/
      output1r[i*4  ] = a*e + c*f;
      output1r[i*4+1] = b*e + d*f;
      output1r[i*4+2] = a*g + c*h;
      output1r[i*4+3] = b*g + d*h;
      
      /*fill in the imaginary part of the output matrix*/
      output1i[i*4  ] = ai*e + ci*f;
      output1i[i*4+1] = bi*e + di*f;
      output1i[i*4+2] = ai*g + ci*h;
      output1i[i*4+3] = bi*g + di*h;
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
      
      /*matrix 2*/
      e  = input2r[i*4  ]; f  = input2r[i*4+1]; g  = input2r[i*4+2]; h  = input2r[i*4+3];
      ei = input2i[i*4  ]; fi = input2i[i*4+1]; gi = input2i[i*4+2]; hi = input2i[i*4+3];
    
      /*fill in the real part of the output matrix*/
      output1r[i*4  ] = (a*e-ai*ei) + (c*f-ci*fi);
      output1r[i*4+1] = (b*e-bi*ei) + (d*f-di*fi);
      output1r[i*4+2] = (a*g-ai*gi) + (c*h-ci*hi);
      output1r[i*4+3] = (b*g-bi*gi) + (d*h-di*hi);
      
      /*fill in the imaginary part of the output matrix*/
      output1i[i*4  ] = (ai*e+a*ei) + (ci*f+c*fi);
      output1i[i*4+1] = (bi*e+b*ei) + (di*f+d*fi);
      output1i[i*4+2] = (ai*g+a*gi) + (ci*h+c*hi);
      output1i[i*4+3] = (bi*g+b*gi) + (di*h+d*hi);
    }
    return;
  }
}
