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
  double a,  b,  c,  d,  e,  f,  g,  h, absa, absb, absc, absd, offr, offi;
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
    
      output1r[i*4  ] = e*a*a + 2*f*a*c       + h*c*c;
      output1r[i*4+1] = e*a*b + f*a*d + f*b*c + h*c*d;
      output1r[i*4+2] = e*a*b + f*b*c + f*a*d + h*c*d;
      output1r[i*4+3] = e*b*b + 2*f*b*d       + h*d*d;
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
      
      /*compute some quantities only once*/
      absa = a*a;
      absb = b*b;
      absc = c*c;
      absd = (d*d+di*di);
      offr = e*a*b + f*a*d + f*b*c + h*c*d;
      offi = -ei*a*b - fi*a*d + fi*b*c - hi*c*d;
      
      /*fill in the real part of the output matrix*/
      output1r[i*4  ] = e*absa + 2*f*a*c + h*absc;
      output1r[i*4+1] = offr;
      output1r[i*4+2] = offr;
      output1r[i*4+3] = e*absb + 2*f*b*d + h*absd;
    
      /*fill in the imaginary part of the output matrix*/
      output1i[i*4  ] = ei*absa + hi*absc;
      output1i[i*4+1] = -offi;
      output1i[i*4+2] =  offi;
      output1i[i*4+3] = ei*absb + hi*absd;
    
    }
  }
  /*else if (input2i == NULL)*/
    
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
    
      /*compute some quantities only once*/
      absa = (a*a+ai*ai);
      absb = (b*b+bi*bi);
      absc = (c*c+ci*ci);
      absd = (d*d+di*di);
      offr = (e*a*b+e*ai*bi-ei*a*bi+ei*ai*b) + (f*a*d+f*ai*di-fi*a*di+fi*ai*d) + (f*b*c+f*bi*ci-fi*b*ci+fi*bi*c) + (h*c*d+h*ci*di-hi*c*di+hi*ci*d);
      offi = (-ei*ai*bi-e*a*bi+e*ai*b-ei*a*b) + (-fi*ai*di-f*a*di+f*ai*d-fi*a*d) + (fi*bi*ci+f*b*ci-f*bi*c+fi*b*c) + (-hi*ci*di-h*c*di+h*ci*d-hi*c*d);
      
      /*fill in the real part of the output matrix*/
      output1r[i*4  ] = e*absa + 2*(f*a*c+f*ai*ci-fi*a*ci+fi*ai*c) + h*absc;
      output1r[i*4+1] = offr;
      output1r[i*4+2] = offr;
      output1r[i*4+3] = e*absb + 2*(f*b*d+f*bi*di-fi*b*di+fi*bi*d) + h*absd;
    
      /*fill in the imaginary part of the output matrix*/
      output1i[i*4  ] = ei*absa + hi*absc;
      output1i[i*4+1] = -offi;
      output1i[i*4+2] =  offi;
      output1i[i*4+3] = ei*absb + hi*absd;
    }
    return;
  }
}

/*in the following c and b are swapped with respect to the above
%a b     e f'   a' c'
%c d     f h    b' d'
%
%a*e+b*f  a*f'+b*h  a' c'
%c*e+d*f  c*f'+d*h  b' d'
%
%a*e*a'+b*f*a'+a*f'*b'+b*h*b' a*e*c'+b*f*c'+a*f'*d'+b*h*d'
%c*e*a'+d*f*a'+c*f'*b'+d*h*b' c*e*c'+d*f*c'+c*f'*d'+d*h*d'
%
%e*abs(a)^2 + f*(a'*b) + f'*(a*b') + h*abs(b)^2    e*a*c'    + f*b*c'   + f'*a*d'   + h*b*d'
%e*a'*c    + f*a'*d   + f'*b'*c   + h*b'*d       e*abs(c)^2 + f*(c'*d) + f'*(c*d') + h*abs(d)^2*/

