#include <math.h>
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /*declare variables*/
  const mwSize *dims;
  mwIndex indx;
  int i, numdims, indx1, indx2, indx3;
  int numelin;
  mxClassID classid;
  double *input1r, *input1i, *input2r, *input2i, *output1r, *output1i;
  double a[3][3],  b[3][3], c[3][3];
  double ai[3][3], bi[3][3], ci[3][3];
  
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
    for (i=0; i<numelin/9; i++)
    {
      /* real-valued case*/
      for (indx2=0; indx2<3; ++indx2)
        for (indx1=0; indx1<3; ++indx1)
          {
           a[indx1][indx2] = input1r[i*9+indx1+indx2*3];
           b[indx1][indx2] = input2r[i*9+indx1+indx2*3];
           
           c[indx1][indx2] = 0;
           output1r[i*9+indx1+indx2*3] = 0;
          }
      
      for (indx3=0; indx3<3; ++indx3)
        for (indx2=0; indx2<3; ++indx2)
          for (indx1=0; indx1<3; ++indx1)
            {
             c[indx2][indx3] += a[indx2][indx1] * b[indx1][indx3];
            }

      for (indx3=0; indx3<3; ++indx3)
        for (indx2=0; indx2<3; ++indx2)
          for (indx1=0; indx1<3; ++indx1)
            {
             output1r[i*9+indx2+indx3*3] += c[indx2][indx1] * a[indx3][indx1];
            }
    }
    return;
  }
  else if (input1i == NULL && input2i != NULL) 
  {
    for (i=0; i<numelin/9; i++)
    {
      /* first input real-valued case, second input complex*/
      for (indx2=0; indx2<3; ++indx2)
        for (indx1=0; indx1<3; ++indx1)
          {
           a[indx1][indx2]  = input1r[i*9+indx1+indx2*3];
           b[indx1][indx2]  = input2r[i*9+indx1+indx2*3];
           bi[indx1][indx2] = input2i[i*9+indx1+indx2*3];
          
           c[indx1][indx2] = 0;
           ci[indx1][indx2] = 0;
           output1r[i*9+indx1+indx2*3] = 0;
           output1i[i*9+indx1+indx2*3] = 0;
          }
      
      for (indx3=0; indx3<3; ++indx3)
        for (indx2=0; indx2<3; ++indx2)
          for (indx1=0; indx1<3; ++indx1)
            {
             c[indx2][indx3]  += a[indx2][indx1] * b[indx1][indx3];
             ci[indx2][indx3] += a[indx2][indx1] * bi[indx1][indx3];
            }

      for (indx3=0; indx3<3; ++indx3)
        for (indx2=0; indx2<3; ++indx2)
          for (indx1=0; indx1<3; ++indx1)
            {
             output1r[i*9+indx2+indx3*3] += c[indx2][indx1] * a[indx3][indx1];
             output1i[i*9+indx2+indx3*3] += ci[indx2][indx1] * a[indx3][indx1];
            }
    }
    return;  
    
  }
  else if (input1i != NULL && input2i == NULL)
  {
    for (i=0; i<numelin/9; i++)
    {
      /* first input complex-valued, second input real-valued*/
      for (indx2=0; indx2<3; ++indx2)
        for (indx1=0; indx1<3; ++indx1)
          {
           a[indx1][indx2]  = input1r[i*9+indx1+indx2*3];
           b[indx1][indx2]  = input2r[i*9+indx1+indx2*3];
           ai[indx1][indx2] = input1i[i*9+indx1+indx2*3];
          
           c[indx1][indx2] = 0;
           ci[indx1][indx2] = 0;
           output1r[i*9+indx1+indx2*3] = 0;
           output1i[i*9+indx1+indx2*3] = 0;
          }
      
      for (indx3=0; indx3<3; ++indx3)
        for (indx2=0; indx2<3; ++indx2)
          for (indx1=0; indx1<3; ++indx1)
            {
             c[indx2][indx3] += a[indx2][indx1] * b[indx1][indx3];
             ci[indx2][indx3] += ai[indx2][indx1] * b[indx1][indx3];
            }

      for (indx3=0; indx3<3; ++indx3)
        for (indx2=0; indx2<3; ++indx2)
          for (indx1=0; indx1<3; ++indx1)
            {
             output1r[i*9+indx2+indx3*3] += c[indx2][indx1] * a[indx3][indx1] + ci[indx2][indx1] * ai[indx3][indx1];
             output1i[i*9+indx2+indx3*3] += ci[indx2][indx1] * a[indx3][indx1] - c[indx2][indx1] * ai[indx3][indx1];
            }
    }
    return;  
  }  
  else
  {  
    for (i=0; i<numelin/9; i++)
    {
      /* both inputs complex-valued*/
      for (indx2=0; indx2<3; ++indx2)
        for (indx1=0; indx1<3; ++indx1)
          {
           a[indx1][indx2]  = input1r[i*9+indx1+indx2*3];
           b[indx1][indx2]  = input2r[i*9+indx1+indx2*3];
           ai[indx1][indx2] = input1i[i*9+indx1+indx2*3];
           bi[indx1][indx2] = input2i[i*9+indx1+indx2*3];
          
           c[indx1][indx2] = 0;
           ci[indx1][indx2] = 0;
           output1r[i*9+indx1+indx2*3] = 0;
           output1i[i*9+indx1+indx2*3] = 0;
          }
      
      for (indx3=0; indx3<3; ++indx3)
        for (indx2=0; indx2<3; ++indx2)
          for (indx1=0; indx1<3; ++indx1)
            {
             c[indx2][indx3] += a[indx2][indx1] * b[indx1][indx3] - ai[indx2][indx1] * bi[indx1][indx3];
             ci[indx2][indx3] += ai[indx2][indx1] * b[indx1][indx3] + a[indx2][indx1] * bi[indx1][indx3];
            }

      for (indx3=0; indx3<3; ++indx3)
        for (indx2=0; indx2<3; ++indx2)
          for (indx1=0; indx1<3; ++indx1)
            {
             output1r[i*9+indx2+indx3*3] += c[indx2][indx1] * a[indx3][indx1] + ci[indx2][indx1] * ai[indx3][indx1];
             output1i[i*9+indx2+indx3*3] += ci[indx2][indx1] * a[indx3][indx1] - c[indx2][indx1] * ai[indx3][indx1];
            }
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

