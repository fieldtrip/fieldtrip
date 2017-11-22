#include <math.h>
#include <matrix.h>
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /*declare variables*/
  const mwSize *dims;
  mwSize  *dimsout;
  mwIndex indx;
  int i, j, k, numdims;
  int numelin;
  mxClassID classid;
  double *input1r, *input1i, *output1r, *output1i;
  double x[3][3], adjx[3][3];
  double xi[3][3], adjxi[3][3];
  double D,Di,Dabs;

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
      for (j=0; j<3; j++)
        for (k=0; k<3; k++)
        {x[k][j]=input1r[i*9+j*3+k];}
          
      adjx[0][0] =  x[1][1]*x[2][2]-x[1][2]*x[2][1];
      adjx[1][0] = -x[1][0]*x[2][2]+x[1][2]*x[2][0];
      adjx[2][0] =  x[1][0]*x[2][1]-x[1][1]*x[2][0];
      adjx[0][1] = -x[0][1]*x[2][2]+x[0][2]*x[2][1];
      adjx[1][1] =  x[0][0]*x[2][2]-x[0][2]*x[2][0];
      adjx[2][1] = -x[0][0]*x[2][1]+x[0][1]*x[2][0];
      adjx[0][2] =  x[0][1]*x[1][2]-x[0][2]*x[1][1];
      adjx[1][2] = -x[0][0]*x[1][2]+x[0][2]*x[1][0];
      adjx[2][2] =  x[0][0]*x[1][1]-x[0][1]*x[1][0];
      
      D = adjx[0][0]*x[0][0]+adjx[1][0]*x[0][1]+adjx[2][0]*x[0][2];
      
      output1r[i*9   ] = adjx[0][0]/D;
      output1r[i*9+1 ] = adjx[1][0]/D;
      output1r[i*9+2 ] = adjx[2][0]/D;
      output1r[i*9+3 ] = adjx[0][1]/D;
      output1r[i*9+4 ] = adjx[1][1]/D;
      output1r[i*9+5 ] = adjx[2][1]/D;
      output1r[i*9+6 ] = adjx[0][2]/D;
      output1r[i*9+7 ] = adjx[1][2]/D;
      output1r[i*9+8 ] = adjx[2][2]/D;
    }
    return;
  }
  else
  {  
    for (i=0; i<numelin/9; i++)
    {
      for (j=0; j<3; j++)
        for (k=0; k<3; k++)
        {
          x[k][j]  = input1r[i*9+j*3+k];
          xi[k][j] = input1i[i*9+j*3+k];
        }
      
      adjx[0][0] =  x[1][1]*x[2][2]-x[1][2]*x[2][1]-xi[1][1]*xi[2][2]+xi[1][2]*xi[2][1];
      adjx[1][0] = -x[1][0]*x[2][2]+x[1][2]*x[2][0]+xi[1][0]*xi[2][2]-xi[1][2]*xi[2][0];
      adjx[2][0] =  x[1][0]*x[2][1]-x[1][1]*x[2][0]-xi[1][0]*xi[2][1]+xi[1][1]*xi[2][0];
      adjxi[0][0] =  x[1][1]*xi[2][2]-x[1][2]*xi[2][1]+xi[1][1]*x[2][2]-xi[1][2]*x[2][1];
      adjxi[1][0] = -x[1][0]*xi[2][2]+x[1][2]*xi[2][0]-xi[1][0]*x[2][2]+xi[1][2]*x[2][0];
      adjxi[2][0] =  x[1][0]*xi[2][1]-x[1][1]*xi[2][0]+xi[1][0]*x[2][1]-xi[1][1]*x[2][0];
      
      D  = adjx[0][0]*x[0][0]+adjx[1][0]*x[0][1]+adjx[2][0]*x[0][2]-adjxi[0][0]*xi[0][0]-adjxi[1][0]*xi[0][1]-adjxi[2][0]*xi[0][2];
      Di = adjx[0][0]*xi[0][0]+adjx[1][0]*xi[0][1]+adjx[2][0]*xi[0][2]+adjxi[0][0]*x[0][0]+adjxi[1][0]*x[0][1]+adjxi[2][0]*x[0][2];
      Dabs = D*D+Di*Di;
      
      adjx[0][1] = -x[0][1]*x[2][2]+x[0][2]*x[2][1]+xi[0][1]*xi[2][2]-xi[0][2]*xi[2][1];
      adjx[1][1] =  x[0][0]*x[2][2]-x[0][2]*x[2][0]-xi[0][0]*xi[2][2]+xi[0][2]*xi[2][0];
      adjx[2][1] = -x[0][0]*x[2][1]+x[0][1]*x[2][0]+xi[0][0]*xi[2][1]-xi[0][1]*xi[2][0];
      adjx[0][2] =  x[0][1]*x[1][2]-x[0][2]*x[1][1]-xi[0][1]*xi[1][2]+xi[0][2]*xi[1][1];
      adjx[1][2] = -x[0][0]*x[1][2]+x[0][2]*x[1][0]+xi[0][0]*xi[1][2]-xi[0][2]*xi[1][0];
      adjx[2][2] =  x[0][0]*x[1][1]-x[0][1]*x[1][0]-xi[0][0]*xi[1][1]+xi[0][1]*xi[1][0];
      
      adjxi[0][1] = -x[0][1]*xi[2][2]+x[0][2]*xi[2][1]-xi[0][1]*x[2][2]+xi[0][2]*x[2][1];
      adjxi[1][1] =  x[0][0]*xi[2][2]-x[0][2]*xi[2][0]+xi[0][0]*x[2][2]-xi[0][2]*x[2][0];
      adjxi[2][1] = -x[0][0]*xi[2][1]+x[0][1]*xi[2][0]-xi[0][0]*x[2][1]+xi[0][1]*x[2][0];
      adjxi[0][2] =  x[0][1]*xi[1][2]-x[0][2]*xi[1][1]+xi[0][1]*x[1][2]-xi[0][2]*x[1][1];
      adjxi[1][2] = -x[0][0]*xi[1][2]+x[0][2]*xi[1][0]-xi[0][0]*x[1][2]+xi[0][2]*x[1][0];
      adjxi[2][2] =  x[0][0]*xi[1][1]-x[0][1]*xi[1][0]+xi[0][0]*x[1][1]-xi[0][1]*x[1][0];
      
      output1r[i*9   ] = (D*adjx[0][0]+Di*adjxi[0][0])/Dabs;
      output1r[i*9+1 ] = (D*adjx[1][0]+Di*adjxi[1][0])/Dabs;
      output1r[i*9+2 ] = (D*adjx[2][0]+Di*adjxi[2][0])/Dabs;
      output1r[i*9+3 ] = (D*adjx[0][1]+Di*adjxi[0][1])/Dabs;
      output1r[i*9+4 ] = (D*adjx[1][1]+Di*adjxi[1][1])/Dabs;
      output1r[i*9+5 ] = (D*adjx[2][1]+Di*adjxi[2][1])/Dabs;
      output1r[i*9+6 ] = (D*adjx[0][2]+Di*adjxi[0][2])/Dabs;
      output1r[i*9+7 ] = (D*adjx[1][2]+Di*adjxi[1][2])/Dabs;
      output1r[i*9+8 ] = (D*adjx[2][2]+Di*adjxi[2][2])/Dabs;
      
      output1i[i*9   ] = (D*adjxi[0][0]-Di*adjx[0][0])/Dabs;
      output1i[i*9+1 ] = (D*adjxi[1][0]-Di*adjx[1][0])/Dabs;
      output1i[i*9+2 ] = (D*adjxi[2][0]-Di*adjx[2][0])/Dabs;
      output1i[i*9+3 ] = (D*adjxi[0][1]-Di*adjx[0][1])/Dabs;
      output1i[i*9+4 ] = (D*adjxi[1][1]-Di*adjx[1][1])/Dabs;
      output1i[i*9+5 ] = (D*adjxi[2][1]-Di*adjx[2][1])/Dabs;
      output1i[i*9+6 ] = (D*adjxi[0][2]-Di*adjx[0][2])/Dabs;
      output1i[i*9+7 ] = (D*adjxi[1][2]-Di*adjx[1][2])/Dabs;
      output1i[i*9+8 ] = (D*adjxi[2][2]-Di*adjx[2][2])/Dabs;
              
    }
    return;
  }
}
