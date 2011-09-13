#include <mex.h>

#include "FIL_Vista.h"
#include <stdio.h>
#include <string>

// main function
void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
 	// check number of arguments
    int nargin = nrhs;
    if(nargin != 3)
    {
        mexErrMsgTxt("Not enough arguments. Usage: <dim>, <seg>, <filename>\n");
    }
    else if(!mxIsDouble(prhs[0]))
    {
        mexErrMsgTxt("Wrong type of argument. <dim> must be of type 'double'.\n");
    }
    else if(!mxIsDouble(prhs[1]))
    {
        mexErrMsgTxt("Wrong type of argument. <seg> must be of type 'double'.\n");
    }
    else if(!mxIsChar(prhs[2]))
    {
        mexErrMsgTxt("Wrong type of argument. <filename> must be of type 'char'.\n");
    }
    else if(!((mxGetM(prhs[0])==1)&&(mxGetN(prhs[0])==3)))
    {
        mexErrMsgTxt("Wrong dimension. <dim> is supposed to be a 1x3 array.\n");
    }
    
    double *dim, *seg;
    dim = mxGetPr(prhs[0]);
    seg = mxGetPr(prhs[1]);
    
    /*int nlhs1, nrhs1;
    mxArray *plhs1[1];
    const mxArray *prhs1[1];
    nlhs1 = 1;
    nrhs1 = 1;
    prhs1[0] = prhs[1];
    mexCallMATLAB(nlhs1,plhs1,nrhs1,prhs1,"size");
    
    double *dimMat;
    dimMat = mxGetPr(prhs1);

    if (!((dim[0]==dimMat[0])&&(dim[1]==dimMat[1])&&(dim[2]==dimMat[2])))
    {
        mexErrMsgTxt("Wrong dimension. <seg> is supposed to be of dimension <dim>.\n");
    }*/
    
	// get filename from arguments
   std::string filename(mxArrayToString(prhs[2]));
   FILE *outf;
   outf = fopen (filename.c_str(),"w");

   VImage dst;
   dst = VCreateImage(dim[0],dim[1],dim[2], VUByteRepn);
   
   char voxel_string[100];
   float voxelx = 1.0;
   float voxely = 1.0;
   float voxelz = 1.0;
   sprintf(voxel_string, "%f %f %f", voxelx, voxely, voxelz);
   VSetAttr(VImageAttrList(dst), "voxel", NULL, VStringRepn, voxel_string);
   
	 for(int band=0; band<dim[2]; band++)
	   for(int row=0; row<dim[1]; row++) 
	     for(int col=0; col<dim[0]; col++)
         {
    		 VPixel(dst, band, row, col, VUByte)=seg[col+((row+band*(int)dim[1])*(int)dim[0])];
	     }
    
    /*write out the resulting image*/
    VWriteImages(outf, NULL, 1, &dst) ; 
fclose(outf);
}