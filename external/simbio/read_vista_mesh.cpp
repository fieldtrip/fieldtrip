#include <mex.h>

#include "vistaprimitive.h"

#include <iostream>

// main function
void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// check number of arguments
    int nargin = nrhs;
    if(nargin != 1)
    {
        mexErrMsgTxt("Not enough arguments. Usage: <filename>\n");
    }
    else if(!mxIsChar(prhs[0]))
    {
	mexErrMsgTxt("Wrong type of argument. <filename> must be of type 'char'.");
    }

	// get filename from arguments
	std::string filename(mxArrayToString(prhs[0]));

	mexPrintf("Trying to read FE mesh from file %s.\n", filename.c_str());

	// create new cl_ReadVistaFEM object
	bool bError = false;
	cl_ReadVistaFEM *pclFEMesh = NULL;
	if(!bError)
	{
		try {
			pclFEMesh = new cl_ReadVistaFEM; }
		catch (...) {
			pclFEMesh = NULL;
			bError = true;
		}
	}

	// read FE mesh
	if(!bError)
	{
		if(pclFEMesh->_ReadMesh(filename) != 0)
			bError = true;
		//else
		//	mexPrintf("Successfully read FE mesh.\n");
	}

	// query mesh
	std::vector<std::valarray<double> > vvdNodes;
	std::vector<std::valarray<int> > vviElements;
	std::valarray<int> viLabels;
	if(!bError)
	{
		if(pclFEMesh->_GetMesh(vvdNodes, vviElements, viLabels) != 0)
		{
			std::cerr << "Error while trying to query mesh from cl_ReadVistaFEM object." << std::endl;
			bError = true;
		}
	}

	// safely delete cl_ReadVistaFEM object
	if(pclFEMesh)
		delete pclFEMesh;
	pclFEMesh = NULL;

	// cpy mesh to Matlab arrays
	// allox Matlab matrices and vectors
	plhs[0] = mxCreateDoubleMatrix(vvdNodes.size(),3,mxREAL);
	double *matNodes = mxGetPr(plhs[0]);
	plhs[1] = mxCreateDoubleMatrix(vviElements.size(),vviElements[0].size(),mxREAL);
	double *matElements = mxGetPr(plhs[1]);
	plhs[2] = mxCreateDoubleMatrix(viLabels.size(),1,mxREAL);
	double *matLabels = mxGetPr(plhs[2]);

	// cpy nodes
	for(int i=0; i<vvdNodes.size(); ++i)
		for(int k=0; k<3; ++k)
			matNodes[k * vvdNodes.size() + i] = vvdNodes[i][k];

	// cpy elements
	for(int i=0; i<vviElements.size(); ++i)
		for(int k=0; k<vviElements[i].size(); ++k)
			matElements[k*vviElements.size() + i] = vviElements[i][k] + 1; // convert to matlab style indices

    // cpy labels
	for(int i=0; i<viLabels.size(); ++i)
		matLabels[i] = viLabels[i];

}
