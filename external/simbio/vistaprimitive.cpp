#include "vistaprimitive.h"

#include <iostream>

// #include <VGraph.h>
// #include <VImage.h>
#include "FIL_Vista.h"

cl_ReadVistaFEM::cl_ReadVistaFEM(void)
{
	// initialize variables to default values
	_InitVariables();

	return;
}

cl_ReadVistaFEM::~cl_ReadVistaFEM()
{
	// clean up this instance
	_CleanUp();

	return;
}

int cl_ReadVistaFEM::_ReadMesh(std::string sFilename)
{
	int iReturn = 0;
	
	// try to open given file
	FILE *fpInFile = NULL;
	if(iReturn == 0)
	{
		fpInFile = fopen(sFilename.c_str(), "rb");
		if(!fpInFile)
		{
			iReturn = 1;
			std::cerr << "Could not open file " << sFilename << " for reading." << std::endl;
		}
	}

	printf("Read VISTA file.\n");
	// NOTE read vista file using vista lib taking GRI_GridGeneratorVista_c as an example
	// get VISTA attribute list
	VAttrList valList = NULL;
	if(iReturn == 0)
	{
		valList = VReadFile(fpInFile, 0);
		if(valList == NULL)
		{
			iReturn = 1;
			std::cerr << "Could not read VISTA attribute list from file " << sFilename << "." << std::endl;
		}
	}

	// read vertices as first and primitives as second graph element
	VAttrListPosn valpCurrAttribute;
	VGraph vgVertices = 0;
	VGraph vgPrimitives = 0;
	if(iReturn == 0)
	{
		for(VFirstAttr(valList, &valpCurrAttribute); VAttrExists(&valpCurrAttribute);
				VNextAttr(&valpCurrAttribute))
		{
			// check if current attribute represents a graph
			if(VGetAttrRepn(&valpCurrAttribute) != VGraphRepn)
				continue;

			// vertices are read first, so first check if vertices have been
			// read already
			if(vgVertices == 0)
				VGetAttrValue(&valpCurrAttribute, 0, VGraphRepn, &vgVertices);
			else if(vgPrimitives == 0) // read primitives in second place
				VGetAttrValue(&valpCurrAttribute, 0, VGraphRepn, &vgPrimitives);
		}
	}
	if(iReturn == 0)
	{
		// all data has been read at this point, close file
		fclose(fpInFile);
		fpInFile = NULL;
	}

	// get number of partition using lastUsed method, abort if more than 1 partition given
	if(iReturn == 0)
	{
		VLong vlNrPartitions = 0;
		vlNrPartitions = VGetAttr(
				VGraphAttrList(vgVertices), "partition(s)", 0, VLongRepn, &vlNrPartitions);
	}

	// get number of vertices and elements using graph's lastUsed method
	long lNrElements = 0, lNrNodes = 0;
	if(iReturn == 0)
	{
		VLong vlNrVertices, vlNrPrimitives;
		vlNrPrimitives = vgPrimitives->lastUsed;
		vlNrVertices = vgVertices->lastUsed;
		lNrNodes = vlNrVertices;
		lNrElements = vlNrPrimitives;
	}

	// copy node positions from graph structure
	VNode vnCurrNode;
	if(iReturn == 0)
	{
		try {
			_vvdNodes.resize(lNrNodes,
					std::valarray<double> (3)); }
		catch (...) {
			iReturn = 1;
			std::cerr << "Could not allocate enough memory to store " << lNrNodes <<
				" nodes." << std::endl;
		}
	}
	if(iReturn == 0)
	{
		for(int i=1; i<=lNrNodes; ++i)
		{
			vnCurrNode = NULL;
			vnCurrNode = VGraphGetNode(vgVertices, i);
			if(vnCurrNode == NULL)
			{
				iReturn = 1;
				std::cerr << "Could not read vertex no. " << i << "." << std::endl;
			}
			
			if(iReturn != 0)
				break;

			// copy node
			VFloat *pvfNode = reinterpret_cast<VFloat *>(vnCurrNode->data);
			for(int k=0; k<3; ++k)
				_vvdNodes[i-1][k] = static_cast<double>(pvfNode[k+1]);
		}
	}

	// copy elements from graph structure
	// resize element matrix
	if(iReturn == 0)
	{
		try {
			_vviElements.resize(lNrElements,
					std::valarray<int>(0)); }
		catch (...) {
			iReturn = 1;
			std::cerr << "Could not allocate enough memory to store " << lNrElements <<
				" elements." << std::endl;
		}
	}
	if(iReturn == 0)
	{
		for(int i=1; i<=lNrElements; ++i)
		{
			VPrimitive *vpCurrElement = reinterpret_cast<VPrimitive *>(VGraphGetNode(vgPrimitives, i));
			if(vpCurrElement == NULL)
			{
				iReturn = -1;
			}

			if(iReturn != 0)
				break;

			// cpy element
			_vviElements[i-1].resize(vpCurrElement->vcnt);
			for(int k=0; k<vpCurrElement->vcnt; ++k)
				_vviElements[i-1][k] = vpCurrElement->id[k] - 1; // convert to c-style index starting at 0
		}
	}

	// check if graph structure for primitives contains attribute component_interp or matprops ???
	VImage viLabels;
	// create VImage object to store labels for elements
	if(iReturn == 0)
	{
		viLabels = VCreateImage(1, 1, lNrElements, VUByteRepn);
		if(viLabels == NULL)
		{
			std::cerr << "Could not create VImage object to store element labels." << std::endl;
			iReturn = 1;
		}
	}
	// extract element labels
	if(iReturn == 0)
	{
		if(VGetAttr(VGraphAttrList(vgPrimitives), "matprops", 0, VImageRepn, &viLabels) != VAttrFound)
		{
			std::cerr << "No material propertiers (image inside of primitives graph with identifier \"matprops\") found." << std::endl;
			iReturn = 1;
		}
	}

	// copy material properties ( i.e. labels)
	// resize valarray to store labels
	if(iReturn == 0)
	{
		try {
			_viLabels.resize(lNrElements); }
		catch (...) {
			iReturn = 1;
			std::cerr << "Could not allocate enough memory to store " << lNrElements <<
				" labels." << std::endl;
		}
	}
	if(iReturn == 0)
	{
		for(int i=0; i<lNrElements; ++i)
		{
			_viLabels[i] = static_cast<int>(VGetPixel(viLabels, 0, 0, i));
		}
	}

	// finalize method
	if(iReturn == 0) // success
	{
		_iNrVertices = lNrNodes;
		_iNrElements = lNrElements;
		_bMeshPresent = true;

		std::cout << "Successfully read mesh with " << lNrNodes << " nodes and " <<
			lNrElements << " elements." << std::endl;
	}
	else // error
	{
		_iNrElements = 0;
		_iNrVertices = 0;
		_vviElements.resize(0);
		_vvdNodes.resize(0);
		_viLabels.resize(0);
	}
	
	return(iReturn);
}

int cl_ReadVistaFEM::_GetMesh(std::vector<std::valarray<double> > &rvvdNodes,
		std::vector<std::valarray<int> > &rvviElements,
		std::valarray<int> &rviLabels)
{
	int iReturn = 0;

	// check if mesh is present, which could be returned
	if(!_bMeshPresent)
		return(1);

	// cpy nodes
	if(iReturn == 0)
	{
		try {
			rvvdNodes.resize(_vvdNodes.size(),
					std::valarray<double>(3)); }
		catch (...) {
			iReturn = 1;
			std::cerr << "Could not allocate memory for matrix in which nodes are returned." << std::endl;
		}
	}
	if(iReturn == 0)
	{
		// cpy nodes
		rvvdNodes.assign(_vvdNodes.begin(), _vvdNodes.end());
	}

	// cpy elements
	// resize element return matrix
	if(iReturn == 0)
	{
		try {
			rvviElements.resize(_vviElements.size()); }
		catch (...) {
			iReturn = 1;
			std::cerr << "Could not allocate memory for matrix in which elements are returned." << std::endl;
		}
	}
	if(iReturn == 0)
	{
		// cpy elements
		for(int i=0; i<_vviElements.size(); ++i)
		{
			rvviElements[i].resize(_vviElements[i].size());
			rvviElements[i] = _vviElements[i];
		}
	}

	// cpy labels
	if(iReturn == 0)
	{
		try {
			rviLabels.resize(_viLabels.size()); }
		catch (...) {
			iReturn = 1;
			std::cerr << "Could not allocate memory for vector in which labels are returned." << std::endl;
		}
	}
	if(iReturn == 0)
	{
		rviLabels = _viLabels;
	}

	// delete return matrices if copying of mesh failed
	if(iReturn != 0)
	{
		rvviElements.resize(0);
		rvvdNodes.resize(0);
		rviLabels.resize(0);
	}

	return(iReturn);
}

void cl_ReadVistaFEM::_InitVariables(void)
{
	// set variables to default values
	_vviElements.resize(0);
	_vvdNodes.resize(0);
	_viLabels.resize(0);

	_iNrElements = 0;
	_iNrVertices = 0;

	// no data has been read initially
	_bMeshPresent = false;
	
	return;
}

void cl_ReadVistaFEM::_CleanUp(void)
{
	// clean up this instance
	_vviElements.resize(0);
	_vvdNodes.resize(0);
	_viLabels.resize(0);

	return;
}
	
