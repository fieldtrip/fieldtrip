#ifndef __VISTAPRIMITIVE_H__
#define __VISTAPRIMITIVE_H__

#include <string>
#include <vector>
#include <valarray>

// #include <VGraph.h>
#include "FIL_Vista.h"

class cl_ReadVistaFEM
{
	public:
		// standard constructor
		cl_ReadVistaFEM(void);

		// destructor
		~cl_ReadVistaFEM();

		// main initialization function
		int _ReadMesh(std::string sFilename);

		// get mesh
		int _GetMesh(std::vector<std::valarray<double> > &rvvdNodes,
				std::vector<std::valarray<int> > &rvviElements,
				std::valarray<int> &rviLabels);

	private:
		// init variables to default values, mark pointers as invalid
		void _InitVariables(void);

		// clean up class
		void _CleanUp(void);

		// status variables
		// was a mesh read successfully and does this class contain mesh data?
		bool _bMeshPresent;

		// mesh data
		// nodes of FE mesh
		std::vector<std::valarray<double> > _vvdNodes;
		// elements of FE mesh
		std::vector<std::valarray<int> > _vviElements;
		// labels of elements
		std::valarray<int> _viLabels;

		// number of vertices
		int _iNrVertices;
		// number of elements
		int _iNrElements;

};

class VPrimitive : public VNodeBaseRec
{
	public:
		VLong vcnt;
		VLong id[8];
		VPrimitive(void)
		{ hops = 0; visited = 0; head = 0; weight = 0; }
};

#endif /* __VISTAPRIMITIVE_H__ */
