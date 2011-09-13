/*
	$Log: FIL_Mesh.h,v $
	Revision 1.1.1.1  2003/08/21 13:09:55  anwander
	
	Final simbio release of the inverse toolbox coupled with NeuroFEM
	
	Revision 1.10  2002/10/01 14:02:03  wollny
	Another big lsit of changes by FK MT
	
	Revision 1.9  2002/09/23 11:11:33  wollny
	big patch by Dr. Kruggel
	
	Revision 1.8  2002/08/19 15:49:52  fingberg
	added nodal displacement fields to conversion routines
	
	Revision 1.7  2002/07/30 12:24:11  wollny
	add function to include a deformation vector image in vertices
	
*/

/*
 *  NeuroFEM license:
 *  =================
 *  Copyright (c) 2007 by 
 *  Dr.Carsten Wolters, Dr.Alfred Anwander, Dr.Matthias Duempelmann, 
 *  Dr.Thomas Knoesche, Dr. F. Kruggel
 *
 *  Permission is hereby granted, free of charge, to any person obtaining 
 *  a copy of the NeuroFEM software and associated documentation files 
 *  (the "Software"), to deal in the Software without restrictions, including 
 *  without limitations the rights to use, copy, modify, merge, publish, 
 *  distribute, sublicense, and/or sell copies of the Software, and to permit 
 *  persons to whom the Software is furnished to do so, subject to the following 
 *  conditions:
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 *  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 *  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
 *  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
 *  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN 
 *  THE SOFTWARE. THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE 
 *  SOFTWARE IS WITH YOU.
 *
 *  The above copyright notice, permission and warranty notices shall be 
 *  included in all copies of the Software. A copyright notice of the 
 *  contributing authors and the above permission and warranty notices 
 *  shall be included in all source code files of the Software. 
 *
 */


#ifndef FIL_MESH_H
#define FIL_MESH_H

#include <stdio.h>
#include <math.h>
#include <map>
#include <list>
#include <assert.h>
#include "FIL_Vista.h"
#include "FIL_Vec3.h"

class vertex : public VNodeBaseRec  {
public:
	VFloat		type;
	VFloat		x, y, z;
        VUByte		color;
	
	vertex(FIL_Vec3 v = 0) : x(v.x), y(v.y), z(v.z)
        		{ hops = 0; visited = 0; head = 0;
        		  weight = 0; type = 1; };
	vertex(double _x, double _y, double _z) : x(_x), y(_y), z(_z)
        		{ hops = 0; visited = 0; head = 0;
        		  weight = 0; type = 1; };
	FIL_Vec3		getPoint()
			{ return FIL_Vec3(x, y, z); };
	void		setPoint(FIL_Vec3 &v)
			{ x = v.x; y = v.y; z = v.z; };
	FIL_Vec3		center(VGraph vertices);
	void		smooth(VGraph vertices, double fac);
	inline bool	operator==(const vertex& v) const;
	void		setType(int t)
			{ type = t; };
	int		getType()
			{ return (int)type; };
	void		setColor(double c);
	double		getColor();
	FIL_Vec3		getRGB();
	void		setRGB(FIL_Vec3 &v);
};

inline bool vertex::operator==(const vertex& v) const
{
	double dx = x-v.x, dy = y-v.y, dz = z-v.z;
	return (dx*dx + dy*dy + dz*dz) < 1e-6;
}

class nvertex : public vertex  {
public:
	VFloat		nx, ny, nz;

	FIL_Vec3		getNormal()
			{ return FIL_Vec3(nx, ny, nz); };
	void		setNormal(FIL_Vec3 &v)
			{ nx = v.x; ny = v.y; nz = v.z; };
};



class fcvertex : public vertex  {
public:
        VFloat          fx, fy, fz;  // initial nodal force
        VFloat          flag;        // nodal contraint

        FIL_Vec3            getForce()
                        { return FIL_Vec3(fx, fy, fz); };
        void            setForce(FIL_Vec3 &v)
                        { fx = v.x; fy = v.y; fz = v.z; };

	int		getFlag()
			{ return (int)flag; };
	void		setFlag(int t)
			{ flag = t; };
};


class face : public VNodeBaseRec  {
public:
	VLong		vcnt;
	VLong		vtx[3];

 	face(int a, int b, int c)
       			{ hops = 0; visited = 0; head = 0; weight = 0;
			  vcnt = 3; vtx[0] = a; vtx[1] = b; vtx[2] = c; };
	VNode		getVertex(VGraph vertices, int i)
			{ return (i >= 0 && i < vcnt)?
			  VGraphGetNode(vertices, vtx[i]): 0; };
	FIL_Vec3		getPoint(VGraph vertices, int i)
			{ vertex *v = (vertex *)getVertex(vertices, i);
			  return v? v->getPoint(): (FIL_Vec3)0; };
	int		getType()
			{ return vcnt; };
        int	    	vertexPosition(int v) const
                	{ if (vtx[0] == v) return 0; if (vtx[1] == v) return 1;
                	  if (vtx[2] == v) return 2; return -1; };
	void		voxelize(VGraph vertices, VImage dst, FIL_Vec3 &dim);
	FIL_Vec3		computeNormal(VGraph vertices);
	void		boundingBox(VGraph vertices, FIL_Vec3 &bmin, FIL_Vec3 &bmax);
	double		cornerAngle(VGraph vertices, int i);
	int		remapVertex(int f, int t);
	double		area(VGraph vertices);
	double		perimeter(VGraph vertices);
	double		compactness (VGraph vertices);
	int		otherVertex (int v0, int v1);
};

class FIL_Mesh {
public:
	VGraph		vertices;
	VGraph		faces;
	VGraph		primitives;
	
        VImage          mp;  // element label vector
        VImage          pe;  // element partitioning vector
        VImage          pn;  // nodal   partitioning vector
        VImage          iv;  // source node vector
        VImage          tI;  // nodal constraints (tag image)
        VImage          fI;  // initial nodal forces (force image)
        VImage          sI;  // stress image
        VImage          dI;  // displacement image
	
	VAttrList	list;

	bool		resizeForNormals();
	bool		resizeForColor();
	void		computeNormals();
	bool		promoteMesh(int t);
	int		vertexType()
			{ vertex *v = (vertex *)VGraphFirstNode(vertices);
			  return v? v->getType(): 0; };
	void		compact();

public:
	FIL_Mesh(FILE *fp);
	FIL_Mesh(FILE *fp, int *myid, int *olevel);
	FIL_Mesh(VImage src, double lim);
	~FIL_Mesh();

        int             myid_;
        int             olevel_;

	bool		isComplete()   
			{ return vertices && faces; };
	void		voxelize(VImage dst);
	void		deform(VImage src, double sigma, double lim, double w1, double w2, double w3, int it);
	unsigned int	optimize(int pl, int wp, double bw, double cr, double mp, int ft);
	void		boundingBox(FIL_Vec3 &bmin, FIL_Vec3 &bmax);
	void		save(FILE *fp);
	void 		savevgraph(FILE *fp);
	void		colorize(VImage src);
	
	/** Add an image to the vertex-graph that describes the deformation of the FIL_Mesh
	    \param defo a 3D Vectorfield in the FIL_Vista format 
	*/
	void            add_defo(VImage defo);
	void		vista2drama(int *dim, int *dist, double* coord, int* enptr, int* conn,
				    int *gnn, double* force, int* constr, VAttrList *list,
                                    int *myid );
								
	void		getdim( int *dim, int *dist, int *myid, int *nprocs );
	void		vtoa( FILE* outf );
	void		vtogmv( FILE* outf );
	void		vtogpp( FILE* outf );
        VGraph          getPrimitives();
        VGraph          getVertices();
	void		smooth(double weight, int iter);
};



extern void linkFaceNeighbours(VGraph vertices, VGraph faces);
extern unsigned int countNodes(VGraph graph);
extern double getPixel(VImage image, FIL_Vec3 &v);

       void atov(FILE* inf, FILE* outf, int seq);
       void gpptov(FILE* inf, FILE* outf);

#endif
