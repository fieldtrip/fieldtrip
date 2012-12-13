#ifndef SIMBIO_MESH_H
#define SIMBIO_MESH_H

/*
  © Copyright 2003, C&C Research Laboratories, NEC Europe Ltd.
  © Copyright 2003, Max-Planck-Institute of Cognitive Neuroscience
  

    This file is part of SimBio-Vgrid.

    SimBio-Vgrid is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SimBio-Vgrid is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SimBio-Vgrid.  If not, see <http://www.gnu.org/licenses/>.


*/

/*
  $Id$
*/



#include <stdio.h>
#include <math.h>
#include <map>
#include <list>
#include <assert.h>
#include <vista.h>
#include <vec3.h>

class vertex : public VNodeBaseRec  {
public:
  VFloat		type;
  VFloat		x, y, z;
  VUByte		color;
	
  vertex(vec3 v = 0) : x(v.x), y(v.y), z(v.z)
  { hops = 0; visited = 0; head = 0;
  weight = 0; type = 1; };
  vertex(double _x, double _y, double _z) : x(_x), y(_y), z(_z)
  { hops = 0; visited = 0; head = 0;
  weight = 0; type = 1; };
  vec3		getPoint() const
  { return vec3(x, y, z); };
  void		setPoint(vec3 const&v)
  { x = v.x; y = v.y; z = v.z; };
  // barycenter of neighbors
  vec3		center(VGraph vertices) const;
  void		smooth(VGraph vertices, double fac);
  // barycenter of neighbors, using only neighbors on surface
  vec3		surface_center(VGraph vertices) const;
  vec3          smooth_surface_coord(VGraph vertices, double fac) const;
  vec3          smooth_volume_coord (VGraph vertices, double fac) const;

  inline bool	operator==(const vertex& v) const;
  void		setType(int t)   { type = t; };
  int		getType() const  { return (int)type; };
  void		setColor(double c);
  double	getColor() const;
  vec3		getRGB() const;
  void		setRGB(vec3 &v);
};

inline bool vertex::operator==(const vertex& v) const
{
  double dx = x-v.x, dy = y-v.y, dz = z-v.z;
  return (dx*dx + dy*dy + dz*dz) < 1e-6;
}

class nvertex : public vertex  {
public:
  VFloat		nx, ny, nz;

  vec3		getNormal()
  { return vec3(nx, ny, nz); };
  void		setNormal(vec3 &v)
  { nx = v.x; ny = v.y; nz = v.z; };
};



class fcvertex : public vertex  {
public:
  VFloat          fx, fy, fz;  // initial nodal force
  VFloat          flag;        // nodal contraint

  vec3            getForce()
  { return vec3(fx, fy, fz); };
  void            setForce(vec3 &v)
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
  vec3		getPoint(VGraph vertices, int i)
  { vertex *v = (vertex *)getVertex(vertices, i);
  return v? v->getPoint(): (vec3)0; };
  int		getType()
  { return vcnt; };
  int	    	vertexPosition(int v) const
  { if (vtx[0] == v) return 0; if (vtx[1] == v) return 1;
  if (vtx[2] == v) return 2; return -1; };
  void		voxelize(VGraph vertices, VImage dst, vec3 &dim);
  vec3		computeNormal(VGraph vertices);
  void		boundingBox(VGraph vertices, vec3 &bmin, vec3 &bmax);
  double		cornerAngle(VGraph vertices, int i);
  int		remapVertex(int f, int t);
  double		area(VGraph vertices);
  double		perimeter(VGraph vertices);
  double		compactness (VGraph vertices);
  int		otherVertex (int v0, int v1);
};

class mesh {
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
  mesh(FILE *fp, int *myid = NULL, int *olevel = NULL);
  // where is this constructor defined???
  mesh(VImage src, double lim);
  mesh(VGraph v, VGraph p, int olevel = 1);

  ~mesh();

  bool            owned; // are the data owned or just referenced?
  int             myid_;
  int             olevel_;

  bool		isComplete()   
  { return vertices && faces; };
  void		voxelize(VImage dst);
  void		boundingBox(vec3 &bmin, vec3 &bmax);
  void		save(FILE *fp);
  void 		savevgraph(FILE *fp);
  void		colorize(VImage src);
	
  /** Add an image to the vertex-graph that describes the deformation of the mesh
      \param defo a 3D Vectorfield in the vista format 
  */
  void          add_defo(VImage defo);
  void		vista2drama(int *dim, int *dist, double* coord, int* enptr, int* conn,
			    int *gnn, double* force, int* constr, VAttrList *list,
			    int *myid );
								
  void		getdim( int *dim, int *dist, int *myid, int *nprocs );
  void		vtoa( FILE* outf );
  void		vtogmv( FILE* outf, bool do_scale = false);
  void		vtodx ( FILE* outf );
  void		vtogpp( FILE* outf );
  VGraph        getPrimitives();
  VGraph        getVertices();
  void		smooth(double weight, int iter);
};



extern void linkFaceNeighbours(VGraph vertices, VGraph faces);
extern unsigned int countNodes(VGraph graph);
extern double getPixel(VImage image, vec3 &v);

extern void atov(FILE* inf, FILE* outf, int seq);
extern void gpptov(FILE* inf, FILE* outf);

#endif
