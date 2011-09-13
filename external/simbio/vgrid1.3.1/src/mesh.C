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
 *
 *  mesh.C: mesh handling routines, mesh I/O
 *
 *  $Id: mesh.C,v 1.29 2005/01/26 12:07:21 berti Exp $
 *
 */

#include <mesh.h>
#include <primitive.h>

#include <assert.h>   

#include <string>

using namespace std; 

extern vec3 getForce(VImage grad, const vec3 &v);


#define DIST(x)		(sqrt(fabs(x)))

unsigned int countNodes(VGraph graph)
{
  unsigned int n = 0;
  for (int i=1; i<=graph->lastUsed; i++)
    if (VGraphGetNode(graph, i)) n++;
  return n;
}

typedef multimap<VLong, VLong, less<VLong> > mmap;
typedef list<VLong> VLongList;

// this is a small helper class to implement a vertex cache
// based on the vertex position

static mmap *makeVertexToFaceMap(VGraph vertices, VGraph faces)
  // constructs a data structure for a faster access from vertices to faces
{
  mmap *vtof = new mmap;
  for (int i = 1; i <= faces->lastUsed; i++)  {
    VNode n = VGraphGetNode(faces, i);
    if (n == NULL) continue;
    VLong vcnt = ((VLong *)n->data)[0];
    for (int j = 1; j <= vcnt; j++)  {
      VLong d1 = ((VLong*)n->data)[j];
      VLong d2 = i;

      // avoiding duplicates
      bool found = false;
      for (mmap::iterator it = vtof->find(d1); it != vtof->end() && it->first == d1; it++)  {
	if (it->second == d2) { found = true; break; }
      }
      if (found == false) vtof->insert(pair<const VLong, VLong>(d1, d2));
    }
  }
  return vtof;
}

void linkFaceNeighbours(VGraph vertices, VGraph faces)
  // link polygons with a shared edge
{
  mmap *vtof = makeVertexToFaceMap(vertices, faces);
  VGraphClearVisit(vertices);
  for (int i = 1; i <= vertices->lastUsed; i++)  {
    VNode n = VGraphGetNode(vertices, i);
    if (n == NULL) continue;
    VNodeSetVisit(n);
    for (VAdjacency adj = n->base.head; adj; adj = adj->next)  {
      if (VNodeTestVisit(VGraphGetNode(vertices, adj->id))) continue;
      VLongList flist;
      for (mmap::iterator it1 = vtof->find(adj->id); (unsigned int)it1->first == adj->id; it1++)  {
	for (mmap::iterator it2 = vtof->find(i); it2->first == i; it2++)  {
	  if (it1->second == it2->second) flist.push_back(it1->second);
	}
      }
      if (flist.size() > 1)  {
	VNode f1 = VGraphGetNode(faces, flist.front());
	bool found = false;
	for (VAdjacency fadj = f1->base.head; fadj; fadj = fadj->next)  {
	  if (fadj->id == (unsigned int)flist.back())  { found = true; break; }
	}
	if (found == false) linkNodes(faces,flist.front(),flist.back());
      }
    }
  }
  delete vtof;
}

static double distanceToTriangle(vec3 &v0, vec3 &v1, vec3 &v2, vec3 &p)
  // given a triangle with vertices v0, v1 and v2, compute distance from p to triangle
  // believe me or not: this is a non-trivial problem
  // a solution was kindly provided by Dave Eberly (http://www.magic-software.com)
{
  vec3 e0 = v1-v0; vec3 e1 = v2-v0;
  vec3 diff = v0-p;
  double A = e0.dot(e0);
  double B = e0.dot(e1);
  double C = e1.dot(e1);
  double D = e0.dot(diff);
  double E = e1.dot(diff);
  double F = diff.dot(diff);
  double det = fabs(A*C-B*B);
  double s = B*E-C*D; double t = B*D-A*E;
  if (s+t <= det) {
    if (s < 0) {
      if (t < 0) {
	if (D < 0) {
	  if (-D >= A) return DIST(A+2*D+F);
	  else return DIST(F-D*D/A);
	} else {
	  if (E >= 0) return DIST(F);
	  else if (-E >= C) return DIST(C+2*E+F);
	  else return DIST(F-E*E/C);
	}
      } else {
	if (E >= 0) return DIST(F);
	else if (-E >= C) return DIST(C+2*E+F);
	else return DIST(F-E*E/C);
      }
    } else if (t < 0) {
      if (D >= 0) return DIST(F);
      else if (-D >= A) return DIST(A+2*D+F);
      else return DIST(F-D*D/A);
    } else {
      double invDet = 1.0/det;
      s *= invDet; t *= invDet;
      return DIST(s*(A*s+B*t+2*D)+t*(B*s+C*t+2*E)+F);
    }
  } else {
    double tmp0, tmp1, numer, denom;
    if (s < 0)  {
      tmp0 = B+D; tmp1 = C+E;
      if (tmp1 > tmp0) {
	numer = tmp1 - tmp0; denom = A-2*B+C;
	if (numer >= denom) return DIST(A+2*D+F);
	else {
	  s = numer/denom; t = 1-s;
	  return DIST(s*(A*s+B*t+2*D)+t*(B*s+C*t+2*E)+F);
	}
      } else {
	if (tmp1 <= 0) return DIST(C+2*E+F);
	else if (E >= 0) return DIST(F);
	else return DIST(F-E*E/C);
      }
    } else if (t < 0) {
      tmp0 = B+E; tmp1 = A+D;
      if (tmp1 > tmp0) {
	numer = tmp1 - tmp0; denom = A-2*B+C;
	if (numer >= denom) return DIST(C+2*E+F);
	else {
	  t = numer/denom; s = 1-t;
	  return DIST(s*(A*s+B*t+2*D)+t*(B*s+C*t+2*E)+F);
	}
      } else {
	if (tmp1 <= 0) return DIST(A+2*D+F);
	else if (D >= 0) return DIST(F);
	else return DIST(F-D*D/A);
      }
    } else {
      numer = C+E-B-D;
      if (numer <= 0) return DIST(C+2*E+F);
      else {
	denom = A-2*B+C;
	if (numer >= denom) return DIST(A+2*D+F);
	else {
	  s = numer/denom; t = 1-s;
	  return DIST(s*(A*s+B*t+2*D)+t*(B*s+C*t+2*E)+F);
	}
      }
    }
  }
}

vec3 vertex::center (VGraph vertices) const
  // computes the center of gravity around a vertex
{
  vec3 c = 0;
  int cnt = 0;
  for (VAdjacency adj = head; adj; adj = adj->next) {
    vertex *v = (vertex *) VGraphGetNode (vertices, adj->id);
    if (v == 0)
      continue;
    c += v->getPoint ();
    cnt++;
  };
  return (cnt ? c / cnt : getPoint ());
}

void vertex::smooth (VGraph vertices, double fac)
  // perform per-vertex smoothing operation
{
  vec3 c = center (vertices);
  vec3 v = (1 - fac) * getPoint () + fac * c;
  setPoint (v);
}

// computes the center of gravity around a vertex, using only surface vertices (color 1)
vec3 vertex::surface_center(VGraph vertices) const
{
  vec3 c = 0;
  int cnt = 0;
  for (VAdjacency adj = head; adj; adj = adj->next) {
    vertex *v = (vertex *) VGraphGetNode (vertices, adj->id);
    if (v == 0 || v->color != 1)
      continue;
    c += v->getPoint ();
    cnt++;
  };
  return (cnt ? c / cnt : getPoint ());
}

// return the new smoother coordinates of the vertex, without changing it
vec3 vertex::smooth_surface_coord(VGraph vertices, double fac) const {
  return (1 - fac) * getPoint () + fac * surface_center(vertices);
}

vec3 vertex::smooth_volume_coord(VGraph vertices, double fac) const {
  return (1 - fac) * getPoint () + fac * center(vertices);
}

double vertex::getColor() const
  // returns field c as a double
{
  int off;
  switch ((int)type)  {
  case 4: off =  4; break;
  case 5: off =  7; break;
  case 6: off = 10; break;
  default:
    return 0.0;
  };
  return (double)(((VFloat *)((VNode)this)->data)[off]);
}

void vertex::setColor(double c)
  // sets field c from a double
{
  int off;
  switch ((int)type)  {
  case 4: off =  4; break;
  case 5: off =  7; break;
  case 6: off = 10; break;
  default:
    return;
  };
  (((VFloat *)((VNode)this)->data)[off]) = (VFloat)c;
}

vec3 vertex::getRGB() const
  // returns fields (cr cg cb) as a vector
{
  int off;
  vec3 v = 0;

  switch ((int)type)  {
  case 7: off = 4; break;
  case 8: off = 7; break;
  default:
    return v;
  };
  v.x = ((VFloat *)((VNode)this)->data)[off];
  v.y = ((VFloat *)((VNode)this)->data)[off+1];
  v.z = ((VFloat *)((VNode)this)->data)[off+2];
  return v;
}

void vertex::setRGB(vec3 &v)
  // sets fields (cr cg cb) from a vector
{
  int off;
  switch ((int)type)  {
  case 7: off = 4; break;
  case 8: off = 7; break;
  default:
    return;
  };
  ((VFloat *)((VNode)this)->data)[off]   = (VFloat)v.x;
  ((VFloat *)((VNode)this)->data)[off+1] = (VFloat)v.y;
  ((VFloat *)((VNode)this)->data)[off+2] = (VFloat)v.z;
}


void face::boundingBox(VGraph vertices, vec3 &bmin, vec3 &bmax)
{
  //bmin = +HUGE; bmax = -HUGE;
  if(vcnt > 0) {
    bmin = getPoint(vertices,0);
    bmax = getPoint(vertices,0);
    for (int i = 0; i < vcnt; i++) {
      vec3 p = getPoint(vertices, i);
      if (p.x < bmin.x) bmin.x = p.x;
      if (p.y < bmin.y) bmin.y = p.y;
      if (p.z < bmin.z) bmin.z = p.z;
      if (p.x > bmax.x) bmax.x = p.x;
      if (p.y > bmax.y) bmax.y = p.y;
      if (p.z > bmax.z) bmax.z = p.z;
    };
  }
}

vec3 face::computeNormal(VGraph vertices)
{
  vec3 v1 = getPoint(vertices, 0);
  vec3 v2 = getPoint(vertices, 1);
  vec3 v3 = getPoint(vertices, 2);
  return ((v1-v2).cross(v1-v3)).normalize();
}

double face::cornerAngle(VGraph vertices, int i)
{
  int i_prev = (i==0)? 2: i-1;
  int i_next = (i==2)? 0: i+1;
  vec3 ci = getPoint(vertices, i);
  vec3 cp = getPoint(vertices, i_prev);
  vec3 cn = getPoint(vertices, i_next);
  vec3 ep = (cp-ci).normalize();
  vec3 en = (cn-ci).normalize();
  return acos(ep.dot(en));
}

int face::remapVertex(int f, int t)
{
  int n = 0;
  for (int i = 0; i < 3; i++)
    if (vtx[i] == f) { vtx[i] = t; n++; }
  return n;
}

void face::voxelize(VGraph vertices, VImage dst, vec3 &dim)
  // voxelize a polygon
{
  vec3 bmin, bmax;
  double delta = dim.scalar()/2;

  boundingBox(vertices, bmin, bmax);
  vec3 v0 = getPoint(vertices, 0);
  vec3 v1 = getPoint(vertices, 1);
  vec3 v2 = getPoint(vertices, 2);
  bmin.x = bmin.x/dim.x; if (bmin.x < 0) bmin.x = 0;
  bmin.y = bmin.y/dim.y; if (bmin.y < 0) bmin.y = 0;
  bmin.z = bmin.z/dim.z; if (bmin.z < 0) bmin.z = 0;
  bmax.x = bmax.x/dim.x; if (bmax.x >= VImageNColumns(dst)) bmax.x = VImageNColumns(dst)-1;
  bmax.y = bmax.y/dim.y; if (bmax.y >= VImageNRows(dst)) bmax.y = VImageNRows(dst)-1;
  bmax.z = bmax.z/dim.z; if (bmax.z >= VImageNBands(dst)) bmax.z = VImageNBands(dst)-1;
  for (int z = (int)bmin.z; z <= (int)bmax.z; z++)  {
    for (int y = (int)bmin.y; y <= (int)bmax.y; y++)  {
      for (int x = (int)bmin.x; x <= (int)bmax.x; x++)  {
	vec3 p(x, y, z);
	p = p*dim;
	if (distanceToTriangle(v0, v1, v2, p) > delta) continue;
	VSetPixel(dst, z, y, x, 1);
      };
    };
  };
}

double face::area(VGraph vertices)
{
  if (vcnt <= 2) return 0;

  if (vcnt == 3)  {
    vec3 v0 = getPoint(vertices, 0);
    vec3 v1 = getPoint(vertices, 1);
    vec3 v2 = getPoint(vertices, 2);
    vec3 a = v1 - v0;
    vec3 b = v2 - v0;
    vec3 c = a.cross(b);
    return fabs(0.5 * c.scalar());
  } else {
    // assuming a planar polygon 
    // taken from comp.graphics.algorithms faq item 2.01
    vec3 c = 0.0, nm = computeNormal(vertices);
    for (int i = 0; i < vcnt; i++)  {
      vec3 v1 = getPoint(vertices, i);
      vec3 v2 = getPoint(vertices, ((i+1) % vcnt));
      c += v1.cross(v2);
    };
    return fabs(0.5 * nm.dot(c));
  };
}

double face::perimeter(VGraph vertices)
{
  double p = 0.0;
  for (int i = 0; i < vcnt; i++)  {
    vec3 v1 = getPoint(vertices, i);
    vec3 v2 = getPoint(vertices, ((i+1) % vcnt));
    p += (v1-v2).scalar();
  };
  return p;
}

double face::compactness (VGraph vertices)
{
  double p = perimeter(vertices);
  double a = area(vertices);
  return a? p*p/(a*4*M_PI): 0;
}

int face::otherVertex (int v0, int v1)
  // given vertices v0 and v1, find the other vertex (of a triangle)
{
  if (vcnt != 3)
    return 0;
  for (int i = 0; i < vcnt; i++)
    if (vtx[i] != v0 && vtx[i] != v1)
      return vtx[i];
  return 0;		// should not happen
}

//  ============== mesh stuff ===========================================

void mesh::smooth (double weight, int iter)
{
	
  for (int i = 0; i < iter; i++) {
    for (vertex * v = (vertex *) VGraphFirstNode (vertices); v;
	 v = (vertex *) VGraphNextNode (vertices)){
      if(v->color == 1 ){
	v->smooth (vertices, weight);
      }
    }
  }
	
}

void mesh::boundingBox(vec3 &bmin, vec3 &bmax)
{
  //bmin = +HUGE; bmax = -HUGE;
  if(VGraphFirstNode(vertices)) {
    vertex *v = (vertex *)VGraphFirstNode(vertices);
    bmin = v->getPoint();
    bmax = v->getPoint();
    for (; v; v = (nvertex *)VGraphNextNode(vertices))  {
      vec3 p = v->getPoint();
      if (p.x < bmin.x) bmin.x = p.x;
      if (p.y < bmin.y) bmin.y = p.y;
      if (p.z < bmin.z) bmin.z = p.z;
      if (p.x > bmax.x) bmax.x = p.x;
      if (p.y > bmax.y) bmax.y = p.y;
      if (p.z > bmax.z) bmax.z = p.z;
    }
  }
}

void mesh::voxelize(VImage dst)
  // voxelize a mesh
{
  VStringConst attr;
  vec3 dim = 1;
  if (VGetAttr(VImageAttrList(dst), "voxel", 0, VStringRepn, &attr) == VAttrFound)
    sscanf(attr, "%lf%lf%lf", &dim.x, &dim.y, &dim.z);

  for (face *f = (face *)VGraphFirstNode(faces); f;
       f = (face *)VGraphNextNode(faces))
    f->voxelize(vertices, dst, dim);
}

mesh::mesh(VGraph v, VGraph p, int olevel)
{
  vertices   = v;
  primitives = p;
  owned = false;
  olevel_ = olevel;

  mp = NULL;
  pe = NULL;
  pn = NULL;
  iv = NULL;
  tI = NULL;
  fI = NULL;
  sI = NULL;
}

mesh::mesh(FILE *fp, int *myid, int *olevel)
  // create a mesh: read and check
{
  owned = true;
  // make this debugging stuff at least compatible with the old mesh version
  myid_ = myid ? *myid : 1;
  olevel_ = olevel ? *olevel : 0;

  if (fp == NULL) return;

  vertices   = NULL;
  primitives = NULL;

  mp = NULL;
  pe = NULL;
  pn = NULL;
  iv = NULL;
  tI = NULL;
  fI = NULL;
  sI = NULL;

  list = VReadFile(fp, 0);
  if (list == NULL) {
    VError(" mesh::mesh : Input file not found");
    return;
  }

  VAttrListPosn p;
  for (VFirstAttr(list, &p); VAttrExists(&p); VNextAttr(&p)) {
    if (VGetAttrRepn(&p) != VGraphRepn) continue;
    if (vertices == NULL)
      VGetAttrValue(&p, 0, VGraphRepn, &vertices);
    else if (primitives == NULL)
      VGetAttrValue(&p, 0, VGraphRepn, &primitives);
  };
  if ( !isComplete() ) return;
}

void mesh::savevgraph(FILE *fp)
{
  VWriteGraphs(fp, NULL, 1, &vertices);
  return;
}

void mesh::save(FILE *fp)
  // save mesh into file
{
  VGraph graph[2];
  //changed by UH (02/05/01)
  //vm & headfem convention
  graph[0] = vertices;
  graph[1] = primitives; 
  VWriteGraphs(fp, NULL, 2, graph);
}

VGraph mesh::getVertices()
{
  VGraph tmp;
  tmp = vertices; 
  return tmp;
}


VGraph mesh::getPrimitives()
{
  VGraph tmp;
  tmp = primitives;
  return tmp;
}


mesh::~mesh()
  // mesh destructor
{
  if(owned) {
    if (olevel_ > 0) {
      fprintf(stdout,"%d entering mesh destructor ... \n",myid_);
      fflush(stdout);

      VAttrList list;
      VAttrRec *a, *a_next;
      list = vertices->attributes;
      fprintf(stdout,"%d Attributes of vertices graph:\n",myid_);
      if (list) {
	fprintf(stdout,"%d name = <%s>\n",myid_,list->name);
	for (a = list->next; a; a = a_next) {
	  a_next = a->next;
	  fprintf(stdout,"%d name = <%s>\n",myid_,a->name);
	}
      }
	    
      list = primitives->attributes;
      fprintf(stdout,"%d Attributes of primitives graph:\n",myid_);
      if (list) {
	fprintf(stdout,"%d name = <%s>\n",myid_,list->name);
	for (a = list->next; a; a = a_next) {
	  a_next = a->next;
	  fprintf(stdout,"%d name = <%s>\n",myid_,a->name);
	}
      }
    }
	  
    if (vertices   != NULL) {
      if (olevel_ > 1) {
	fprintf(stdout,"%d mesh destructor: destroy vertices ... \n",myid_);
	fflush(stdout);
      }
      VDestroyGraph(vertices);
      if (olevel_ > 1) {
	fprintf(stdout,"%d mesh destructor: ... destroyed vertices \n",myid_);
	fflush(stdout);
      }
    }
	  
    if (primitives != NULL) {
      if (olevel_ > 1) {
	fprintf(stdout,"%d mesh destructor: destroy primitives ... \n",myid_);
	fflush(stdout);
      }
      VDestroyGraph(primitives);
      if (olevel_ > 1) {
	fprintf(stdout,"%d mesh destructor: ... destroyed primitives \n",myid_);
	fflush(stdout);
      }
    }
	  
    if (olevel_ > 0) {
      fprintf(stdout,"%d ... returning from mesh destructor. \n",myid_);
      fflush(stdout);
    }
  }
       
}


void mesh::computeNormals()
  // compute normals at vertices using face information
{
  vec3 nm = 0;
	
  // for all faces: compute vertex normals
  for (nvertex *n = (nvertex *)VGraphFirstNode(vertices); n;
       n = (nvertex *)VGraphNextNode(vertices))  {
    n->setNormal(nm);
  };
  for (face *f = (face *)VGraphFirstNode(faces); f;
       f = (face *)VGraphNextNode(faces))  {
    // we silently assume that we have planar polygons
    nm = f->computeNormal(vertices);
    int vcnt = f->getType();
    for (int i = 0; i < vcnt; i++)  {
      nvertex *n = (nvertex *)f->getVertex(vertices, i);
      vec3 ni = nm + n->getNormal();
      n->setNormal(ni);
    };
  };
  for (nvertex *n = (nvertex *)VGraphFirstNode(vertices); n;
       n = (nvertex *)VGraphNextNode(vertices))  {
    nm = -(n->getNormal()).normalize();
    n->setNormal(nm);
  };
}

bool mesh::promoteMesh(int t)
  // change vertex type in a mesh
  //
  // list of vertex types
  // 1 x y z
  // 2 x y z nx ny nz
  // 3 x y z nx ny nz kn kg km
  // 4 x y z c
  // 5 x y z nx ny nz c
  // 6 x y z nx ny nz kn kg km c
  // 7 x y z cr cg cb
  // 8 x y z nx ny nz cr cg cb
  //
{
  // promotion compatibility map
  // 0: illegal, 1: noop, 2: move color, 3: move RGB
  static int comp[][9] = {
    { 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 1, 1, 1, 1, 1, 1, 1, 1 },
    { 0, 0, 1, 1, 0, 1, 1, 0, 1 },
    { 0, 0, 0, 1, 0, 0, 1, 0, 0 },
    { 0, 0, 0, 0, 1, 2, 2, 0, 0 },
    { 0, 0, 0, 0, 0, 1, 2, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 1, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 1, 3 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 1 },
  };
  static int nf[] = { 0, 4, 7, 10, 5, 8, 11, 7, 10 };
  int type = vertexType();

  // incompatible promotion
  if (t > 8 || comp[type][t] == 0)  {
    fprintf(stderr, "illegal operation requested (%d -> %d)\n", type, t);
    return false;
  };

  // identical type: noop
  if (type == t) return true;

  // resize fields
  if (VGraphResizeFields(vertices, nf[t]) == 0) return false;

  // resize only - set new vertex type
  if (comp[type][t] == 1) {
    // move color value
    for (vertex *v = (vertex *)VGraphFirstNode(vertices); v;
	 v = (vertex *)VGraphNextNode(vertices)) v->setType(t);
    return true;
  };

  if (comp[type][t] == 2)  {
    // move color value
    for (vertex *v = (vertex *)VGraphFirstNode(vertices); v;
	 v = (vertex *)VGraphNextNode(vertices))  {
      double c = v->getColor(); v->setType(t); v->setColor(c);
    };
    return true;
  } else if (comp[type][t] == 3)  {
    // move RGB value
    for (vertex *v = (vertex *)VGraphFirstNode(vertices); v;
	 v = (vertex *)VGraphNextNode(vertices))  {
      vec3 c = v->getRGB(); v->setType(t); v->setRGB(c);
    };
    return true;
  };

  // we should not come here
  return false;
}

bool mesh::resizeForNormals()
{
  // entry is old type, return new type
  static int comp[] = { 0, 2, 0, 0, 5, 0, 0, 8, 0 };

  // find new type
  int type = comp[vertexType()];

  // see if we need to work on this...
  return (type? promoteMesh(type): true);
}

bool mesh::resizeForColor()
{
  // entry is old type, return new type
  static int comp[] = { 0, 4, 5, 6, 0, 0, 0, -1, -1 };

  // find new type
  int type = comp[vertexType()];
  if (type < 0)  {
    fprintf(stderr, "resizeForColor: illegal source type (%d)\n", vertexType());
    return false;
  };

  // see if we need to work on this...
  return (type? promoteMesh(type): true);
}



void mesh::getdim( int *dim, int *dist, int *myid, int *nprocs )
  // get array dimensions for DRAMAtool
{
  VLong 		nv=0, ne=0;
  VUByte		npart=1;

  assert(*nprocs >  0);
  assert(*myid   >= 0);
  assert(dim        != NULL);
  assert(dist       != NULL);
  assert(vertices   != NULL);
  assert(primitives != NULL);

  if (VGetAttr(VGraphAttrList(vertices), "nverts", 0, VLongRepn, &nv) != VAttrFound) {
    nv = countNodes(vertices);
    if (nv <= 0) {
      fprintf(stderr,"getdim:  ERROR number of vertices undefined! \n");
      return;
    }
  }
  if (VGetAttr(VGraphAttrList(vertices), "partition(s)", 0, VUByteRepn, &npart) != VAttrFound) {
    npart = 1;
  }
  if (VGetAttr(VGraphAttrList(primitives), "nelems", 0, VLongRepn, &ne) != VAttrFound) {
    ne = countNodes(primitives);
    if (ne <= 0) {
      fprintf(stderr," ERROR: number of primitives undefined! \n");
      return;
    }
  }

  fprintf(stderr,"%d getdim: nv=%d ne=%d \n",*myid, nv,ne);

  //get image pointers
  assert(VGetAttr(VGraphAttrList(primitives),"partelem", 0, VImageRepn, &pe)==VAttrFound);
  assert(VGetAttr(VGraphAttrList(primitives),"matprops", 0, VImageRepn, &mp)==VAttrFound);
  assert(VGetAttr(VGraphAttrList(vertices),  "partnode", 0, VImageRepn, &pn)==VAttrFound);

  if ( pn == NULL || pe == NULL || mp == NULL )  {
    fprintf(stderr, "getdim: ERROR found NULL image pointers.\n");
    return;
  }

  // count # of local nodes & elements
  int local_nv=0;
  int *cnt = new int[*nprocs];
  for (int i=0; i<(*nprocs); i++) {
    cnt[i] = 0;
  }
  for (int j=0; j<nv; j++){
    int pid=(int)VGetPixel(pn, 0, 0, j);
    cnt[ pid ]++;
  }
  local_nv = cnt[ *myid ];
  dist[0] = 0;
  for (int i=0; i<(*nprocs); i++) {
    dist[i+1] = dist[i] + cnt[i];
  }
  delete [] cnt;

  // determine :
  // 1. number of local elements
  // 2. min and max material ID
  // 3. size of connectivity array
  int local_ne = 0;
  int gmax = 0;
  int gmin = 0;
  dim[1] = 0;
  if (ne > 0) {
    int mid = (int)VGetPixel(mp, 0, 0, 0);
    gmax = gmin = mid;
		
    for (int j=0; j<ne; j++){
      int pid = (int)VGetPixel(pe, 0, 0, j);
      if (pid == (*myid)) local_ne++;
      int mid = (int)VGetPixel(mp, 0, 0, j);
      if (mid > gmax) gmax=mid;
      if (mid < gmin) gmin=mid;

      primitive *p = (primitive *)VGraphGetNode(primitives, j+1);
      if (p != NULL) {
	if (pid == (*myid)) dim[1] += p->vcnt;
      }

    }
  }

  dim[0] = local_ne;
  dim[2] = local_ne+1;
  dim[3] = local_nv;
  dim[4] = 3;
  dim[5] = gmax-gmin+1;
  dim[6] = gmin;
	
  return;
}


void atov(FILE *inf, FILE *outf, int seq)
  // ascii to vista format conversion
  // not part of mesh class!
{
  int  np=1, ne, nn, n, p, i, in, m;
  char str1[16], str2[16], str3[16];
  bool implicit_links=true;

  fprintf(stdout," entering atov ... \n");
  fflush(stdout);

  //read ASCII mesh header
  fscanf(inf,"%i %i %i", &nn, &ne, &np);

  fscanf(inf,"%s %s %s", str1, str2, str3);

  if (strstr(str1,"version") != NULL) {
    fprintf(stdout,"%s %s %s\n", str1, str2, str3);
    fflush(stdout);
  }
  else {
    rewind( inf );
    fscanf(inf,"%i %i %i", &nn, &ne, &np);
  }

  if (seq > 0) np = 1;

  //generate  graphs:
  //          vertices   --> nodes    (4 fields: type=1, x,y,z)
  //          primitives --> elements (9 fields)

  VGraph vertices = VCreateGraph(nn, 4, VFloatRepn, false);
  if (vertices == NULL)  {
    fprintf(stderr, "atov: ERROR cannot create vertices graph.\n");
    return;
  }

  VGraph primitives = VCreateGraph(ne, 9, VLongRepn,  false);
  if (primitives == NULL)  {
    fprintf(stderr, "atov: ERROR cannot create primitives graph.\n");
    return;
  }

  //create images
  VImage mp = VCreateImage( 1, 1, ne, VUByteRepn );  // element label vector
  VImage pe = VCreateImage( 1, 1, ne, VUByteRepn );  // element partitioning vector
  VImage sI = VCreateImage( 6, 1, ne, VFloatRepn );  // stress image
  VImage pn = VCreateImage( 1, 1, nn, VUByteRepn );  // nodal   partitioning vector
  VImage iv = VCreateImage( 1, 1, nn, VBitRepn   );  // source node vector
  VImage tI = VCreateImage( 1, 1, nn, VUByteRepn );  // nodal constraints (tag image)
  VImage fI = VCreateImage( 3, 1, nn, VFloatRepn );  // initial nodal forces (force image)
  VImage dI = VCreateImage( 3, 1, nn, VFloatRepn );  // nodal displacements

  if ( mp == 0 ||  iv == 0 || pn == 0 || pe == 0 || tI == 0 || fI == 0 || sI == 0 )  {
    fprintf(stderr, "atov: ERROR cannot create attribute images\n");
    return;
  }

  VImageNFrames( mp )     = 1;
  VImageNComponents( mp ) = 1;

  VImageNFrames( pe )     = 1;
  VImageNComponents( pe ) = 1;

  VImageNFrames( sI )     = 1;
  VImageNComponents( sI ) = 6;

  VImageNFrames( pn )     = 1;
  VImageNComponents( pn ) = 1;

  VImageNFrames( iv )     = 1;
  VImageNComponents( iv ) = 1;

  VImageNFrames( tI )     = 1;
  VImageNComponents( tI ) = 1;

  VImageNFrames( fI )     = 1;
  VImageNComponents( fI ) = 3;

  VImageNFrames( dI )     = 1;
  VImageNComponents( dI ) = 3;

  //
  // nodes
  //

  fprintf(stdout," %d nodes ... \n",nn);
  fflush(stdout);

  for(i=0; i<nn; i++){
    float  x,y,z, fx,fy,fz;
    int    flag;

    fscanf(inf,"%i %i %i %g %g %g %i %g %g %g",
	   &n, &p, &in,
	   &x, &y, &z,
	   &flag,
	   &fx, &fy, &fz);

    if (seq > 0) p = 0;

    assert(i == n);                            // assume vertices are ordered
    // by global node number.
    VSetPixel(pn, 0, 0, i, (VUByte)p   );      // node partition
    VSetPixel(iv, 0, 0, i, (VBit)in    );      // src nodes
    VSetPixel(tI, 0, 0, i, (VUByte)flag);      // nodal constraints
    VPixel(fI, 0, 0, i, VFloat) = (VFloat)fx;  // initial nodal forces
    VPixel(fI, 1, 0, i, VFloat) = (VFloat)fy;
    VPixel(fI, 2, 0, i, VFloat) = (VFloat)fz;

    vec3 d = 0;
    d.x = x;
    d.y = y;
    d.z = z;
    vertex u( d );

    VGraphAddAndGrow(vertices, (VNode)&u, i+1);
  }

  //
  // elements
  //

  fprintf(stdout," %d elements ... \n",ne);
  fflush(stdout);

  for (i=0; i<ne; i++) {
    primitive e;
    int       id[8];

    fscanf(inf,"%i %i %i", &n, &p, &m );

    if (seq > 0) p = 0;

    //fprintf(stdout," %i %i %i \n",n,p,m);
    //fflush(stdout);

    // FE-type: m < 100 -> tet    m >= 100 -> hex
    if (m < 100) {
      fscanf(inf,"%i %i %i %i",
	     &id[0], &id[1], &id[2], &id[3]);

      e.vcnt  = 4;
      e.id[0] = id[0]+1;
      e.id[1] = id[1]+1;
      e.id[2] = id[2]+1;
      e.id[3] = id[3]+1;

      if (! implicit_links) {
	linkNodes(vertices, e.id[0], e.id[1]);
	linkNodes(vertices, e.id[0], e.id[2]);
	linkNodes(vertices, e.id[0], e.id[3]);
	linkNodes(vertices, e.id[1], e.id[2]);
	linkNodes(vertices, e.id[1], e.id[3]);
	linkNodes(vertices, e.id[2], e.id[3]);
      }
    }
    else {                     // internal connectivity: Hughes convention
      fscanf(inf,"%i %i %i %i %i %i %i %i",
	     &id[0], &id[1], &id[2], &id[3],
	     &id[4], &id[5], &id[6], &id[7]);

      e.vcnt  = 8;
      e.id[0] = id[0]+1;
      e.id[1] = id[1]+1;
      e.id[2] = id[2]+1;
      e.id[3] = id[3]+1;
      e.id[4] = id[4]+1;
      e.id[5] = id[5]+1;
      e.id[6] = id[6]+1;
      e.id[7] = id[7]+1;

      if (! implicit_links) {
	linkNodes(vertices, e.id[0], e.id[1]);
	linkNodes(vertices, e.id[1], e.id[2]);
	linkNodes(vertices, e.id[2], e.id[3]);
	linkNodes(vertices, e.id[3], e.id[0]);
	linkNodes(vertices, e.id[0], e.id[4]);
	linkNodes(vertices, e.id[1], e.id[5]);
	linkNodes(vertices, e.id[2], e.id[6]);
	linkNodes(vertices, e.id[3], e.id[7]);
	linkNodes(vertices, e.id[4], e.id[5]);
	linkNodes(vertices, e.id[5], e.id[6]);
	linkNodes(vertices, e.id[6], e.id[7]);
	linkNodes(vertices, e.id[7], e.id[4]);
      }
    }

    VGraphAddAndGrow(primitives, (VNode)&e, i+1);

    VSetPixel( mp, 0, 0, i, (VUByte)m );    // element label
    VSetPixel( pe, 0, 0, i, (VUByte)p );    // element partition
  }

  //
  // stress
  //
  if (feof(inf) == 0) {
    fscanf(inf,"%s", str1 );
    if (str1 == "stress") {
      fprintf(stdout," %d stresses ... \n",ne);
      fflush(stdout);
      for (i=0; i<ne; i++) {
	int   n;
	float sxx,syy,szz,sxy,sxz,syz;

	fscanf(inf,"%i %g %g %g %g %g %g",
	       &n, &sxx, &syy, &szz, &sxy, &sxz, &syz);

	assert( i == n );

	VPixel(sI, 0, 0, i, VFloat) = (VFloat)sxx;
	VPixel(sI, 1, 0, i, VFloat) = (VFloat)syy;
	VPixel(sI, 2, 0, i, VFloat) = (VFloat)szz;
	VPixel(sI, 3, 0, i, VFloat) = (VFloat)sxy;
	VPixel(sI, 4, 0, i, VFloat) = (VFloat)sxz;
	VPixel(sI, 5, 0, i, VFloat) = (VFloat)syz;
      }
    }
  }

  //
  // displacement
  //
  if (feof(inf) == 0) {
    fscanf(inf,"%s", str1 );
    if (str1 == "displacement") {
      fprintf(stdout," %d displacements ... \n",nn);
      fflush(stdout);
      for (i=0; i<nn; i++) {
	int   n;
	float dx,dy,dz;

	fscanf(inf,"%i %g %g %g", &n, &dx, &dy, &dz);

	assert( i == n );

	VPixel(dI, 0, 0, i, VFloat) = (VFloat)dx;
	VPixel(dI, 1, 0, i, VFloat) = (VFloat)dy;
	VPixel(dI, 2, 0, i, VFloat) = (VFloat)dz;
      }
    }
  }

  // patient name
  VSetAttr(VGraphAttrList(vertices), "patient", NULL, VStringRepn, "ascii to vista");

  // examination date
  VSetAttr(VGraphAttrList(vertices), "date", NULL, VStringRepn, "19 Aug 2002");

  // copy convention
  VSetAttr(VGraphAttrList(vertices), "convention", NULL, VStringRepn, "natural");

  // copy orientation
  VSetAttr(VGraphAttrList(vertices), "orientation", NULL, VStringRepn, "axial");

  // see section 4.2
  VSetAttr(VGraphAttrList(vertices),   "component_interp", NULL, VStringRepn, "vertex"   );
  VSetAttr(VGraphAttrList(primitives), "component_interp", NULL, VStringRepn, "primitive");

  // see section 4.7
  VSetAttr(VGraphAttrList(primitives), "primitive_interp", NULL, VStringRepn, "volume");

  // additional attribute # of partitions
  VSetAttr(VGraphAttrList(vertices),   "partition(s)",     NULL, VUByteRepn, np);

  // implicit links
  if(implicit_links)
    VSetAttr(VGraphAttrList(primitives), "implicit_links", NULL, VStringRepn, "true");
  else
    VSetAttr(VGraphAttrList(primitives), "implicit_links", NULL, VStringRepn, "false");

  //image attributes
  VAppendAttr(VImageAttrList(sI), "component_repn",   0, VStringRepn, "tensor6");
  VAppendAttr(VImageAttrList(sI), "component_interp", 0, VStringRepn, "stress E");

  VAppendAttr(VImageAttrList(fI), "component_repn",   0, VStringRepn, "vector3");
  VAppendAttr(VImageAttrList(fI), "component_interp", 0, VStringRepn, "force N");

  VAppendAttr(VImageAttrList(dI), "component_repn",   0, VStringRepn, "vector3");
  VAppendAttr(VImageAttrList(dI), "component_interp", 0, VStringRepn, "deformation");

  VAppendAttr(VImageAttrList(tI), "component_repn",   0, VStringRepn, "scalar" );
  VAppendAttr(VImageAttrList(tI), "component_interp", 0, VStringRepn, "tag"    );

  VAppendAttr(VImageAttrList(iv), "component_repn",   0, VStringRepn, "scalar" );
  VAppendAttr(VImageAttrList(iv), "component_interp", 0, VStringRepn, "src node tag");

  VAppendAttr(VImageAttrList(pn), "component_repn",   0, VStringRepn, "scalar" );
  VAppendAttr(VImageAttrList(pn), "component_interp", 0, VStringRepn, "nodal partition");

  VAppendAttr(VImageAttrList(mp), "component_repn",   0, VStringRepn, "scalar" );
  VAppendAttr(VImageAttrList(mp), "component_interp", 0, VStringRepn, "element label");

  VAppendAttr(VImageAttrList(pe), "component_repn",   0, VStringRepn, "scalar" );
  VAppendAttr(VImageAttrList(pe), "component_interp", 0, VStringRepn, "element partition");

  //graph attributes
  VAppendAttr(VGraphAttrList(vertices),   "image",      NULL, VImageRepn, tI);
  VAppendAttr(VGraphAttrList(vertices),   "forceimage", NULL, VImageRepn, fI);
  VAppendAttr(VGraphAttrList(vertices),   "deformation",NULL, VImageRepn, dI);
  VAppendAttr(VGraphAttrList(vertices),   "srcnodes",   NULL, VImageRepn, iv);
  VAppendAttr(VGraphAttrList(vertices),   "partnode",   NULL, VImageRepn, pn);
  VAppendAttr(VGraphAttrList(primitives), "matprops",   NULL, VImageRepn, mp);
  VAppendAttr(VGraphAttrList(primitives), "partelem",   NULL, VImageRepn, pe);
  VAppendAttr(VGraphAttrList(primitives), "stressimage",NULL, VImageRepn, sI);

  //save graphs
  VGraph g[2];
  //make sure that the graph node numbering is continous
  assert(vertices->lastUsed   == (int)countNodes(vertices));
  assert(primitives->lastUsed == (int)countNodes(primitives));
  //dirty --- cut off the unused tail
  VWarning("gridGenerator::write : dirty --- cut off the unused tail of VISTA graph!");
  vertices->size     = vertices->lastUsed;
  vertices->nnodes   = vertices->lastUsed;
  primitives->size   = primitives->lastUsed;
  primitives->nnodes = primitives->lastUsed;
  //vm & headfem convention
  g[0] = vertices;
  g[1] = primitives;
  VWriteGraphs(outf, NULL, 2, g);

  VDestroyImage(pn);
  VDestroyImage(pe);
  VDestroyImage(mp);
  VDestroyImage(iv);
  VDestroyImage(tI);
  VDestroyImage(fI);
  VDestroyImage(sI);
  VDestroyImage(dI);

  fprintf(stdout," ... done atov \n");
  fflush(stdout);

  return;

};



void gpptov(FILE *inf, FILE *outf)
  // gpp to vista format conversion
  // not part of mesh class!
{
  int  np=1, ne, n1, n2, n3, n4, nn;
  bool implicit_links=true;

  fprintf(stdout," entering gpptov ... \n");
  fflush(stdout);

  //read ASCII mesh header
  fscanf(inf,"%i %i %i %i %i", &nn, &n1, &n2, &n3, &n4);
  ne = n3 + n4;

  //generate  graphs:
  //          vertices   --> nodes    (4 fields: type=1, x,y,z)
  //          primitives --> elements (9 fields)

  VGraph vertices = VCreateGraph(nn, 4, VFloatRepn, false);
  if (vertices == NULL)  {
    fprintf(stderr, "gpptov: ERROR cannot create vertices graph.\n");
    return;
  }

  VGraph primitives = VCreateGraph(ne, 9, VLongRepn,  false);
  if (primitives == NULL)  {
    fprintf(stderr, "gpptov: ERROR cannot create primitives graph.\n");
    return;
  }

  //create images
  VImage mp = VCreateImage( 1, 1, ne, VUByteRepn );  // element label vector
  VImage pe = VCreateImage( 1, 1, ne, VUByteRepn );  // element partitioning vector
  VImage pn = VCreateImage( 1, 1, nn, VUByteRepn );  // nodal   partitioning vector
  VImage iv = VCreateImage( 1, 1, nn, VBitRepn   );  // source node vector
  VImage tI = VCreateImage( 1, 1, nn, VUByteRepn );  // nodal constraints (tag image)
  VImage fI = VCreateImage( 3, 1, nn, VFloatRepn );  // initial nodal forces (force image)

  if ( mp == 0 ||  iv == 0 || pn == 0 || pe == 0 || tI == 0 || fI == 0 )  {
    fprintf(stderr, "gpptov: ERROR cannot create attribute images\n");
    return;
  }

  VImageNFrames( mp )     = 1;
  VImageNComponents( mp ) = 1;

  VImageNFrames( pe )     = 1;
  VImageNComponents( pe ) = 1;

  VImageNFrames( pn )     = 1;
  VImageNComponents( pn ) = 1;

  VImageNFrames( iv )     = 1;
  VImageNComponents( iv ) = 1;

  VImageNFrames( tI )     = 1;
  VImageNComponents( tI ) = 1;

  VImageNFrames( fI )     = 1;
  VImageNComponents( fI ) = 3;

  //
  // nodes
  //

  fprintf(stdout," %d nodes ... \n",ne);
  fflush(stdout);

  for (int i=0; i<nn; i++){
    float  x,y,z;
    VFloat fx=0.0,fy=0.0,fz=0.0;
    VUByte flag=0, in=0;
    VUByte p=0;

    fscanf(inf,"%g %g %g", &x, &y, &z);

    VSetPixel(pn, 0, 0, i, p );        // node partition
    VSetPixel(iv, 0, 0, i, in  );      // src nodes
    VSetPixel(tI, 0, 0, i, flag);      // nodal constraints
    VPixel(fI, 0, 0, i, VFloat) = fx;  // initial nodal forces
    VPixel(fI, 1, 0, i, VFloat) = fy;
    VPixel(fI, 2, 0, i, VFloat) = fz;

    vec3 d = 0;
    d.x = x;
    d.y = y;
    d.z = z;
    vertex u( d );

    VGraphAddAndGrow(vertices, (VNode)&u, i+1);
  }

  //
  // tetrahedral elements
  //

  if (n3 > 0) {
    fprintf(stdout," %d tetrahedral elements ... \n",n3);
    fflush(stdout);
  }

  for (int i=0; i<n3; i++){
    primitive e;
    int       id[4];

    fscanf(inf,"%i %i %i %i", &id[0], &id[1], &id[2], &id[3]);

    e.vcnt  = 4;
    e.id[0] = id[0]+1;
    e.id[1] = id[1]+1;
    e.id[2] = id[2]+1;
    e.id[3] = id[3]+1;

    if (! implicit_links) {
      linkNodes(vertices, e.id[0], e.id[1]);
      linkNodes(vertices, e.id[0], e.id[2]);
      linkNodes(vertices, e.id[0], e.id[3]);
      linkNodes(vertices, e.id[1], e.id[2]);
      linkNodes(vertices, e.id[1], e.id[3]);
      linkNodes(vertices, e.id[2], e.id[3]);
    }

    VGraphAddAndGrow(primitives, (VNode)&e, i+1);

    VUByte m=0, p=0;
    VSetPixel( mp, 0, 0, i, (VUByte)m );    // element label
    VSetPixel( pe, 0, 0, i, (VUByte)p );    // element partition
  }

  //
  // hexahedral elements
  //

  if (n4 > 0) {
    fprintf(stdout," %d hexahedral elements ... \n",n4);
    fflush(stdout);
  }

  for (int i=0; i<n4; i++){
    primitive e;
    int       id[8];

    fscanf(inf,"%i %i %i %i %i %i %i %i",
	   &id[0], &id[1], &id[2], &id[3],
	   &id[4], &id[5], &id[6], &id[7]);

    e.vcnt  = 8;
    e.id[0] = id[0]+1;
    e.id[1] = id[1]+1;
    e.id[2] = id[2]+1;
    e.id[3] = id[3]+1;
    e.id[4] = id[4]+1;
    e.id[5] = id[5]+1;
    e.id[6] = id[6]+1;
    e.id[7] = id[7]+1;

    if (! implicit_links) {
      linkNodes(vertices, e.id[0], e.id[1]);
      linkNodes(vertices, e.id[1], e.id[2]);
      linkNodes(vertices, e.id[2], e.id[3]);
      linkNodes(vertices, e.id[3], e.id[0]);
      linkNodes(vertices, e.id[0], e.id[4]);
      linkNodes(vertices, e.id[1], e.id[5]);
      linkNodes(vertices, e.id[2], e.id[6]);
      linkNodes(vertices, e.id[3], e.id[7]);
      linkNodes(vertices, e.id[4], e.id[5]);
      linkNodes(vertices, e.id[5], e.id[6]);
      linkNodes(vertices, e.id[6], e.id[7]);
      linkNodes(vertices, e.id[7], e.id[4]);
    }
  

    int ii = i + n3;

    VGraphAddAndGrow(primitives, (VNode)&e, ii+1);

    VUByte m=100, p=0;
    VSetPixel( mp, 0, 0, ii, (VUByte)m );    // element label
    VSetPixel( pe, 0, 0, ii, (VUByte)p );    // element partition
  }

  // patient name
  VSetAttr(VGraphAttrList(vertices), "patient", NULL, VStringRepn, "gpp to vista");

  // examination date
  VSetAttr(VGraphAttrList(vertices), "date", NULL, VStringRepn, "23 Oct 2002");

  // copy convention
  VSetAttr(VGraphAttrList(vertices), "convention", NULL, VStringRepn, "natural");

  // copy orientation
  VSetAttr(VGraphAttrList(vertices), "orientation", NULL, VStringRepn, "axial");

  // see section 4.2
  VSetAttr(VGraphAttrList(vertices),   "component_interp", NULL, VStringRepn, "vertex"   );
  VSetAttr(VGraphAttrList(primitives), "component_interp", NULL, VStringRepn, "primitive");

  // see section 4.7
  VSetAttr(VGraphAttrList(primitives), "primitive_interp", NULL, VStringRepn, "volume");

  // additional attribute # of partitions
  VSetAttr(VGraphAttrList(vertices),   "partition(s)",     NULL, VUByteRepn, np);

  // implicit links
  if(implicit_links)
    VSetAttr(VGraphAttrList(primitives), "implicit_links", NULL, VStringRepn, "true");
  else
    VSetAttr(VGraphAttrList(primitives), "implicit_links", NULL, VStringRepn, "false");

  //image attributes
  VAppendAttr(VImageAttrList(fI), "component_repn",   0, VStringRepn, "vector3");
  VAppendAttr(VImageAttrList(fI), "component_interp", 0, VStringRepn, "force N");

  VAppendAttr(VImageAttrList(tI), "component_repn",   0, VStringRepn, "scalar" );
  VAppendAttr(VImageAttrList(tI), "component_interp", 0, VStringRepn, "tag"    );

  VAppendAttr(VImageAttrList(iv), "component_repn",   0, VStringRepn, "scalar" );
  VAppendAttr(VImageAttrList(iv), "component_interp", 0, VStringRepn, "src node tag");

  VAppendAttr(VImageAttrList(pn), "component_repn",   0, VStringRepn, "scalar" );
  VAppendAttr(VImageAttrList(pn), "component_interp", 0, VStringRepn, "nodal partition");

  VAppendAttr(VImageAttrList(mp), "component_repn",   0, VStringRepn, "scalar" );
  VAppendAttr(VImageAttrList(mp), "component_interp", 0, VStringRepn, "element label");

  VAppendAttr(VImageAttrList(pe), "component_repn",   0, VStringRepn, "scalar" );
  VAppendAttr(VImageAttrList(pe), "component_interp", 0, VStringRepn, "element partition");

  //graph attributes
  VAppendAttr(VGraphAttrList(vertices),   "image",      NULL, VImageRepn, tI);
  VAppendAttr(VGraphAttrList(vertices),   "forceimage", NULL, VImageRepn, fI);
  VAppendAttr(VGraphAttrList(vertices),   "srcnodes",   NULL, VImageRepn, iv);
  VAppendAttr(VGraphAttrList(vertices),   "partnode",   NULL, VImageRepn, pn);
  VAppendAttr(VGraphAttrList(primitives), "matprops",   NULL, VImageRepn, mp);
  VAppendAttr(VGraphAttrList(primitives), "partelem",   NULL, VImageRepn, pe);

  //save graphs
  VGraph g[2];
  //make sure that the graph node numbering is continous
  assert(vertices->lastUsed   == (int)countNodes(vertices));
  assert(primitives->lastUsed == (int)countNodes(primitives));
  //dirty --- cut off the unused tail
  VWarning("gridGenerator::write : dirty --- cut off the unused tail of VISTA graph!");
  vertices->size     = vertices->lastUsed;
  vertices->nnodes   = vertices->lastUsed;
  primitives->size   = primitives->lastUsed;
  primitives->nnodes = primitives->lastUsed;
  //vm & headfem convention
  g[0] = vertices;
  g[1] = primitives;
  VWriteGraphs(outf, NULL, 2, g);

  VDestroyImage(pn);
  VDestroyImage(pe);
  VDestroyImage(mp);
  VDestroyImage(iv);
  VDestroyImage(tI);
  VDestroyImage(fI);

  fprintf(stdout," ... done gpptov \n");
  fflush(stdout);

  return;

};


void mesh::vtoa( FILE *outf )
  // vista to ascii format conversion
{
  VLong 		nv=0, ne=0;
  VUByte          npart=1;

  fprintf(stdout," entering vtoa ... \n");
  fflush(stdout);

  assert(vertices   != NULL);
  assert(primitives != NULL);

  if (VGetAttr(VGraphAttrList(vertices), "nverts", 0, VLongRepn, &nv) != VAttrFound) {
    nv = countNodes(vertices);
    if (nv <= 0) {
      fprintf(stderr,"vtoa: ERROR number of vertices undefined! \n");
      return;
    }
  }
  if (VGetAttr(VGraphAttrList(vertices), "partition(s)", 0, VUByteRepn, &npart) != VAttrFound) {
    npart = 1;
  }
  if (VGetAttr(VGraphAttrList(primitives), "nelems", 0, VLongRepn, &ne) != VAttrFound) {
    ne = countNodes(primitives);
    if (ne <= 0) {
      fprintf(stderr,"vtoa: ERROR number of primitives undefined! \n");
      return;
    }
  }

  fprintf(outf,"%i %i %i     version date: 19.8.2002 \n",(int)nv,(int)ne,(int)npart);
  fflush(stdout);

  //attribute image pointers
  assert(VGetAttr(VGraphAttrList(primitives), "partelem", 0, VImageRepn, &pe)==VAttrFound);
  assert(VGetAttr(VGraphAttrList(vertices),   "partnode", 0, VImageRepn, &pn)==VAttrFound);
  assert(VGetAttr(VGraphAttrList(primitives), "matprops", 0, VImageRepn, &mp)==VAttrFound);

  if (VGetAttr(VGraphAttrList(primitives), "stressimage", 0, VImageRepn, &sI)!=VAttrFound) {
    sI = NULL;
  }
  if (VGetAttr(VGraphAttrList(vertices),   "srcnodes",    0, VImageRepn, &iv)!=VAttrFound) {
    iv = NULL;
  }
  if (VGetAttr(VGraphAttrList(vertices),   "image",       0, VImageRepn, &tI)!=VAttrFound) {
    tI = NULL;
  }
  if (VGetAttr(VGraphAttrList(vertices),   "forceimage",  0, VImageRepn, &fI)!=VAttrFound) {
    fI = NULL;
  }
  if (VGetAttr(VGraphAttrList(vertices),   "deformation", 0, VImageRepn, &dI)!=VAttrFound) {
    dI = NULL;
  }

  if ( mp == NULL || pn == NULL || pe == NULL )  {
    fprintf(stderr, "vtoa: ERROR found NULL image pointers.\n");
    return;
  }

  fprintf(stdout," %d nodes ... \n",nv);
  fflush(stdout);

  //write nodes
  for (int i=0; i<nv; i++) {
    vertex *v = (vertex *)VGraphGetNode(vertices, i+1);
    if (v != NULL) {
      int    mic, fid;
      double fx, fy, fz;

      double x = (double)v->x;
      double y = (double)v->y;
      double z = (double)v->z;

      int pid = (int)VGetPixel(pn, 0, 0, i);
      if (iv != NULL) {
	mic = (int)VGetPixel(iv, 0, 0, i);
      }
      else {
	mic = 0;
      }
      if (tI != NULL) {
	fid = (int)VGetPixel(tI, 0, 0, i);
      }
      else {
	fid = 0;
      }

      if (fI != NULL) {
	fx = (double)VGetPixel(fI, 0, 0, i);
	fy = (double)VGetPixel(fI, 1, 0, i);
	fz = (double)VGetPixel(fI, 2, 0, i);
      }
      else {
	fx = 0.0;
	fy = 0.0;
	fz = 0.0;
      }

      fprintf(outf,"%i %i %i  %15.7le %15.7le %15.7le  %i  %15.7le %15.7le %15.7le\n",
	      i,pid,mic, x, y, z,
	      fid, fx, fy, fz);
    }
  }

  fprintf(stdout," %d elements ... \n",ne);
  fflush(stdout);

  //write elements
  for (int i=0; i<ne; i++) {
    primitive *p = (primitive *)VGraphGetNode(primitives, i+1);
    if (p != NULL) {
      int pid = (int)VGetPixel(pe, 0, 0, i);
      int mid = (int)VGetPixel(mp, 0, 0, i);

      fprintf(outf,"%i %i %i  ", i, pid, mid);
      for (int j=0; j<p->vcnt; j++) {
	fprintf(outf,"%i ", p->id[j] - 1);
      }
      fprintf(outf,"\n");
    }
  }

  if (sI != NULL) {
    fprintf(stdout," %d stresses ...\n",ne);
    fflush(stdout);

    fprintf(outf,"stress\n");
    for (int i=0; i<ne; i++) {
      float sxx = VGetPixel(sI, 0, 0, i);
      float syy = VGetPixel(sI, 1, 0, i);
      float szz = VGetPixel(sI, 2, 0, i);
      float sxy = VGetPixel(sI, 3, 0, i);
      float sxz = VGetPixel(sI, 4, 0, i);
      float syz = VGetPixel(sI, 5, 0, i);
      fprintf(outf," %d %15.7le %15.7le %15.7le %15.7le %15.7le %15.7le\n",
	      i,sxx,syy,szz,sxy,sxz,syz);
    }
  }

  if (dI != NULL) {
    fprintf(stdout," %d displacements ...\n",nv);
    fflush(stdout);

    fprintf(outf,"displacement\n");
    for (int i=0; i<nv; i++) {
      float dx = VGetPixel(dI, 0, 0, i);
      float dy = VGetPixel(dI, 1, 0, i);
      float dz = VGetPixel(dI, 2, 0, i);
      fprintf(outf," %d %15.7le %15.7le %15.7le\n", i,dx,dy,dz);
    }
  }

  fprintf(stdout," ... done vtoa \n");
  fflush(stdout);

  return;
}






void mesh::vtodx( FILE *outf )
  // vista to DX format conversion
{
  VLong 		nv=0, ne=0;
  VUByte          npart=1;

  fprintf(stdout," entering vtodx ... \n");
  fflush(stdout);

  assert(vertices   != NULL);
  assert(primitives != NULL);

  if (VGetAttr(VGraphAttrList(vertices), "nverts", 0, VLongRepn, &nv) != VAttrFound) {
    nv = countNodes(vertices);
    if (nv <= 0) {
      fprintf(stderr,"vtoa: ERROR number of vertices undefined! \n");
      return;
    }
  }
  if (VGetAttr(VGraphAttrList(vertices), "partition(s)", 0, VUByteRepn, &npart) != VAttrFound) {
    npart = 1;
  }
  if (VGetAttr(VGraphAttrList(primitives), "nelems", 0, VLongRepn, &ne) != VAttrFound) {
    ne = countNodes(primitives);
    if (ne <= 0) {
      fprintf(stderr,"vtoa: ERROR number of primitives undefined! \n");
      return;
    }
  }


  //attribute image pointers
  assert(VGetAttr(VGraphAttrList(primitives), "partelem", 0, VImageRepn, &pe)==VAttrFound);
  assert(VGetAttr(VGraphAttrList(vertices),   "partnode", 0, VImageRepn, &pn)==VAttrFound);
  assert(VGetAttr(VGraphAttrList(primitives), "matprops", 0, VImageRepn, &mp)==VAttrFound);

  if (VGetAttr(VGraphAttrList(primitives), "stressimage", 0, VImageRepn, &sI)!=VAttrFound) {
    sI = NULL;
  }
  if (VGetAttr(VGraphAttrList(vertices),   "srcnodes",    0, VImageRepn, &iv)!=VAttrFound) {
    iv = NULL;
  }
  if (VGetAttr(VGraphAttrList(vertices),   "image",       0, VImageRepn, &tI)!=VAttrFound) {
    tI = NULL;
  }
  if (VGetAttr(VGraphAttrList(vertices),   "forceimage",  0, VImageRepn, &fI)!=VAttrFound) {
    fI = NULL;
  }
  if (VGetAttr(VGraphAttrList(vertices),   "deformation", 0, VImageRepn, &dI)!=VAttrFound) {
    dI = NULL;
  }

  if ( mp == NULL || pn == NULL || pe == NULL )  {
    fprintf(stderr, "vtodx: ERROR found NULL image pointers.\n");
    return;
  }

  fprintf(stdout," %d nodes ... \n",nv);
  fflush(stdout);

  int objcnt = 1;
  //write nodes
  fprintf(outf, "Object %d class array type float rank 1 shape 3 items %d data follows\n", objcnt, nv);

  for (int i=0; i<nv; i++) {
    vertex *v = (vertex *)VGraphGetNode(vertices, i+1);
    if (v != NULL) {
      fprintf(outf,"%15.7le %15.7le %15.7le\n",v->x, v->y, v->z);
    }
  }
  fprintf(outf, "attribute \"dep\" string \"positions\"\n");
  ++objcnt;

  fprintf(stdout," %d elements ... \n",ne);
  fflush(stdout);

  // assuming the same element type througout
  int nodes_per_element = ((primitive *)VGraphGetNode(primitives, 1))->vcnt;
  std::string element_type(nodes_per_element == 4 ? std::string("tetrahedra") : std::string("cubes"));
  fprintf(outf, "Object %d class array type int rank 1 shape %d  items %d  data follows\n", 
	  objcnt, nodes_per_element, ne);
  //write elements
  for (int i=0; i<ne; i++) {
    primitive *p = (primitive *)VGraphGetNode(primitives, i+1);
    if (p != NULL) {
      assert(p->vcnt == nodes_per_element);
      if(nodes_per_element == 4)  // tet
	for (int j=0; j<p->vcnt; j++) {
	  fprintf(outf,"%i ", p->id[j] - 1);
	}
      else { // cube
	int * v = p->id;
	fprintf(outf, "%d %d %d %d %d %d %d %d", v[0]-1, v[1]-1, v[3]-1, v[2]-1, v[4]-1, v[5]-1, v[7]-1, v[6]-1);
      }
      fprintf(outf,"\n");
    }
  }
  fprintf(outf, "attribute \"element type\" string \"%s\"\n", element_type.c_str());
  fprintf(outf, "attribute \"ref\" string \"positions\"\n");
  ++objcnt;

  // boundary condition tags
  if(tI != NULL) {
    fprintf(outf, "Object \"boundary_conditions\" class array  type float rank 0   items %d data follows\n", 
	    nv);
    for (int i=0; i<nv; i++) {
      fprintf(outf, "%d ", (int)VGetPixel(tI, 0, 0, i));
    }
    fprintf(outf, "\nattribute \"dep\" string \"positions\"\n");
  }

  // prescribed forces or displacements 
  if(fI != NULL) {
    fprintf(outf, "Object \"forces\" class array  type float rank 1 shape 3 items %d data follows\n", 
	    nv);
    for (int i=0; i<nv; i++) {
      fprintf(outf, "%f %f %f\n", 
	      (float)VGetPixel(fI,0,0,i), VGetPixel(fI,1,0,i), VGetPixel(fI,2,0,i));
    }  
    fprintf(outf, "\nattribute \"dep\" string \"positions\"\n");
  }

  // materials
  if(mp != NULL) {
    fprintf(outf, "Object \"materials\" class array type int rank 0 item %d  data follows\n", ne);
    for(int i=0; i < ne; ++i) {
      fprintf(outf, "%d ", (int)VGetPixel(mp, 0,0,i));
    }
    fprintf(outf, "\nattribute \"dep\" string \"connections\"\n");
  }

  // cell partitions
  if(pe != NULL) {
       fprintf(outf, "Object \"cell_partitions\" class array type int rank 0 item %d  data follows\n", ne);
    for(int i=0; i < ne; ++i) {
      fprintf(outf, "%d ", (int)VGetPixel(pe, 0,0,i));
    }
    fprintf(outf, "\nattribute \"dep\" string \"connections\"\n");
   }


  fprintf(outf, "Object \"mesh\" class field\n");
  fprintf(outf, "component \"positions\"   value 1\n");
  fprintf(outf, "component \"connections\" value 2\n");
  if(tI != NULL) fprintf(outf,"component \"boundary_conditions\" value \"boundary_conditions\"\n");
  if(fI != NULL) fprintf(outf,"component \"forces\" value \"forces\"\n");
  if(mp != NULL) fprintf(outf,"component \"materials\" value \"materials\"\n");
  if(fI != NULL) fprintf(outf,"component \"cell_partitions\" value \"cell_partitions\"\n");

  fprintf(outf,"end\n");

  return;
}








void mesh::vtogmv( FILE *outf, bool do_scale)
  // vista to gmv format conversion
{
  VLong 		nv=0, ne=0;
  VUByte          npart=1;
  VImage          eI;

  if(olevel_ > 0) {
    fprintf(stdout," entering vtogmv ... \n");
    fflush(stdout);
  }

  assert(vertices   != NULL);
  assert(primitives != NULL);

  if (VGetAttr(VGraphAttrList(vertices), "nverts", 0, VLongRepn, &nv) != VAttrFound) {
    nv = countNodes(vertices);
    if (nv <= 0) {
      fprintf(stderr,"vtogmv: ERROR number of vertices undefined! \n");
      return;
    }
  }
  if (VGetAttr(VGraphAttrList(vertices), "partition(s)", 0, VUByteRepn, &npart) != VAttrFound) {
    npart = 1;
  }
  if (VGetAttr(VGraphAttrList(primitives), "nelems", 0, VLongRepn, &ne) != VAttrFound) {
    ne = countNodes(primitives);
    if (ne <= 0) {
      fprintf(stderr,"vtogmv: ERROR number of primitives undefined! \n");
      return;
    }
  }
  // get bounding box
  double scale = 1, transx=0, transy=0, transz=0;
  if(do_scale) {
    vertex* v = (vertex *)VGraphGetNode(vertices, 1);
    double xmax=v->x, xmin=v->x, ymax=v->y, ymin=v->y, zmax=v->z, zmin=v->z;
    for (int i=0; i<nv; i++) {
      vertex *v = (vertex *)VGraphGetNode(vertices, i+1);
      assert(v != NULL);
      xmax = ( v->x > xmax ? v->x : xmax);
      ymax = ( v->y > ymax ? v->y : ymax);
      zmax = ( v->z > zmax ? v->z : zmax);
      xmin = ( v->x < xmin ? v->x : xmin);
      ymin = ( v->y < ymin ? v->y : ymin);
      zmin = ( v->z < zmin ? v->z : zmin);
	    
    }
    double max_bb_edge = max(xmax-xmin,max(ymax-ymin, zmax-zmin));
    scale = 2.0/max_bb_edge;
    transx = -0.5*(xmin + xmax);
    transy = -0.5*(ymin + ymax);
    transz = -0.5*(zmin + zmax);
  }

	

  fprintf(outf,"gmvinput ascii \n");
  fprintf(outf,"comments \n");
  fprintf(outf," version date: 19.8.2002 \n");
  fprintf(outf,"endcomm \n");

  //attribute images
  assert(VGetAttr(VGraphAttrList(primitives), "matprops", 0, VImageRepn, &mp)==VAttrFound);

  //assert(VGetAttr(VGraphAttrList(primitives), "partelem", 0, VImageRepn, &pe)==VAttrFound);
  //assert(VGetAttr(VGraphAttrList(vertices),   "partnode", 0, VImageRepn, &pn)==VAttrFound);

  if (VGetAttr(VGraphAttrList(primitives), "partelem", 0, VImageRepn, &pe)!=VAttrFound) {
    pe = NULL;
  }
  if (VGetAttr(VGraphAttrList(vertices),   "partnode", 0, VImageRepn, &pn)!=VAttrFound) {
    pn = NULL;
  }

  if (VGetAttr(VGraphAttrList(primitives), "stressimage", 0, VImageRepn, &sI)!=VAttrFound) {
    sI = NULL;
  }
  if (VGetAttr(VGraphAttrList(vertices),   "srcnodes",    0, VImageRepn, &iv)!=VAttrFound) {
    iv = NULL;
  }
  if (VGetAttr(VGraphAttrList(vertices),   "image",       0, VImageRepn, &tI)!=VAttrFound) {
    tI = NULL;
  }
  if (VGetAttr(VGraphAttrList(vertices),   "forceimage",  0, VImageRepn, &fI)!=VAttrFound) {
    fI = NULL;
  }
  if (VGetAttr(VGraphAttrList(vertices),   "deformation", 0, VImageRepn, &dI)!=VAttrFound) {
    dI = NULL;
  }
  if (VGetAttr(VGraphAttrList(vertices),   "error", 0, VImageRepn, &eI)!=VAttrFound) {
    eI = NULL;
  }

  if ( mp == NULL ) { // || pn == NULL || pe == NULL )  {
    fprintf(stderr, "vtogmv: ERROR found NULL image pointers.\n");
    return;
  }

  if(olevel_ > 0) {
    fprintf(stdout," %d nodes ...\n",nv);
    fflush(stdout);
  }

  //write nodes
  fprintf(outf,"nodes %i\n",nv);

  for (int i=0; i<nv; i++) {
    vertex *v = (vertex *)VGraphGetNode(vertices, i+1);
    assert(v != NULL);
    if(!do_scale)
      fprintf(outf," %g",v->x);
    else
      fprintf(outf," %g",(v->x +transx)*scale);
  }
  fprintf(outf,"\n");

  for (int i=0; i<nv; i++) {
    vertex *v = (vertex *)VGraphGetNode(vertices, i+1);
    assert(v != NULL);
    if(!do_scale)
      fprintf(outf," %g",v->y);
    else
      fprintf(outf," %g",(v->y +transy)*scale);
  }
  fprintf(outf,"\n");

  for (int i=0; i<nv; i++) {
    vertex *v = (vertex *)VGraphGetNode(vertices, i+1);
    assert(v != NULL);

    if(!do_scale)
      fprintf(outf," %g",v->z);
    else
      fprintf(outf," %g",(v->z +transy)*scale);

  }
  fprintf(outf,"\n");


  if(olevel_ > 0) {
    fprintf(stdout," %d elements ...\n",ne);
    fflush(stdout);
  }

  //write elements
  fprintf(outf,"cells %i\n",ne);

  int matmin=(int)VGetPixel(mp, 0, 0, 0);
  int matmax=matmin;
  for (int i=0; i<ne; i++) {
    primitive *p = (primitive *)VGraphGetNode(primitives, i+1);
    assert(p != NULL);

    int mid = (int)VGetPixel(mp, 0, 0, i);
    if (mid < matmin) matmin=mid;
    if (mid > matmax) matmax=mid;

    if (p->vcnt == 8) fprintf(outf,"hex 8 \n");
    if (p->vcnt == 4) fprintf(outf,"tet 4 \n");

    for (int j=0; j<p->vcnt; j++) {
      fprintf(outf," %i", p->id[j]);   //    1 ... nv
    }
    fprintf(outf,"\n");
  }

  //write element labels
  if(olevel_ > 0) {
    fprintf(stdout," element labels in [%d %d] ...\n",matmin,matmax);
    fflush(stdout);
  }

  //                                              0-cell data
  fprintf(outf,"material %i %i\n",matmax-matmin+1,0);

  for (int i=matmin; i<=matmax; i++) {
    fprintf(outf,"mat%i \n", i-matmin+1);
  }
  for (int i=0; i<ne; i++) {
    int mid = (int)VGetPixel(mp, 0, 0, i);
    fprintf(outf," %i", mid-matmin+1);
  }
  fprintf(outf,"\n");

  //field variables: 0-cell data, 1-node-data, 2-face data
  fprintf(outf,"variable\n");

  if(pe != NULL) {
    fprintf(outf,"parte 0\n");
    for (int i=0; i<ne; i++) {
      int pid = (int)VGetPixel(pe, 0, 0, i);
      fprintf(outf," %g", (float)pid );
    }
    fprintf(outf,"\n");
  }
  if(pn != NULL) {
    fprintf(outf,"partn 1\n");
    for (int i=0; i<nv; i++) {
      int pid = (int)VGetPixel(pn, 0, 0, i);
      fprintf(outf," %g", (float)pid );
    }
    fprintf(outf,"\n");
  }
  if (dI != NULL) {
    fprintf(stdout," displacements ...\n");
    fflush(stdout);

    fprintf(outf,"displx 1\n");
    for (int i=0; i<nv; i++) {
      float dx = VGetPixel(dI, 0, 0, i);
      fprintf(outf," %g",dx);
    }
    fprintf(outf,"\n");

    fprintf(outf,"disply 1\n");
    for (int i=0; i<nv; i++) {
      float dy = VGetPixel(dI, 1, 0, i);
      fprintf(outf," %g",dy);
    }
    fprintf(outf,"\n");

    fprintf(outf,"displz 1\n");
    for (int i=0; i<nv; i++) {
      float dz = VGetPixel(dI, 2, 0, i);
      fprintf(outf," %g",dz);
    }
    fprintf(outf,"\n");

    fprintf(outf,"displ 1\n");
    for (int i=0; i<nv; i++) {
      float dx = VGetPixel(dI, 0, 0, i);
      float dy = VGetPixel(dI, 1, 0, i);
      float dz = VGetPixel(dI, 2, 0, i);
      fprintf(outf," %g",sqrt(dx*dx+dy*dy+dz*dz));
    }
    fprintf(outf,"\n");
  }

  if (fI != NULL) {

    if(olevel_ > 0) {
      fprintf(stdout," forces ...\n");
      fflush(stdout);
    }
    fprintf(outf,"forcex 1\n");
    for (int i=0; i<nv; i++) {
      float fx = VGetPixel(fI, 0, 0, i);
      fprintf(outf," %g",fx);
    }
    fprintf(outf,"\n");

    fprintf(outf,"forcey 1\n");
    for (int i=0; i<nv; i++) {
      float fy = VGetPixel(fI, 1, 0, i);
      fprintf(outf," %g",fy);
    }
    fprintf(outf,"\n");

    fprintf(outf,"forcez 1\n");
    for (int i=0; i<nv; i++) {
      float fz = VGetPixel(fI, 2, 0, i);
      fprintf(outf," %g",fz);
    }
    fprintf(outf,"\n");
  }

  if (eI != NULL) {
    fprintf(stdout," errors ...\n");
    fflush(stdout);

    fprintf(outf,"errorx 1\n");
    for (int i=0; i<nv; i++) {
      float dx = VGetPixel(eI, 0, 0, i);
      fprintf(outf," %g",dx);
    }
    fprintf(outf,"\n");

    fprintf(outf,"errory 1\n");
    for (int i=0; i<nv; i++) {
      float dy = VGetPixel(eI, 1, 0, i);
      fprintf(outf," %g",dy);
    }
    fprintf(outf,"\n");

    fprintf(outf,"errorz 1\n");
    for (int i=0; i<nv; i++) {
      float dz = VGetPixel(eI, 2, 0, i);
      fprintf(outf," %g",dz);
    }
    fprintf(outf,"\n");

    fprintf(outf,"error 1\n");
    for (int i=0; i<nv; i++) {
      float dx = VGetPixel(eI, 0, 0, i);
      float dy = VGetPixel(eI, 1, 0, i);
      float dz = VGetPixel(eI, 2, 0, i);
      fprintf(outf," %g",sqrt(dx*dx+dy*dy+dz*dz));
    }
    fprintf(outf,"\n");
  }



  if (sI != NULL) {
    if(olevel_ > 0) {
      fprintf(stdout," stresses ...\n");
      fflush(stdout);
    }

    fprintf(outf,"stressxx 0\n");
    for (int i=0; i<ne; i++) {
      float s = VGetPixel(sI, 0, 0, i);
      fprintf(outf," %g",s);
    }
    fprintf(outf,"\n");

    fprintf(outf,"stressyy 0\n");
    for (int i=0; i<ne; i++) {
      float s = VGetPixel(sI, 1, 0, i);
      fprintf(outf," %g",s);
    }
    fprintf(outf,"\n");

    fprintf(outf,"stresszz 0\n");
    for (int i=0; i<ne; i++) {
      float s = VGetPixel(sI, 2, 0, i);
      fprintf(outf," %g",s);
    }
    fprintf(outf,"\n");

    fprintf(outf,"stressxy 0\n");
    for (int i=0; i<ne; i++) {
      float s = VGetPixel(sI, 3, 0, i);
      fprintf(outf," %g",s);
    }
    fprintf(outf,"\n");

    fprintf(outf,"stressxz 0\n");
    for (int i=0; i<ne; i++) {
      float s = VGetPixel(sI, 4, 0, i);
      fprintf(outf," %g",s);
    }
    fprintf(outf,"\n");

    fprintf(outf,"stressyz 0\n");
    for (int i=0; i<ne; i++) {
      float s = VGetPixel(sI, 5, 0, i);
      fprintf(outf," %g",s);
    }
    fprintf(outf,"\n");
  }

  fprintf(outf,"endvars\n");


  //write flags
  if(tI != 0 || iv != 0) {
    if(olevel_ > 0) {
      fprintf(stdout," flags ...\n");
      fflush(stdout);
    }

    fprintf(outf,"flags\n");

    if(tI != 0) {
      fprintf(outf,"tagtype 8 1\n");
      fprintf(outf,"free xfix yfix xyfix zfix xzfix yzfix xyzfix\n");
	    
      for (int i=0; i<nv; i++) {
	int fid;
	if (tI != NULL) {
	  fid = (int)VGetPixel(tI, 0, 0, i);
	}
	else {
	  fid = 0;
	}
	fprintf(outf," %d",fid+1);
      }
      fprintf(outf,"\n");
    }
    if( iv != 0) {
      fprintf(outf,"srctype 2 1\n");
      fprintf(outf,"off on\n");
	    
      for (int i=0; i<nv; i++) {
	int sid;
	if (iv != NULL) {
	  sid = (int)VGetPixel(iv, 0, 0, i);
	}
	else {
	  sid = 0;
	}
	fprintf(outf," %d",sid+1);
      }
      fprintf(outf,"\n");
    }
    fprintf(outf,"endflag\n");
  }

  fprintf(outf, "endgmv\n");

  if(olevel_ > 0) {
    fprintf(stdout," ... done vtogmv \n");
    fflush(stdout);
  }
  return;
}



void mesh::vtogpp( FILE *outf )
  // vista to gpp format conversion
{
  VLong 		nv=0, ne=0;
  VUByte          npart=1;

  fprintf(stdout," entering vtogpp ... \n");
  fflush(stdout);

  assert(vertices   != NULL);
  assert(primitives != NULL);

  if (VGetAttr(VGraphAttrList(vertices), "nverts", 0, VLongRepn, &nv) != VAttrFound) {
    nv = countNodes(vertices);
    if (nv <= 0) {
      fprintf(stderr,"vtogpp: ERROR number of vertices undefined! \n");
      return;
    }
  }
  if (VGetAttr(VGraphAttrList(vertices), "partition(s)", 0, VUByteRepn, &npart) != VAttrFound) {
    npart = 1;
  }
  if (VGetAttr(VGraphAttrList(primitives), "nelems", 0, VLongRepn, &ne) != VAttrFound) {
    ne = countNodes(primitives);
    if (ne <= 0) {
      fprintf(stderr,"vtogpp: ERROR number of primitives undefined! \n");
      return;
    }
  }

  //attribute image pointers
  assert(VGetAttr(VGraphAttrList(primitives), "partelem",   0, VImageRepn, &pe)==VAttrFound);
  assert(VGetAttr(VGraphAttrList(vertices),   "partnode",   0, VImageRepn, &pn)==VAttrFound);
  assert(VGetAttr(VGraphAttrList(primitives), "matprops",   0, VImageRepn, &mp)==VAttrFound);

  if (VGetAttr(VGraphAttrList(primitives), "stressimage", 0, VImageRepn, &sI)!=VAttrFound) {
    sI = NULL;
  }
  if (VGetAttr(VGraphAttrList(vertices),   "srcnodes",    0, VImageRepn, &iv)!=VAttrFound) {
    iv = NULL;
  }
  if (VGetAttr(VGraphAttrList(vertices),   "image",       0, VImageRepn, &tI)!=VAttrFound) {
    tI = NULL;
  }
  if (VGetAttr(VGraphAttrList(vertices),   "forceimage",  0, VImageRepn, &fI)!=VAttrFound) {
    fI = NULL;
  }

  if ( mp == NULL || pn == NULL || pe == NULL )  {
    fprintf(stderr, "vtogpp: ERROR found NULL image pointers.\n");
    return;
  }

  int matmin=(int)VGetPixel(mp, 0, 0, 0);
  int matmax=matmin;
  int ntet=0, nhex=0;
  for (int i=0; i<ne; i++) {
    int mid = (int)VGetPixel(mp, 0, 0, i);
    if (mid < matmin) matmin=mid;
    if (mid > matmax) matmax=mid;
    if (mid < 100) {
      ntet++;
    }
    else {
      nhex++;
    }
  }
  assert( ne == ntet+nhex );

  fprintf(outf,"%d %d %d %d %d\n",nv,0,0,ntet,nhex);

  fprintf(stdout," %d nodes ... \n",nv);
  fflush(stdout);

  //write nodes
  for (int i=0; i<nv; i++) {
    vertex *v = (vertex *)VGraphGetNode(vertices, i+1);
    assert( v != NULL );

    fprintf(outf," %e %e %e\n",v->x,v->y,v->z);
  }

  //
  // elements ordered by type: 1. tetra   2. hexa
  //
  fprintf(stdout," %d elements ... \n",ne);
  fflush(stdout);

  //write tetrahedral elements
  if (ntet > 0) {
    for (int i=0; i<ne; i++) {
      primitive *p = (primitive *)VGraphGetNode(primitives, i+1);
      assert( p != NULL );

      int mid = (int)VGetPixel(mp, 0, 0, i);
  
      if (mid < 100) {
	assert( p->vcnt == 4 );
	for (int j=0; j<4; j++) {
	  fprintf(outf," %i", p->id[j]-1);   //    id in [1 nv]
	  assert( p->id[j] >   0 );
	  assert( p->id[j] <= nv );
	}
	fprintf(outf,"\n");
      }
    }
  }

  //write hexahedral  elements
  if (nhex > 0) {
    for (int i=0; i<ne; i++) {
      primitive *p = (primitive *)VGraphGetNode(primitives, i+1);
      assert( p != NULL );

      int mid = (int)VGetPixel(mp, 0, 0, i);
   
      if (mid >= 100) {
	assert( p->vcnt == 8 );
	for (int j=0; j<8; j++) {
	  fprintf(outf,"%i ", p->id[j]-1);   //    id in [1 nv]
	  assert( p->id[j] >   0 );
	  assert( p->id[j] <= nv );
	}
	fprintf(outf,"\n");
      }
    }
  }

  fprintf(stdout," ... done vtogpp \n");
  fflush(stdout);

  return;
}



void mesh::vista2drama(int *dim,int *dist,double* coord,int* enptr,int* con,
                       int* gnn,double* force,int* constr,VAttrList *list,
                       int *myid)
{
  VLong 		nv=0, ne=0;
  VUByte		npart=1;

  assert(*myid      >= 0   );
  assert(dim        != NULL);
  assert(dist       != NULL);
  assert(coord      != NULL);
  assert(enptr      != NULL);
  assert(con        != NULL);
  assert(gnn        != NULL);
  assert(force      != NULL);
  assert(constr     != NULL);
  assert(vertices   != NULL);
  assert(primitives != NULL);

  // copy vertices attribute list to list
  *list = VCopyAttrList(VGraphAttrList(vertices));

  if (VGetAttr(VGraphAttrList(vertices), "nverts", 0, VLongRepn, &nv) != VAttrFound) {
    nv = countNodes(vertices);
    if (nv <= 0) {
      fprintf(stderr,"vista2drama: ERROR number of vertices undefined! \n");
      return;
    }
  }
  if (VGetAttr(VGraphAttrList(vertices), "partition(s)", 0, VUByteRepn, &npart) != VAttrFound) {
    npart = 1;
  }
  if (VGetAttr(VGraphAttrList(primitives), "nelems", 0, VLongRepn, &ne) != VAttrFound) {
    ne = countNodes(primitives);
    if (ne <= 0) {
      fprintf(stderr,"vista2drama: ERROR number of primitives undefined! \n");
      return;
    }
  }

  //attribute image pointers
  assert(VGetAttr(VGraphAttrList(primitives), "partelem",   0, VImageRepn, &pe)==VAttrFound);
  assert(VGetAttr(VGraphAttrList(vertices),   "partnode",   0, VImageRepn, &pn)==VAttrFound);
  assert(VGetAttr(VGraphAttrList(primitives), "matprops",   0, VImageRepn, &mp)==VAttrFound);

  if (VGetAttr(VGraphAttrList(primitives), "stressimage", 0, VImageRepn, &sI)!=VAttrFound) {
    sI = NULL;
  }
  if (VGetAttr(VGraphAttrList(vertices),   "srcnodes",    0, VImageRepn, &iv)!=VAttrFound) {
    iv = NULL;
  }
  if (VGetAttr(VGraphAttrList(vertices),   "image",       0, VImageRepn, &tI)!=VAttrFound) {
    tI = NULL;
  }
  if (VGetAttr(VGraphAttrList(vertices),   "forceimage",  0, VImageRepn, &fI)!=VAttrFound) {
    fI = NULL;
  }

  if ( mp == NULL || pn == NULL || pe == NULL )  {
    fprintf(stderr, "vista2drama: ERROR found NULL image pointers.\n");
    return;
  }

  // nodes
  int count=0;
  for(int j=0; j<nv; j++){
    int pid=(int)VGetPixel(pn, 0, 0, j);
    if (pid == (*myid)) {
      int id;
      if (iv != NULL) {	
	id = (int)VGetPixel(iv, 0, 0, j);
      }
      else {
	id = 0;
      }
      gnn[count]=id;
			
      vertex *v = (vertex *)VGraphGetNode(vertices, j+1);
      assert(v != NULL);
			
      coord[3*count    ] = (double)v->x;
      coord[3*count + 1] = (double)v->y;
      coord[3*count + 2] = (double)v->z;
	
      if (fI != NULL) {		
	force[3*count    ] = (double)VGetPixel(fI, 0, 0, j);
	force[3*count + 1] = (double)VGetPixel(fI, 1, 0, j);
	force[3*count + 2] = (double)VGetPixel(fI, 2, 0, j);
      }
      else {
	force[3*count    ] = 0.0;
	force[3*count + 1] = 0.0;
	force[3*count + 2] = 0.0;
      }

      if (tI != NULL) {
	constr[count] = (int)VGetPixel(tI, 0, 0, j);
      }
      else {
	constr[count] = 0;
      }
			
      count++;		
			
    }	
  }

  // fill enptr field
  int ecount =0;
  int lastptr=1;
  int min=dim[6];
  if (ne > 0) {	
    for (int j=0; j<ne; j++) {
      int ped=(int)VGetPixel(pe, 0, 0, j);
      if (ped == (*myid)){
	primitive *p = (primitive *)VGraphGetNode(primitives, j+1);
	assert( p != NULL );

	int mid=(int)VGetPixel(mp, 0, 0, j);
	enptr[2*ecount    ] = lastptr;
	enptr[2*ecount + 1] = mid-min+1;
				
	lastptr+=p->vcnt;
	ecount++;	
				
      }
    }
    if (dim[0] != ecount) {
      fprintf(stderr,"P%d mesh::vista2drama WARNING inconsistent number of elements!\n",*myid);
    }
		
    enptr[2*dim[0]    ] = lastptr;
    enptr[2*dim[0] + 1] = -1;
  }
  else{
    enptr[2*dim[0]    ] = lastptr;
    enptr[2*dim[0] + 1] = -1;
  }
	
  // fill con field
  ecount = 0;
  for(int j=0; j<ne; j++){
    int pid=(int)VGetPixel(pe, 0, 0, j);
    if (pid == (*myid)){

      primitive *p = (primitive *)VGraphGetNode(primitives, j+1);
      assert(p != NULL);

      for(int i=0; i<p->vcnt; i++){
	int gid = p->id[i];
	int processor=0;
	while (gid > dist[processor]){
	  processor++;
	}
					
	con[2*ecount    ] = gid - dist[processor-1];
	con[2*ecount + 1] = processor-1;
	ecount++;
      }	
    }
  }
	
  return;
}



void mesh::colorize(VImage src)
{
  resizeForColor();

  VStringConst attr;
  vec3 dim = 1;
  if (VGetAttr(VImageAttrList(src), "voxel", 0, VStringRepn, &attr) == VAttrFound)
    sscanf(attr, "%lf%lf%lf", &dim.x, &dim.y, &dim.z);

  for (vertex *v = (vertex *)VGraphFirstNode(vertices); v;
       v = (vertex *)VGraphNextNode(vertices))  {
    vec3 p = v->getPoint()/dim;
    v->setColor(getPixel(src, p));
  };
}


void mesh::add_defo(VImage defo)
{
  VStringConst attr;
  vec3 dim = 1;
  if (VGetAttr(VImageAttrList(defo), "voxel", 0, VStringRepn, &attr) == VAttrFound)
    sscanf(attr, "%lf%lf%lf", &dim.x, &dim.y, &dim.z);

  VImage deformation = VCreateImage(3, 1, countNodes(vertices),VFloatRepn);
  VImageNFrames(deformation) = 1;
  VImageNComponents(deformation) = 3;
		
  int i = 0; 
	
  for (vertex *v = (vertex *)VGraphFirstNode(vertices); v;
       v = (vertex *)VGraphNextNode(vertices))  {
    vec3 p = v->getPoint();
    vec3 pv = p/dim;
		
    // this is rape
    vec3 f = getForce(defo, pv);
    VPixel(deformation, 0, 0, i, VFloat) = f.x; 
    VPixel(deformation, 1, 0, i, VFloat) = f.y; 
    VPixel(deformation, 2, 0, i, VFloat) = f.z; 		

    ++i; 
  }
	
  VSetAttr(VImageAttrList(deformation),"component_repn", 0, VStringRepn, "vector3");
  VSetAttr(VImageAttrList(deformation),"component_interp",0, VStringRepn, "deformation"); 
	
  VSetAttr(VGraphAttrList(vertices), "deformation", 0, VImageRepn, deformation);
}




void mesh::compact()
  // compact graphs
{
  VGraph nvertices = NULL, nfaces = NULL;
  map< VLong , VLong , less<VLong> > vertexmap, facemap;

  // remove all faces with illegal vertices
  if (faces) {
    for (int i = 1; i <= faces->lastUsed; i++)  {
      VNode n = VGraphGetNode(faces, i);
      if (n == NULL) continue;
      int vcnt = ((VLong*)n->data)[0];
      for (int j = 1; j <= vcnt; j++)  {
	if (!VGraphGetNode(vertices, ((VLong*)n->data)[j]))  {
	  VDestroyNode(faces, i);
	  break;
	};
      }
    }
  }

  // allocate new graphs
  if (vertices) {
    unsigned int n = countNodes(vertices);
    if (n == 0)  {
      fprintf(stderr, "mesh does not contain any vertices.\n");
      return;
    };
    nvertices = VCreateGraph(n, vertices->nfields, vertices->node_repn, vertices->useWeights);
    if (nvertices == 0)  {
      fprintf(stderr, "cannot create new vertex graph.\n");
      return;
    };
    if (VGraphAttrList(vertices))
      VGraphAttrList(nvertices) = VCopyAttrList(VGraphAttrList(vertices));
  }
  if (faces) {
    unsigned int n = countNodes(faces);
    if (n == 0)  {
      fprintf(stderr, "mesh does not contain any faces.\n");
      return;
    };
    nfaces = VCreateGraph(n, faces->nfields, faces->node_repn, faces->useWeights);
    if (nfaces == 0)  {
      fprintf(stderr, "cannot create new face graph.\n");
      return;
    };
    if (VGraphAttrList(faces))
      VGraphAttrList(nfaces) = VCopyAttrList(VGraphAttrList(faces));
  }

  // copy vertices and faces
  if (vertices) {
    int j = 1;
    for (int i = 1; i <= vertices->lastUsed; i++)  {
      VNode n = VGraphGetNode(vertices, i);
      if (n == NULL) continue;
      VGraphAddNodeAt(nvertices, n, j);
      vertexmap[i] = j; j++;
    };
  };
  if (faces) {
    int j = 1;
    for (int i = 1; i <= faces->lastUsed; i++)  {
      VNode n = VGraphGetNode(faces, i);
      if (n == NULL) continue;
      VGraphAddNodeAt(nfaces, n, j);
      facemap[i] = j; j++;
    };
  };

  // add neighbours of vertices and faces
  if (vertices) {
    for (int i = 1; i <= vertices->lastUsed; i++)  {
      VNode n = VGraphGetNode(vertices, i);
      if (n == NULL) continue;
      for (VAdjacency adj = n->base.head; adj; adj = adj->next)
	VGraphLinkNodes(nvertices, vertexmap[i], vertexmap[adj->id]);
    }
  };
  if (faces) {
    for (int i = 1; i <= faces->lastUsed; i++)  {
      VNode n = VGraphGetNode(faces, i);
      if (n == NULL) continue;
      for (VAdjacency adj = n->base.head; adj; adj = adj->next)
	VGraphLinkNodes(nfaces, facemap[i], facemap[adj->id]);

      // update indices in face
      VNode nn = VGraphGetNode(nfaces, facemap[i]);
      VLong vcnt = ((VLong *)n->data)[0];
      for (int j = 1; j <= vcnt; j++)
	((VLong *)nn->data)[j] = vertexmap[((VLong *)n->data)[j]];
    };
  };

  // cleanup
  if (vertices)  { VDestroyGraph(vertices); vertices = nvertices; };
  if (faces)     { VDestroyGraph(faces); faces = nfaces; };
  vertexmap.erase(vertexmap.begin(), vertexmap.end());
  facemap.erase(facemap.begin(), facemap.end());
}  


