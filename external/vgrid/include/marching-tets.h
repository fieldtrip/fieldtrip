#ifndef SIMBIO_MARCHING_TETS_H
#define SIMBIO_MARCHING_TETS_H

/*
  © Copyright 2003, C&C Research Laboratories, NEC Europe Ltd.
  

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

/*! \file 

    Support for the generalized volumetric marching tetrahedra algorithm

    \author Guntram Berti <berti@ccrl-nece.de>

     $Id$

*/


#include "pre-post-conditions.h"


#include <map>
#include <vector>
#include <algorithm> // sort, min/max


/*! Contains classes for implementing the volumetric marching tets algorithm

    \author Guntram Berti (berti@ccrl-nece.de)
*/
namespace marching_tets {
 
  class edge {
    typedef edge self;
    int v[2];
  public:
    edge() {}
    edge(int i, int j) { if (i<j) { v[0]=i;v[1]=j;} else { v[0]=j; v[1]=i;} }
    int operator[](int i)  const { REQUIRE( 0 <= i && i <= 1, "i = " << i, 1);  return v[i];} 

    bool operator< (self const& rhs) const 
    { return v[0] < rhs.v[0] || (v[0] == rhs.v[0] && v[1] < rhs.v[1]);}
    bool operator==(self const& rhs) const 
    { return rhs.v[0] == v[0] && rhs.v[1]==v[1];}
    bool operator!=(self const& rhs) const { return !(*this == rhs);}
 };

  struct hash_edge {
    typedef edge key_type;
    typedef key_type argument_type;
    typedef size_t result_type;
    result_type operator()(key_type const& e) const { return 10*e[1]+e[0];}
  };
}


namespace marching_tets {

  class triangle {
    typedef triangle self;
    int  v[3];
  public:
    triangle() {}
    triangle(int t0, int t1, int t2) 
    { v[0] = t0; v[1] = t1; v[2] = t2; }
    
    int  operator[](int i) const { return v[i];}
    int& operator[](int i)       { return v[i];}

  };

  inline  std::ostream& operator<<(std::ostream& out, triangle const& t) {
    return(out << t[0] << ' ' << t[1] << ' ' << t[2]);
  }

  class tet {
    typedef tet self;
    int  v[4];
  public:
    tet() {}
    tet(int t0, int t1, int t2, int t3) 
    { v[0] = t0; v[1] = t1; v[2] = t2; v[3] = t3;}
    
    int  operator[](int i) const { return v[i];}
    int& operator[](int i)       { return v[i];}

  };

  inline  std::ostream& operator<<(std::ostream& out, tet const& t) {
    return(out << t[0] << ' ' << t[1] << ' ' << t[2] << ' ' << t[3]);
  }

  //typedef STDHASH::hash_map<edge, int, hash_edge> cut_vertex_map;
  typedef std::map<edge,int> cut_vertex_map;

  /*! \brief Subdivision of a tetrahedron according to sign pattern of vertices

     \todo Self-check, for example check whether sign of tets is consistent with signs
           of vertices
  */

  class subdivision 
  {
  private:
    int np; /* number of positive vertices, $n_+$ */
    int n0; /* number of zero     vertices, $n_0$ */
    int nm; /* number of negative vertices, $n_-$ */
    int nc; /* number of cut      vertices  $n_c$ */

    std::vector<int> v;    // ordered vertices of subdivided tet, order:
                      // v^+_1 < ... < v^+_{np},
                      // v^0_1 < ... < v^0_{n0},
                      // v^-_1 < ... < v^-_{nm},
                      // v^c_1 , ... , v^c_{nc}
                      // order on original vertices (v^+, v^0, v^-) is defined by 
                      // global vertex id order,
                      // order on cut vertices v^c does not matter

    int  vp_, v0_, vm_, vc_; // offsets in v
    std::vector<tet>  tets;  // tetrahedral cells subdividing tet t
    std::vector<int>  sign;  // sign of tets, in {+1,-1}

    std::vector<triangle> tris;  // surface triangles
    std::vector<int>      tri_sign;
    cut_vertex_map const* cutv; // cut vertices on edges of t (or the triangulation containing t)

    tet t; // master tet (to be subdivided)
  public:
    subdivision(int             pat[4],      // in {+1,0,-1}^4 : Function value signs on vertices of t
		tet const&      t_,          // the tet to be subdivided
		cut_vertex_map  const& cv);  // the cut vertices on the edges of t
    //! \f$ i \mapsto v^+_i \f$
    int vp(int i) const { REQUIRE((0 <=i && i < np), "i = " << i, 1);  return v[i+vp_];}
    //! \f$ i \mapsto v^-_i \f$
    int vm(int i) const { REQUIRE((0 <=i && i < nm), "i = " << i, 1);  return v[i+vm_];}
    // \f$ i \mapsto v^0_i \f$
    int v0(int i) const { REQUIRE((0 <=i && i < n0), "i = " << i, 1);  return v[i+v0_];}

    //! cut vertex on edge between vertices with global handles i and j
    int vc(int i, int j) const {
      cut_vertex_map::const_iterator e = cutv->find(edge(i,j));
      REQUIRE(e != cutv->end(), " i=" << i << " j=" << j, 1); 
      return e->second;
    }

 
    int NumOfTets()     const { return tets.size();}
    int NumOfCells()    const { return NumOfTets();}
    int NumOfVertices() const { return v.size();}
    int NumOfSurfaceTriangles() const { return tris.size();}

    typedef std::vector<tet>::const_iterator tet_iterator;
    tet_iterator FirstTet() const { return tets.begin();}
    tet_iterator EndTet  () const { return tets.end();}
    int          Sign(tet_iterator const& t) const { return sign[t-tets.begin()];}
    int          Sign(int v) const;

    typedef std::vector<triangle>::const_iterator  triangle_iterator;
    triangle_iterator FirstSurfaceTriangle() const { return tris.begin();}
    triangle_iterator EndSurfaceTriangle  () const { return tris.end  ();}
    // sign = 0 for internal surface triangle, +1 for boundary triangle of positive tet, -1 for neg. tet.
    // can be used to differentiate internal from external surface triangles, 
    // and to choose exactly one of external triangles (occuring from two different tets, one +1, one -1)
    int Sign(triangle_iterator const& t) const { return tri_sign[t-tris.begin()];}

    /*
    typedef vertex_iterator_int<self,0> VertexIterator;
    typedef VertexIterator              Vertex;
    typedef cell_iterator_int  <self,0> CellIterator;
    typedef CellIterator                Cell;
    typedef
    */ 

  private:
    void add_tet(tet const& t, int s)
    { 
      tets.push_back(t);
      sign.push_back(s);
    }
    void add_triangle(triangle const& t, int s) 
    { 
      tris    .push_back(t); 
      tri_sign.push_back(s);
    }
    void create_tets();
  };
}


#endif
