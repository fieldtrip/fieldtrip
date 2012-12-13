#ifndef SIMBIO_GENERATE_GRID_UTIL_H
#define SIMBIO_GENERATE_GRID_UTIL_H

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


/*! \file generateGrid-util.h
    \author Guntram Berti
    
     Some utility types for mesh generation.
     
     $Id$
 */


class edge {
 int v[2];
public:
  edge() { v[0] = v[1] = -1;}
  edge(int v0, int v1) { 
    if(v0 < v1) {
      v[0] = v0; v[1] = v1;
    } else {
      v[0] = v1; v[1] = v0;
    }
  }
  int V(int i) const { assert (0 <= i && i <= 1); return v[i];}  
};

inline bool operator==(edge const& lhs, edge const& rhs)
{ return (lhs.V(0) == rhs.V(0) && lhs.V(1) == rhs.V(1));}


class triangle {
 int v[3];
public:
  triangle() { v[0] = v[1] = v[2] = -1;}
  triangle(int v0, int v1, int v2)
    { 
      v[0] = v0; v[1] = v1; v[2] = v2; 
    }
  int V(int i) const { assert (0 <= i && i <= 2);  return v[i];}
};


class tetra {
 int v[4];
public:
  tetra() { v[0] = v[1] = v[2] = v[3] = -1;}
  tetra(int v0, int v1, int v2, int v3)
    { 
      v[0] = v0; v[1] = v1; v[2] = v2; v[3] = v3;
    }
  int V(int i) const { assert (0 <= i && i <= 3);  return v[i];}

  edge E(int i) const { 
    assert (0 <= i && i <= 5);
    return edge(v[edgetable[i][0]], v[edgetable[i][1]]);
  }
private:
  typedef int edgevec[2]; 
  static edgevec edgetable[6]; // = { {0,1}, {2,3}, {0,2}, {1,3}, {0,3}, {1,2}};
};


class gridGenerator;

class octcell {
  gridGenerator const* g;
  int x_,y_,z_, d_;
public:
  octcell() : g(0) { x_ = y_ =z_ = d_ = -1;}
  octcell(gridGenerator const& gg,
	  int x, int y, int z,  int d)
    : g(&gg), x_(x), y_(y), z_(z), d_(d) {}

  int x() const   { return x_;}
  int y() const   { return y_;}
  int z() const   { return z_;}
  int dim() const { return d_;}
  /*
  int handle() const { return x*(g->ny)*(g->nz) +y*(g->nz) + z;}
  int material() const { return g->getMaterial(*this);}
  
  struct hasher_type {
    unsigned operator()(octcell const& c) const { 
      assert (c.handle() >= 0); return c.handle();
    }
  };
  */
};

inline bool operator==(octcell const& l, octcell const& r)
{ return l.x() == r.x() && l.y() == r.y() && l.z() == r.z() && l.dim() == r.dim();}


class octcell_tesselation {
private:
  int                            midpoint;
  octcell                        c;
  std::vector<int>               materials;
  std::vector<std::vector<int> > elements;
  
public:
  octcell_tesselation() : midpoint(-1) {}
  octcell_tesselation(octcell const& cc, int mp = -1) : midpoint(mp), c(cc) {}

  void addMidpoint(int mp) { midpoint=mp;}
  bool hasMidpoint() const { return midpoint != -1;}
  int  getMidpoint() const { return midpoint;}

  void addElement(std::vector<int> const& e, int m) 
  { 
    elements .push_back(e); 
    materials.push_back(m); 
  }
  // convenience functions
  void addTetra(int v1, int v2, int v3, int v4, int m) 
  {
    elements.push_back(std::vector<int>(4));
    elements.back()[0] = v1; 
    elements.back()[1] = v2; 
    elements.back()[2] = v3; 
    elements.back()[3] = v4;
    materials.push_back(m); 
  }
  void addCube(int* pl, int m)
  {
    elements.push_back(std::vector<int>(pl,pl+8));
    materials.push_back(m);
  }
  octcell const& TheCell() const { return c;}  

  bool empty() const { return elements.empty();}
  typedef std::vector<std::vector<int> >::const_iterator element_iterator;
  unsigned NumOfElements() const { return elements.size();}
  element_iterator begin_elements() const { return elements.begin();}
  element_iterator end_elements  () const { return elements.end();}

  int getMaterial(element_iterator e) const { return materials[e-elements.begin()];}

};


#endif
