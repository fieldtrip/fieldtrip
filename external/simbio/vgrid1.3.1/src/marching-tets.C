#include "marching-tets.h"

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

/*! \brief Support for the generalized volumetric marching tetrahedra algorithm

    \author Guntram Berti <berti@ccrl-nece.de>

     $Id: marching-tets.C,v 1.7 2004/02/16 14:49:27 berti Exp $

*/


namespace marching_tets {
  using namespace std;

  subdivision::subdivision(int             pat[4],      // pat[i] in {+1,0,-1}^4 : Function value signs on vertices of t
			   tet             const& t_,   // the tet to be subdivided
			   cut_vertex_map  const& cv)   // the cut vertices on the edges of t
    {
      cutv = &cv;
      t    = t_;
      np = 0; n0 = 0; nm = 0; nc = 0;
      // ns = 0; nv = 0;
      for(int i = 0; i < 4; ++i) {
	np += (pat[i] == +1 ? 1 : 0);
	n0 += (pat[i] ==  0 ? 1 : 0);
	nm += (pat[i] == -1 ? 1 : 0);
      }
      vp_ = 0; 
      v0_ = vp_ + np; 
      vm_ = v0_ + n0;
      vc_ = vm_ + nm;
 
     // sort vertices according to +,0,- and < within each group,
      // keep track of permutation
      int ip = 0, i0 = 0, im = 0;
      v.resize(4);
      for(int i = 0; i < 4; ++i) {
	if(pat[i] == +1) { v[      ip] = t[i]; ++ip;}
	if(pat[i] ==  0) { v[np   +i0] = t[i]; ++i0;}
	if(pat[i] == -1) { v[np+n0+im] = t[i]; ++im;}
      }
      int* b = &v[0];
      std::sort(b,       b+np);
      std::sort(b+np,    b+np+n0);
      std::sort(b+np+n0, b+np+n0+nm);
      
      for(int i = 0; i < 4; ++i)
	for(int j = i+1; j < 4; ++j)
	  if(cutv->find(edge(t[i],t[j])) != cutv->end()) {
	    ++nc;
	    v.push_back( (*(cutv->find(edge(t[i],t[j])))).second);
	  }
	    
      // sanity checks
      // we cannot have more then 4 edges cut
      REQUIRE( nc <= 4, " nc= " << nc,1);
      if(np == 0 || nm == 0)
	// there cannot be cuts when no vertices of opposite sign exist
	REQUIRE( nc == 0, " nc=" << nc, 1);

      create_tets();     
    }    

  int subdivision::Sign(int gv) const 
  {
    int lv = std::find(v.begin(), v.end(), gv) - v.begin();
    REQUIRE(lv != (int)v.size(), "gv = ", 1);
    int s = 77;
    if      (vp_ <= lv && lv < vp_ +np) s = +1;
    else if (vm_ <= lv && lv < vm_ +nm) s = -1;
    else                                s = 0;
    return s;
  }

   void subdivision::create_tets()  {
      if(np == 0 || nm == 0) {
	// nothing to do: just take original tet.
	if(np == 0 && nm == 0) 
	  // all 4 vertices have sign 0: no decision possible for tet.
	  // User has to handle this case
	  add_tet(tet(v[0],v[1],v[2],v[3]),  0);
	else if(nm > 0) {
	  add_tet(tet(v[0],v[1],v[2],v[3]), -1);
	  if(n0 == 3)	
	    add_triangle(triangle(v0(0),v0(1),v0(2)), -1);
	}
	else  { // np > 0 
	  add_tet(tet(v[0],v[1],v[2],v[3]), +1);	  
	  if(n0 == 3)	
	    add_triangle(triangle(v0(0),v0(1),v0(2)), +1);
	}
      }
      // (A) n0 == 0
      else if( np == 1 && n0 == 0 && nm == 3) {
	// case (1,0,3)
	add_tet(tet(vp(0),vc(vp(0),vm(0)),vc(vp(0),vm(1)),vc(vp(0),vm(2))), +1);
	add_tet(tet(vm(0),vc(vp(0),vm(0)),vc(vp(0),vm(1)),vc(vp(0),vm(2))), -1);
	add_tet(tet(vm(0),vm(1),          vm(2),          vc(vp(0),vm(2))), -1);
	add_tet(tet(vm(0),vm(1),          vc(vp(0),vm(1)),vc(vp(0),vm(2))), -1);
	
	add_triangle(triangle(vc(vp(0),vm(0)),vc(vp(0),vm(1)),vc(vp(0),vm(2))), 0);
      }
      else if( np == 3 && n0 == 0 && nm == 1) {
	// case (3,0,1): symmetric to (1,0,3), vm <-> vp, +1 <-> -1
	add_tet(tet(vm(0),vc(vm(0),vp(0)),vc(vm(0),vp(1)),vc(vm(0),vp(2))), -1);
	add_tet(tet(vp(0),vc(vm(0),vp(0)),vc(vm(0),vp(1)),vc(vm(0),vp(2))), +1);
	add_tet(tet(vp(0),vp(1),          vp(2),          vc(vm(0),vp(2))), +1);
	add_tet(tet(vp(0),vp(1),          vc(vm(0),vp(1)),vc(vm(0),vp(2))), +1);
	
	add_triangle(triangle(vc(vm(0),vp(0)),vc(vm(0),vp(1)),vc(vm(0),vp(2))), 0);
      }
      else if( np == 2 && n0 == 0 && nm == 2) {
	// case (2,0,2)
	add_tet(tet(vp(0), vp(1),            vc(vp(1),vm(0)), vc(vp(1), vm(1))), +1);
	add_tet(tet(vp(0), vc(vp(0), vm(1)), vc(vp(1),vm(0)), vc(vp(0), vm(0))), +1);
	add_tet(tet(vp(0), vc(vp(0), vm(1)), vc(vp(1),vm(0)), vc(vp(1), vm(1))), +1);

	add_tet(tet(vm(0), vm(1),            vc(vm(1),vp(0)), vc(vm(1), vp(1))), -1);
	add_tet(tet(vm(0), vc(vp(0), vm(1)), vc(vp(1),vm(0)), vc(vp(0), vm(0))), -1);
	add_tet(tet(vm(0), vc(vp(0), vm(1)), vc(vp(1),vm(0)), vc(vp(1), vm(1))), -1);

	add_triangle(triangle(vc(vp(0), vm(1)), vc(vp(1),vm(0)), vc(vp(0), vm(0))), 0);
	add_triangle(triangle(vc(vp(0), vm(1)), vc(vp(1),vm(0)), vc(vp(1), vm(1))), 0);
      }

      // (B) n0 == 1
      else if( np == 1 && n0 == 1 && nm == 2) {
        add_tet(tet(vp(0), v0(0), vc(vp(0), vm(0)), vc(vp(0), vm(1))), +1);
	add_tet(tet(vm(0), v0(0), vc(vp(0), vm(0)), vc(vp(0), vm(1))), -1);
	add_tet(tet(vm(0), v0(0), vm(1),            vc(vp(0), vm(1))), -1);
	
	add_triangle(triangle(v0(0), vc(vp(0), vm(0)), vc(vp(0), vm(1))), 0);
      }
      else if( np == 2 && n0 == 1 && nm == 1) {
        add_tet(tet(vm(0), v0(0), vc(vm(0), vp(0)), vc(vm(0), vp(1))), -1);
	add_tet(tet(vp(0), v0(0), vc(vm(0), vp(0)), vc(vm(0), vp(1))), +1);
	add_tet(tet(vp(0), v0(0), vp(1),            vc(vm(0), vp(1))), +1);

	add_triangle(triangle(v0(0), vc(vm(0), vp(0)), vc(vm(0), vp(1))), 0);
      }
      // other cases caught by np == 0 || nm == 0

      // (C) n0 == 2
      else if( np == 1 && n0 == 2 && nm == 1) {
        add_tet(tet(vp(0), v0(0), v0(1), vc(vp(0), vm(0))), +1);
        add_tet(tet(vm(0), v0(0), v0(1), vc(vp(0), vm(0))), -1);
	
	add_triangle(triangle(v0(0), v0(1), vc(vp(0), vm(0))), 0);
      }
      // other cases caught by np == 0 || nm == 0

      // n0 == 3
      // caught by np == 0 || nm == 0

      // n0 == 4
      // caught by np == 0 || nm == 0



      // other case should not occur
      else
	cerr << "Undefined pattern: " << np << ',' << n0 << ',' << nm << endl; 
  }



}
