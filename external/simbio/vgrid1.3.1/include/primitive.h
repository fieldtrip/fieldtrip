#ifndef SIMBIO_PRIMITIVE_H
#define SIMBIO__PRIMITIVE_H

/*
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

#include "vista.h"

class primitive : public VNodeBaseRec  {
	friend class 	gridGenerator;
	friend class 	delaunay;
	friend class 	mesh;
       
public:
		VLong		vcnt;
        	VLong		id[8];
        	primitive() 	{ hops = 0; visited = 0; head = 0; weight = 0; }
};


#endif


