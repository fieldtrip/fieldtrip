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

#ifndef SIMBIO_WP1_FIL_PRIMITIVE_H
#define SIMBIO_WP1_FIL_PRIMITIVE_H

#include "FIL_Vista.h"

class FIL_Primitive : public VNodeBaseRec  {
	friend class 	gridGenerator;
	friend class 	delaunay;
	friend class 	FIL_Mesh;
       
public:
		VLong		vcnt;
        	VLong		id[8];
        	FIL_Primitive() 	{ hops = 0; visited = 0; head = 0; weight = 0; }
};


#endif


