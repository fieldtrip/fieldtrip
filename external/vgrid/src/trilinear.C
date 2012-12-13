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
 * trilinear interpolation on Vista images
 *
 * $Id: trilinear.C,v 1.1 2004/02/16 12:56:48 berti Exp $
 */ 


#include <vista.h>
#include <vec3.h>

static inline vec3 get(VImage grad, int x, int y, int z)
{
	vec3 v;
	v.x = (double)VPixel(grad, z*3+0, y, x, VFloat);
	v.y = (double)VPixel(grad, z*3+1, y, x, VFloat);
	v.z = (double)VPixel(grad, z*3+2, y, x, VFloat);
	return v;
}

vec3 getForce(VImage grad, const vec3 &v)
// trilinear vec3 interpolation
{
	int nx = VImageNColumns(grad);
	int ny = VImageNRows(grad);
	int nz = VImageNFrames(grad);

	// dont bother with outside and border pixels
	int x = (int)v.x; int y = (int)v.y; int z = (int)v.z;
       	if (x < 1 || x >= nx-1) return 0;
       	if (y < 1 || y >= ny-1) return 0;
       	if (z < 1 || z >= nz-1) return 0;

	// but now we know all coordinates are inside our volume
	double dx = v.x-x; double dy = v.y-y; double dz = v.z-z;
	double xm = 1-dx;  double ym = 1-dy;  double zm = 1-dz;
	vec3 p;
	p  = (xm * ym * zm) * get(grad, x+0, y+0, z+0);
	p += (dx * ym * zm) * get(grad, x+1, y+0, z+0);
	p += (xm * dy * zm) * get(grad, x+0, y+1, z+0);
	p += (dx * dy * zm) * get(grad, x+1, y+1, z+0);
	p += (xm * ym * dz) * get(grad, x+0, y+0, z+1);
	p += (dx * ym * dz) * get(grad, x+1, y+0, z+1);
	p += (xm * dy * dz) * get(grad, x+0, y+1, z+1);
	p += (dx * dy * dz) * get(grad, x+1, y+1, z+1);
	return p; 
}

double getPixel(VImage image, vec3 &v)
// trilinear intensity interpolation
{
	// dont bother with outside and border pixels
	int x = (int)v.x; int y = (int)v.y; int z = (int)v.z;
	int nx = VImageNColumns(image);
	int ny = VImageNRows(image);
	int nz = VImageNBands(image);
       	if (x < 1 || x >= nx-1) return 0;
       	if (y < 1 || y >= ny-1) return 0;
       	if (z < 1 || z >= nz-1) return 0;

	// but now we know all coordinates are inside our volume
	double dx = v.x-x; double dy = v.y-y; double dz = v.z-z;
	double xm = 1-dx;  double ym = 1-dy;  double zm = 1-dz;
	double p;
	p  = xm * ym * zm * VPixel(image, z+0, y+0, x+0, VUByte);
	p += dx * ym * zm * VPixel(image, z+0, y+0, x+1, VUByte);
	p += xm * dy * zm * VPixel(image, z+0, y+1, x+0, VUByte);
	p += dx * dy * zm * VPixel(image, z+0, y+1, x+1, VUByte);
	p += xm * ym * dz * VPixel(image, z+1, y+0, x+0, VUByte);
	p += dx * ym * dz * VPixel(image, z+1, y+0, x+1, VUByte);
	p += xm * dy * dz * VPixel(image, z+1, y+1, x+0, VUByte);
	p += dx * dy * dz * VPixel(image, z+1, y+1, x+1, VUByte);
	return p; 
}
