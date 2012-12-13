#ifndef SIMBIO_VEC3_H
#define SIMBIO_VEC3_H

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


/* \file
 *  
 *   3D vector operations
 *
 *   $Id$
 *
 *   \author F. Kruggel (E-mail kruggel@cns.mpg.de)
 *
 */

#include <math.h> 
#include <stdio.h> 

struct vec3 {
	double		x, y, z;

	vec3() : x(0), y(0), z(0)
			{  };
	vec3(const double v) : x(v), y(v), z(v)
			{ };
	vec3(const double _x, const double _y, const double _z) :
			  x(_x), y(_y), z(_z) { };
	double		scalar() const
			{ return sqrt(x * x + y * y + z * z); };
	double		dot(const vec3& b)
			{ return (x * b.x + y * b.y + z * b.z); };
	vec3		by(const vec3& b)
			{ return vec3(x*b.x, y*b.y, z*b.z); };
	vec3		cross(const vec3& b)
			{ return vec3(y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x); };
	vec3&		operator=(const vec3& b)
			{ memcpy(this, &b, sizeof(vec3)); return *this; };
	vec3&		operator=(const double b)
			{ x = b; y = b; z = b; return *this; };
	double		max()
			{ double d = (x > y? x: y); return (d > z? d: z); };
	double		min()
			{ double d = (x < y? x: y); return (d < z? d: z); };
	vec3		max(const vec3 &b)
			{ return vec3((x>b.x? x: b.x), (y>b.y? y: b.y), (z>b.z? z: b.z)); };
	vec3		min(const vec3 &b)
			{ return vec3((x<b.x? x: b.x), (y<b.y? y: b.y), (z<b.z? z: b.z)); };
	vec3		max(const double b)
			{ return vec3((x>b? x: b),(y>b? y: b),(z>b? z: b)); };
	vec3		min(const double b)
			{ return vec3((x<b? x: b),(y<b? y: b),(z<b? z: b)); };
	void		operator+=(const double b)
			{ x += b; y += b; z += b; };
	void		operator+=(const vec3 &b)
			{ x += b.x; y += b.y; z += b.z; };
	void		operator-=(const vec3 &b)
			{ x -= b.x; y -= b.y; z -= b.z; };
	void		operator-=(const double b)
			{ x -= b; y -= b; z -= b; };
	void		operator*=(const double b)
			{ x *= b; y *= b; z *= b; };
	void		operator/=(const double b)
			{ x /= b; y /= b; z /= b; };
	vec3		normalize()
			{ double sc = scalar();
                          return (sc? vec3(x/sc, y/sc, z/sc): *this); };
	vec3		toPolar()
			{ if (x == 0 || z == 0)  {
				fprintf(stderr, "illegal input for toPolar() operation.\n");
				return vec3();
			  };
			  return vec3(scalar(), atan(y/x), atan(sqrt(x*x+y*y)/z)); };
	vec3		toCartesian()
			{ return vec3(sin(z)*cos(y)*x, sin(z)*sin(y)*x, cos(z)*x); };
	double		which(int i)
			{ return (i == 0? x: (i == 1? y: z)); };
};

inline vec3	operator+(const vec3& a, const vec3& b)
			{ return vec3(a.x+b.x, a.y+b.y, a.z+b.z); }
inline vec3	operator+(const vec3& a, const double b)
			{ return vec3(a.x+b, a.y+b, a.z+b); }
inline vec3	operator+(const double a, const vec3& b)
			{ return vec3(a+b.x, a+b.y, a+b.z); }
inline vec3	operator-(const vec3& a)
			{ return vec3(-a.x, -a.y, -a.z); }
inline vec3	operator-(const vec3& a, const vec3& b)
			{ return vec3(a.x-b.x, a.y-b.y, a.z-b.z); }
inline vec3	operator-(const vec3& a, const double b)
			{ return vec3(a.x-b, a.y-b, a.z-b); }
inline vec3	operator-(const double a, const vec3& b)
			{ return vec3(a-b.x, a-b.y, a-b.z); }
inline vec3	operator*(const double a, const vec3& b)
			{ return vec3(a*b.x, a*b.y, a*b.z); }
inline vec3	operator*(const vec3& a, const double b)
			{ return vec3(a.x*b, a.y*b, a.z*b); }
inline vec3	operator*(const vec3& a, const vec3& b)
			{ return vec3(a.x*b.x, a.y*b.y, a.z*b.z); }
inline vec3	operator/(const vec3& a, const double b)
			{ return vec3(a.x/b, a.y/b, a.z/b); }
inline vec3	operator/(const vec3& a, const vec3 &b)
			{ return vec3(a.x/b.x, a.y/b.y, a.z/b.z); }
inline bool	operator== (const vec3& a, const vec3& b)
			{ return (a.x==b.x && a.y==b.y && a.z==b.z); }
inline bool	operator!= (const vec3& a, const vec3& b)
			{ return (a.x!=b.x || a.y!=b.y || a.z!=b.z); }
inline bool	operator== (const vec3& a, const double b)
			{ return (a.x==b && a.y==b && a.z==b); }
inline bool	operator!= (const vec3& a, const double b)
			{ return (a.x!=b || a.y!=b || a.z!=b); }
inline bool	operator> (const vec3& a, const double b)
			{ return (a.x>b && a.y>b && a.z>b); }
inline bool	operator>= (const vec3& a, const double b)
			{ return (a.x>=b && a.y>=b && a.z>=b); }
inline bool	operator< (const vec3& a, const double b)
			{ return (a.x<b && a.y<b && a.z<b); }
inline bool	operator<= (const vec3& a, const double b)
			{ return (a.x<=b && a.y<=b && a.z<=b); }
inline bool	operator> (const vec3& a, const vec3& b)
			{ return (a.x>b.x && a.y>b.y && a.z>b.z); }
inline bool	operator>= (const vec3& a, const vec3& b)
			{ return (a.x>=b.x && a.y>=b.y && a.z>=b.z); }
inline bool	operator< (const vec3& a, const vec3& b)
			{ return (a.x<b.x && a.y<b.y && a.z<b.z); }
inline bool	operator<= (const vec3& a, const vec3& b)
			{ return (a.x<=b.x && a.y<=b.y && a.z<=b.z); }

#endif
