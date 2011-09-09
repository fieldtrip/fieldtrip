#ifndef VEC3_H
#define VEC3_H

/*
 * vec3: 3D vec3 operations
 */

/*  NeuroFEM license:
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
//    double	  max()
//  		  { double d = (x > y? x: y); return (d > z? d: z); };
//    double	  min()
//  		  { double d = (x < y? x: y); return (d < z? d: z); };
//    vec3  	  max(const vec3 &b)
//  		  { return vec3((x>b.x? x: b.x), (y>b.y? y: b.y), (z>b.z? z: b.z)); };
//    vec3  	  min(const vec3 &b)
//  		  { return vec3((x<b.x? x: b.x), (y<b.y? y: b.y), (z<b.z? z: b.z)); };
//    vec3  	  max(const double b)
//  		  { return vec3((x>b? x: b),(y>b? y: b),(z>b? z: b)); };
//    vec3  	  min(const double b)
//  		  { return vec3((x<b? x: b),(y<b? y: b),(z<b? z: b)); };
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
