/*
 * $Id: spm_invdef.c 938 2007-10-12 19:09:31Z john $
 * John Ashburner
 */

/*  Invert a deformation field.

   FORMAT: [IY1,IY2,IY3] = spm_invdef(Y1,Y2,Y3,dim,M1,M2);

	IY1,IY2,IY3 - inverted deformation field of size dim(1)*dim(2)*dim(3).
	Y1,Y2,Y3    - original deformation field.
	M1          - An affine mapping from mm to voxels in the co-ordinate
	              system of the inverse deformation field.
	M2          - An affine mapping from voxels to mm in the co-ordinate
	              system of the forward deformation field.

    The field is assumed to consist of a piecewise affine transformations, whereby
    each cube jointing 8 neighbouring voxels contains 8 tetrahedra.  The mapping
    within each tetrahedron is assumed to be affine.
    More documentation can be found in:
	J. Ashburner, J. Andersson and K. J. Friston (2000).
	"Image Registration using a Symmetric Prior - in Three-Dimensions".
	Human Brain Mapping 9(4):212-225.
*/

#include <math.h>
#include "mex.h"

#define MAXV 2048
#define REAL float

static void invertX(REAL X[4][3], REAL IX[4][4])
/* X is a matrix containing the co-ordinates of the four vertices of a tetrahedron.
   IX = inv([X ; 1 1 1 1]);  */
{
	REAL id;
	id = X[0][0]*(X[3][1]*(X[1][2]-X[2][2])+X[1][1]*(X[2][2]-X[3][2])+X[2][1]*(X[3][2]-X[1][2]))+
	     X[1][0]*(X[3][2]*(X[0][1]-X[2][1])+X[0][2]*(X[2][1]-X[3][1])+X[2][2]*(X[3][1]-X[0][1]))+
	     X[2][0]*(X[0][1]*(X[1][2]-X[3][2])+X[3][1]*(X[0][2]-X[1][2])+X[1][1]*(X[3][2]-X[0][2]))+
	     X[3][0]*(X[1][2]*(X[2][1]-X[0][1])+X[0][2]*(X[1][1]-X[2][1])+X[2][2]*(X[0][1]-X[1][1]));
	id = 1/id;
	IX[0][0] = id*(X[1][1]*(X[2][2]-X[3][2])+X[2][1]*(X[3][2]-X[1][2])+X[3][1]*(X[1][2]-X[2][2]));
	IX[0][1] = id*(X[0][1]*(X[3][2]-X[2][2])+X[2][1]*(X[0][2]-X[3][2])+X[3][1]*(X[2][2]-X[0][2]));
	IX[0][2] = id*(X[0][1]*(X[1][2]-X[3][2])+X[1][1]*(X[3][2]-X[0][2])+X[3][1]*(X[0][2]-X[1][2]));
	IX[0][3] = id*(X[0][1]*(X[2][2]-X[1][2])+X[1][1]*(X[0][2]-X[2][2])+X[2][1]*(X[1][2]-X[0][2]));
	IX[1][0] = id*(X[1][0]*(X[3][2]-X[2][2])+X[2][0]*(X[1][2]-X[3][2])+X[3][0]*(X[2][2]-X[1][2]));
	IX[1][1] = id*(X[0][0]*(X[2][2]-X[3][2])+X[2][0]*(X[3][2]-X[0][2])+X[3][0]*(X[0][2]-X[2][2]));
	IX[1][2] = id*(X[0][0]*(X[3][2]-X[1][2])+X[1][0]*(X[0][2]-X[3][2])+X[3][0]*(X[1][2]-X[0][2]));
	IX[1][3] = id*(X[0][0]*(X[1][2]-X[2][2])+X[1][0]*(X[2][2]-X[0][2])+X[2][0]*(X[0][2]-X[1][2]));
	IX[2][0] = id*(X[1][0]*(X[2][1]-X[3][1])+X[2][0]*(X[3][1]-X[1][1])+X[3][0]*(X[1][1]-X[2][1]));
	IX[2][1] = id*(X[0][0]*(X[3][1]-X[2][1])+X[2][0]*(X[0][1]-X[3][1])+X[3][0]*(X[2][1]-X[0][1]));
	IX[2][2] = id*(X[0][0]*(X[1][1]-X[3][1])+X[1][0]*(X[3][1]-X[0][1])+X[3][0]*(X[0][1]-X[1][1]));
	IX[2][3] = id*(X[0][0]*(X[2][1]-X[1][1])+X[1][0]*(X[0][1]-X[2][1])+X[2][0]*(X[1][1]-X[0][1]));
	IX[3][0] = id*(X[1][0]*(X[2][2]*X[3][1]-X[3][2]*X[2][1])+
	               X[2][0]*(X[3][2]*X[1][1]-X[1][2]*X[3][1])+
		       X[3][0]*(X[1][2]*X[2][1]-X[2][2]*X[1][1]));
	IX[3][1] = id*(X[0][0]*(X[3][2]*X[2][1]-X[2][2]*X[3][1])+
		       X[2][0]*(X[0][2]*X[3][1]-X[3][2]*X[0][1])+
		       X[3][0]*(X[2][2]*X[0][1]-X[0][2]*X[2][1]));
	IX[3][2] = id*(X[0][0]*(X[1][2]*X[3][1]-X[3][2]*X[1][1])+
		       X[1][0]*(X[3][2]*X[0][1]-X[0][2]*X[3][1])+
		       X[3][0]*(X[0][2]*X[1][1]-X[1][2]*X[0][1]));
	IX[3][3] = id*(X[0][0]*(X[2][2]*X[1][1]-X[1][2]*X[2][1])+
		       X[1][0]*(X[0][2]*X[2][1]-X[2][2]*X[0][1])+
		       X[2][0]*(X[1][2]*X[0][1]-X[0][2]*X[1][1]));
}

static void getM(REAL Y[4][3], REAL IX[4][4], REAL M[4][3], int i, int j, int k)
/* Determines the affine transform (M) mapping from
   [X+repmat([i j k]', 1,4) ; 1 1 1 1] to [Y ; 1 1 1 1], where IX = inv([X ; 1 1 1 1]);
   This is given by:
        M = Y*inv([X+repmat([i j k]', 1,4) ; 1 1 1 1]);
   or more efficiently by:
        M = Y*(IX - IX(:,1:3)*[i j k]'); */
{
	REAL ix30, ix31, ix32, ix33;

	ix30 = IX[3][0] - (i*IX[0][0] + j*IX[1][0] + k*IX[2][0]);
	ix31 = IX[3][1] - (i*IX[0][1] + j*IX[1][1] + k*IX[2][1]);
	ix32 = IX[3][2] - (i*IX[0][2] + j*IX[1][2] + k*IX[2][2]);
	ix33 = IX[3][3] - (i*IX[0][3] + j*IX[1][3] + k*IX[2][3]);

	M[0][0] = Y[0][0]*IX[0][0] + Y[1][0]*IX[0][1] + Y[2][0]*IX[0][2] + Y[3][0]*IX[0][3];
	M[0][1] = Y[0][1]*IX[0][0] + Y[1][1]*IX[0][1] + Y[2][1]*IX[0][2] + Y[3][1]*IX[0][3];
	M[0][2] = Y[0][2]*IX[0][0] + Y[1][2]*IX[0][1] + Y[2][2]*IX[0][2] + Y[3][2]*IX[0][3];

	M[1][0] = Y[0][0]*IX[1][0] + Y[1][0]*IX[1][1] + Y[2][0]*IX[1][2] + Y[3][0]*IX[1][3];
	M[1][1] = Y[0][1]*IX[1][0] + Y[1][1]*IX[1][1] + Y[2][1]*IX[1][2] + Y[3][1]*IX[1][3];
	M[1][2] = Y[0][2]*IX[1][0] + Y[1][2]*IX[1][1] + Y[2][2]*IX[1][2] + Y[3][2]*IX[1][3];

	M[2][0] = Y[0][0]*IX[2][0] + Y[1][0]*IX[2][1] + Y[2][0]*IX[2][2] + Y[3][0]*IX[2][3];
	M[2][1] = Y[0][1]*IX[2][0] + Y[1][1]*IX[2][1] + Y[2][1]*IX[2][2] + Y[3][1]*IX[2][3];
	M[2][2] = Y[0][2]*IX[2][0] + Y[1][2]*IX[2][1] + Y[2][2]*IX[2][2] + Y[3][2]*IX[2][3];

	M[3][0] = Y[0][0]*ix30     + Y[1][0]*ix31     + Y[2][0]*ix32     + Y[3][0]*ix33;
	M[3][1] = Y[0][1]*ix30     + Y[1][1]*ix31     + Y[2][1]*ix32     + Y[3][1]*ix33;
	M[3][2] = Y[0][2]*ix30     + Y[1][2]*ix31     + Y[2][2]*ix32     + Y[3][2]*ix33;
}

static void mulMM(REAL A[4][3], REAL B[4][3], REAL C[4][3])
/* [A ; 0 0 0 1] = [B ; 0 0 0 1]*[C ; 0 0 0 1]; */
{
	A[0][0] = B[0][0]*C[0][0] + B[1][0]*C[0][1] + B[2][0]*C[0][2];
	A[0][1] = B[0][1]*C[0][0] + B[1][1]*C[0][1] + B[2][1]*C[0][2];
	A[0][2] = B[0][2]*C[0][0] + B[1][2]*C[0][1] + B[2][2]*C[0][2];

	A[1][0] = B[0][0]*C[1][0] + B[1][0]*C[1][1] + B[2][0]*C[1][2];
	A[1][1] = B[0][1]*C[1][0] + B[1][1]*C[1][1] + B[2][1]*C[1][2];
	A[1][2] = B[0][2]*C[1][0] + B[1][2]*C[1][1] + B[2][2]*C[1][2];

	A[2][0] = B[0][0]*C[2][0] + B[1][0]*C[2][1] + B[2][0]*C[2][2];
	A[2][1] = B[0][1]*C[2][0] + B[1][1]*C[2][1] + B[2][1]*C[2][2];
	A[2][2] = B[0][2]*C[2][0] + B[1][2]*C[2][1] + B[2][2]*C[2][2];

	A[3][0] = B[0][0]*C[3][0] + B[1][0]*C[3][1] + B[2][0]*C[3][2] + B[3][0];
	A[3][1] = B[0][1]*C[3][0] + B[1][1]*C[3][1] + B[2][1]*C[3][2] + B[3][1];
	A[3][2] = B[0][2]*C[3][0] + B[1][2]*C[3][1] + B[2][2]*C[3][2] + B[3][2];
}

static void mulMX(REAL A[4][3], REAL B[4][3], REAL C[4][3])
/* [A ; 1 1 1 1] = [B ; 0 0 0 1]*[C ; 1 1 1 1]; */
{
	A[0][0] = B[0][0]*C[0][0] + B[1][0]*C[0][1] + B[2][0]*C[0][2] + B[3][0];
	A[0][1] = B[0][1]*C[0][0] + B[1][1]*C[0][1] + B[2][1]*C[0][2] + B[3][1];
	A[0][2] = B[0][2]*C[0][0] + B[1][2]*C[0][1] + B[2][2]*C[0][2] + B[3][2];

	A[1][0] = B[0][0]*C[1][0] + B[1][0]*C[1][1] + B[2][0]*C[1][2] + B[3][0];
	A[1][1] = B[0][1]*C[1][0] + B[1][1]*C[1][1] + B[2][1]*C[1][2] + B[3][1];
	A[1][2] = B[0][2]*C[1][0] + B[1][2]*C[1][1] + B[2][2]*C[1][2] + B[3][2];

	A[2][0] = B[0][0]*C[2][0] + B[1][0]*C[2][1] + B[2][0]*C[2][2] + B[3][0];
	A[2][1] = B[0][1]*C[2][0] + B[1][1]*C[2][1] + B[2][1]*C[2][2] + B[3][1];
	A[2][2] = B[0][2]*C[2][0] + B[1][2]*C[2][1] + B[2][2]*C[2][2] + B[3][2];

	A[3][0] = B[0][0]*C[3][0] + B[1][0]*C[3][1] + B[2][0]*C[3][2] + B[3][0];
	A[3][1] = B[0][1]*C[3][0] + B[1][1]*C[3][1] + B[2][1]*C[3][2] + B[3][1];
	A[3][2] = B[0][2]*C[3][0] + B[1][2]*C[3][1] + B[2][2]*C[3][2] + B[3][2];
}

static void invertM(REAL M[4][3], REAL IM[4][3])
/* [IM ; 0 0 0 1] = inv([M ; 0 0 0 1]); */
{
	REAL id;
	id = M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])+
	     M[0][1]*(M[1][2]*M[2][0]-M[1][0]*M[2][2])+
	     M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]);

	id = 1.0/id;
	IM[0][0] = (M[1][1]*M[2][2]-M[1][2]*M[2][1])*id;
	IM[0][1] = (M[0][2]*M[2][1]-M[0][1]*M[2][2])*id;
	IM[0][2] = (M[0][1]*M[1][2]-M[0][2]*M[1][1])*id;

	IM[1][0] = (M[1][2]*M[2][0]-M[1][0]*M[2][2])*id;
	IM[1][1] = (M[0][0]*M[2][2]-M[0][2]*M[2][0])*id;
	IM[1][2] = (M[0][2]*M[1][0]-M[0][0]*M[1][2])*id;

	IM[2][0] = (M[1][0]*M[2][1]-M[1][1]*M[2][0])*id;
	IM[2][1] = (M[0][1]*M[2][0]-M[0][0]*M[2][1])*id;
	IM[2][2] = (M[0][0]*M[1][1]-M[0][1]*M[1][0])*id;

	IM[3][0] = (M[1][0]*(M[3][1]*M[2][2]-M[2][1]*M[3][2])+
		    M[1][1]*(M[2][0]*M[3][2]-M[3][0]*M[2][2])+
		    M[1][2]*(M[3][0]*M[2][1]-M[2][0]*M[3][1]))*id;
	IM[3][1] = (M[0][0]*(M[2][1]*M[3][2]-M[3][1]*M[2][2])+
		    M[0][1]*(M[3][0]*M[2][2]-M[2][0]*M[3][2])+
		    M[0][2]*(M[2][0]*M[3][1]-M[3][0]*M[2][1]))*id;
	IM[3][2] = (M[0][0]*(M[3][1]*M[1][2]-M[1][1]*M[3][2])+
		    M[0][1]*(M[1][0]*M[3][2]-M[3][0]*M[1][2])+
		    M[0][2]*(M[3][0]*M[1][1]-M[1][0]*M[3][1]))*id;
}

/****************************************************************************************************/
/* These routines are for locating integer co-ordinates that lie inside a tetrahedron.  See:
	J. Ashburner, J. Andersson and K. J. Friston (2000).
	"Image Registration using a Symmetric Prior - in Three-Dimensions".
	Human Brain Mapping 9(4):212-225. */

static void scan_line(REAL lin[2], int y, int z, int *n, int vox[][3], int maxv)
{
	REAL p[2], t;
	int x, xe;

	/* sort p into ascending order of x */
	p[0] = lin[0]; p[1] = lin[1];
	if (p[1]<p[0]) {t = p[1]; p[1] = p[0]; p[0] = t;}

	/* find voxels where x is integer */
	for(x=ceil(p[0]), xe=floor(p[1]); x<=xe; x++)
	{
		if ((*n)>=maxv-1)
			mexErrMsgTxt("Too many voxels inside a tetrahedron");
		vox[*n][0] = x;
		vox[*n][1] = y;
		vox[*n][2] = z;
		(*n)++;
	}
}

static void scan_triangle(REAL tri[][2], int z, int *n, int vox[][3], int maxv)
{
	REAL *p[3], *t, lin[2];
	REAL x1, x2, y1, y2;
	int y, ye, i;

	/* sort p into ascending order of y */
	p[0] = tri[0]; p[1] = tri[1]; p[2] = tri[2];
	if (p[1][1]<p[0][1]) {t = p[1]; p[1] = p[0]; p[0] = t;}
	if (p[2][1]<p[1][1]) {t = p[2]; p[2] = p[1]; p[1] = t;}
	if (p[1][1]<p[0][1]) {t = p[1]; p[1] = p[0]; p[0] = t;}

	/* find lower lines cutting triangle where y is integer */
	for(y=ceil(p[0][1]), ye=floor(p[1][1]); y<=ye; y++)
	{
		x1 = p[0][0]; y1 = p[0][1];
		for(i=0; i<2; i++)
		{
			x2 = p[i+1][0]; y2 = p[i+1][1];
			if (y2-y1<=0)
				lin[i] = (x1+x2)/2.0;
			else
				lin[i] = (x1*(y2-y)+x2*(y-y1))/(y2-y1);
		}
		scan_line(lin,y,z, n,vox,maxv);
	}

	/* find upper lines cutting triangle where y is integer */
	for(y=ceil(p[1][1]), ye=floor(p[2][1]); y<=ye; y++)
	{
		x2 = p[2][0]; y2 = p[2][1];
		for(i=0; i<2; i++)
		{
			x1 = p[i][0]; y1 = p[i][1];
			if (y2-y1<=0)
				lin[i] = (x1+x2)/2.0;
			else
				lin[i] = (x1*(y2-y)+x2*(y-y1))/(y2-y1);
		}
		scan_line(lin,y,z, n,vox,maxv);
	}
}


static void scan_tetrahedron(REAL Y[4][3], int *n, int vox[][3], int maxv)
/* Y are the vertices of the tetrahedron.  n are the number of located
   integer co-ordinates, vox are the co-ordinates found and maxv are the
   maximum number of co-ordinates allowed. */
{
	REAL *p[4], *t;
	REAL tri[4][2];
	REAL x1, x2, y1, y2, z1, z2;
	int z, ze, i;

	*n = 0;

	/* sort p into ascending order of z */
	p[0] = Y[0]; p[1] = Y[1]; p[2] = Y[2]; p[3] = Y[3];
	if (p[1][2]<p[0][2]) {t = p[1]; p[1] = p[0]; p[0] = t;}
	if (p[2][2]<p[1][2]) {t = p[2]; p[2] = p[1]; p[1] = t;}
	if (p[3][2]<p[2][2]) {t = p[3]; p[3] = p[2]; p[2] = t;}
	if (p[1][2]<p[0][2]) {t = p[1]; p[1] = p[0]; p[0] = t;}
	if (p[2][2]<p[1][2]) {t = p[2]; p[2] = p[1]; p[1] = t;}
	if (p[1][2]<p[0][2]) {t = p[1]; p[1] = p[0]; p[0] = t;}

	/* find lower triangles that intersect tetrahedron where z is integer */
	for(z=ceil(p[0][2]), ze=floor(p[1][2]); z<=ze; z++)
	{
		x1 = p[0][0]; y1 = p[0][1]; z1 = p[0][2];
		for(i=0; i<3; i++)
		{
			x2 = p[i+1][0]; y2 = p[i+1][1]; z2 = p[i+1][2];
			if (z2-z1<=0)
			{
				tri[i][0] = (x1+x2)/2.0;
				tri[i][1] = (y1+y2)/2.0;
			}
			else
			{
				REAL t2 = z2-z, t1 = z-z1, t = z2-z1;
				tri[i][0] = (x1*t2+x2*t1)/t;
				tri[i][1] = (y1*t2+y2*t1)/t;
			}
		}
		scan_triangle(tri,z, n,vox,maxv);
	}

	/* find quadrilaterals that intersect tetrahedron where z is integer */
	/* each quadrilateral divided into two triangles */
	for(z=ceil(p[1][2]), ze=floor(p[2][2]); z<=ze; z++)
	{
		static int ii[] = {0,1,1,0}, jj[] = {3,3,2,2};

		for(i=0; i<4; i++)
		{
			x1 = p[ii[i]][0]; y1 = p[ii[i]][1]; z1 = p[ii[i]][2];
			x2 = p[jj[i]][0]; y2 = p[jj[i]][1]; z2 = p[jj[i]][2];
			if (z2-z1<=0)
			{
				tri[i][0] = (x1+x2)/2.0;
				tri[i][1] = (y1+y2)/2.0;
			}
			else
			{
				REAL t2 = z2-z, t1 = z-z1, t = z2-z1;
				tri[i][0] = (x1*t2+x2*t1)/t;
				tri[i][1] = (y1*t2+y2*t1)/t;
			}
		}
		scan_triangle(tri,z, n,vox,maxv);
		tri[1][0] = tri[3][0];
		tri[1][1] = tri[3][1];
		scan_triangle(tri,z, n,vox,maxv);
	}

	/* find upper triangles that intersect tetrahedron where z is integer */
	for(z=ceil(p[2][2]), ze=floor(p[3][2]); z<=ze; z++)
	{
		x2 = p[3][0]; y2 = p[3][1]; z2 = p[3][2];
		for(i=0; i<3; i++)
		{
			x1 = p[i][0]; y1 = p[i][1]; z1 = p[i][2];
			if (z2-z1<=0)
			{
				tri[i][0] = (x1+x2)/2.0;
				tri[i][1] = (y1+y2)/2.0;
			}
			else
			{
				REAL t2 = z2-z, t1 = z-z1, t = z2-z1;
				tri[i][0] = (x1*t2+x2*t1)/t;
				tri[i][1] = (y1*t2+y2*t1)/t;
			}
		}
		scan_triangle(tri,z, n,vox,maxv);
	}
}

/****************************************************************************************************/

/* Division of a cube into two alternating patterns of tetrahedra.
  This pattern is repeated in a 3D checkerboard pattern, such that
  the whole volume is covered. */
static REAL x[2][5][4][3] = {
{{{ 0,0,0},{ 1,0,1},{ 1,0,0},{ 1,1,0}},
 {{ 0,0,0},{ 1,0,1},{ 0,1,1},{ 0,0,1}},
 {{ 0,0,0},{ 0,1,0},{ 0,1,1},{ 1,1,0}},
 {{ 0,0,0},{ 1,0,1},{ 1,1,0},{ 0,1,1}},
 {{ 1,1,1},{ 1,1,0},{ 0,1,1},{ 1,0,1}},},

{{{ 1,0,0},{ 0,0,1},{ 0,0,0},{ 0,1,0}},
 {{ 1,0,0},{ 0,0,1},{ 1,1,1},{ 1,0,1}},
 {{ 1,0,0},{ 1,1,0},{ 1,1,1},{ 0,1,0}},
 {{ 1,0,0},{ 0,0,1},{ 0,1,0},{ 1,1,1}},
 {{ 0,1,1},{ 0,1,0},{ 1,1,1},{ 0,0,1}},},
};


/* Set up matrices (IX) for each tetrahedron, so that future computations can be
   made faster.  Also set up a few relative file offsets. */
static REAL ix[2][5][4][4];
static int off[2][4][5];
static void setup_consts(int dim[3])
{
	int i, j, k;
	for(k=0; k<2; k++)
		for(i=0; i<5; i++)
		{
			invertX(x[k][i], ix[k][i]);
			for(j=0; j<4; j++)
				off[k][j][i] = x[k][i][j][0]+dim[0]*(x[k][i][j][1]+dim[1]*x[k][i][j][2]);
		}
}

/* Compute the inverse deformation field within a single cube */
static void invert_it(int x0, int x1, int x2, float *y0, float *y1, float *y2,
	int dim_f[3], float *iy0, float *iy1, float *iy2, REAL U[4][3], REAL V[4][3])
{
	int i, j, k, vox[MAXV][3], nvox;
	REAL Y0[4][3], Y[4][3], M[4][3], IM[4][3];

	/* Determine tetrahedral arrangement */
        k = (x0%2)==(x1%2)==(x2%2);
 
	for(i=0; i<5; i++) /* Five tetrahedra within a cube */
	{
		/* Find the vertices (in mm space) */
		Y0[0][0] = y0[off[k][0][i]]; Y0[0][1] = y1[off[k][0][i]]; Y0[0][2] = y2[off[k][0][i]];
		Y0[1][0] = y0[off[k][1][i]]; Y0[1][1] = y1[off[k][1][i]]; Y0[1][2] = y2[off[k][1][i]];
		Y0[2][0] = y0[off[k][2][i]]; Y0[2][1] = y1[off[k][2][i]]; Y0[2][2] = y2[off[k][2][i]];
		Y0[3][0] = y0[off[k][3][i]]; Y0[3][1] = y1[off[k][3][i]]; Y0[3][2] = y2[off[k][3][i]];

		/* Convert vertex co-ordinates to voxels */
		mulMX(Y, U, Y0);

		/* Compute affine transform mapping vertices together */
		getM(Y, ix[k][i], M, x0, x1, x2);

		if (mxIsFinite(M[0][0])) /* Prevent from bombing out when NaNs are encountered */
		{
			/* Find integer co-ordinates within tetrahedron */
			scan_tetrahedron(Y, &nvox, vox, MAXV);

			if (nvox>0)
			{
				/* Invert the affine mapping */
				invertM(M, IM);

				/* Convert the mapping from voxels to mm */
				mulMM(M, V, IM);

				/* Insert the new mappings into each voxel within the tetrahedron */
				for(j=0; j<nvox; j++)
				{
					if ((vox[j][0]>=1) && (vox[j][0]<=dim_f[0]) &&
				 	    (vox[j][1]>=1) && (vox[j][1]<=dim_f[1]) &&
					    (vox[j][2]>=1) && (vox[j][2]<=dim_f[2]))
					{
						int o  = vox[j][0]+dim_f[0]*(vox[j][1]+dim_f[1]*vox[j][2]);

						iy0[o] = M[0][0]*vox[j][0] + M[1][0]*vox[j][1] + M[2][0]*vox[j][2] + M[3][0];
						iy1[o] = M[0][1]*vox[j][0] + M[1][1]*vox[j][1] + M[2][1]*vox[j][2] + M[3][1];
						iy2[o] = M[0][2]*vox[j][0] + M[1][2]*vox[j][1] + M[2][2]*vox[j][2] + M[3][2];
					}
				}
			}
		}
	}
}



static void invert_field(int dim_g[3], float  y0[], float  y1[], float  y2[],
		         int dim_f[3], float iy0[], float iy1[], float iy2[],  REAL U[4][3],  REAL V[4][3])
{
	int x2, x1, x0;

	setup_consts(dim_g);

	/* Convert arrays such that the first voxel is at [1 1 1] rather than [0 0 0].
	   This is a trick used by f2c. */
	y0  -= 1+dim_g[0]*(1 + dim_g[1]);
	y1  -= 1+dim_g[0]*(1 + dim_g[1]);
	y2  -= 1+dim_g[0]*(1 + dim_g[1]);
	iy0 -= 1+dim_f[0]*(1 + dim_f[1]);
	iy1 -= 1+dim_f[0]*(1 + dim_f[1]);
	iy2 -= 1+dim_f[0]*(1 + dim_f[1]);

	/* Loop over all cubes in the deformation field. */
	for(x2=1; x2<dim_g[2]; x2++)
	{
		for(x1=1; x1<dim_g[1]; x1++)
			for(x0=1; x0<dim_g[0]; x0++)
			{
				int o = x0 + dim_g[0]*(x1 + x2*dim_g[1]);
				invert_it(x0, x1, x2, y0+o, y1+o, y2+o, dim_f, iy0, iy1, iy2, U, V);
			}
	}
}

/* Some regions of the inverse deformation may be undefined. */
static void setnan(float *dat, int n)
{
	int j;
	float NaN;
	NaN = mxGetNaN();
	for (j=0; j<n; j++) dat[j] = NaN;
}


static float *get_volume(const mxArray *ptr, int dims[3])
{
	int nd, i;
	const int *ldims;
	if (mxIsStruct(ptr) || !mxIsNumeric(ptr) || mxIsComplex(ptr) ||
                mxIsSparse(ptr) || !mxIsSingle(ptr))
		mexErrMsgTxt("Data must be a single precision floating point multi-dimensional array.");

	nd = mxGetNumberOfDimensions(ptr);
	if (nd>3) mexErrMsgTxt("Too many dimensions in data.");
	ldims = mxGetDimensions(ptr);
	for(i=0; i<nd; i++) dims[i] = ldims[i];
	for(i=nd; i<3; i++) dims[i] = 1;
	return((float *)mxGetPr(ptr));
}

static void get_mat(const mxArray *ptr, REAL M[4][3])
{
	int i, j;
	double *p;

        if (!mxIsNumeric(ptr) || mxIsComplex(ptr) ||
                mxIsComplex(ptr) || !mxIsDouble(ptr) || mxGetM(ptr) != 4 || mxGetN(ptr) != 4)
                mexErrMsgTxt("Affine transform matrix must be 4x4.");
	p = (double *)mxGetPr(ptr);

	for(i=0; i<3; i++)
		for(j=0; j<4; j++)
			M[j][i] = p[i+4*j];

	if (p[3+4*0] != 0.0 || p[3+4*1] != 0.0 || p[3+4*2] != 0.0 || p[3+4*3] != 1.0)
		mexErrMsgTxt("No perspective projections allowed.");
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	float  *y0=0,  *y1=0,  *y2=0, *iy0=0, *iy1=0, *iy2=0;
	int dim_g[3], dim_f[3];
	REAL U[4][3], V[4][3];

        if (nrhs != 6 || nlhs > 3) mexErrMsgTxt("Incorrect usage.");

	y0 = get_volume(prhs[0], dim_g);
	y1 = get_volume(prhs[1], dim_f);
	if (dim_g[0] != dim_f[0] || dim_g[1] != dim_f[1] || dim_g[2] != dim_f[2])
		mexErrMsgTxt("Incompatible dimensions.");
	y2 = get_volume(prhs[2], dim_f);
	if (dim_g[0] != dim_f[0] || dim_g[1] != dim_f[1] || dim_g[2] != dim_f[2])
		mexErrMsgTxt("Incompatible dimensions.");

	if (!mxIsNumeric(prhs[3]) || mxIsComplex(prhs[3]) ||
		mxIsComplex(prhs[3]) || !mxIsDouble(prhs[3]) || mxGetM(prhs[3]) * mxGetN(prhs[3]) != 3)
		mexErrMsgTxt("Output dimensions must be numeric, real, full, double and contain three elements.");

	dim_f[0] = mxGetPr(prhs[3])[0];
	dim_f[1] = mxGetPr(prhs[3])[1];
	dim_f[2] = mxGetPr(prhs[3])[2];

	get_mat(prhs[4],U);
	get_mat(prhs[5],V);

	plhs[0] = mxCreateNumericArray(3, dim_f,mxSINGLE_CLASS,mxREAL);
	plhs[1] = mxCreateNumericArray(3, dim_f,mxSINGLE_CLASS,mxREAL);
	plhs[2] = mxCreateNumericArray(3, dim_f,mxSINGLE_CLASS,mxREAL);

	iy0 = (float *)mxGetPr(plhs[0]);
	iy1 = (float *)mxGetPr(plhs[1]);
	iy2 = (float *)mxGetPr(plhs[2]);

	setnan(iy0, dim_f[0]*dim_f[1]*dim_f[2]);
	setnan(iy1, dim_f[0]*dim_f[1]*dim_f[2]);
	setnan(iy2, dim_f[0]*dim_f[1]*dim_f[2]);

	invert_field(dim_g, y0, y1, y2, dim_f, iy0, iy1, iy2, U, V);
}
