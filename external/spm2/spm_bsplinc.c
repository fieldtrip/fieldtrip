#ifndef lint
static char sccsid[]="@(#)spm_bsplinc.c	2.4 John Ashburner 03/05/07";
#endif
/*
 * This code is a modified version of that of Philippe Thevenaz, which I took from:
 *	http://bigwww.epfl.ch/algorithms.html
 *
 * It has been substantially modified, so blame me (John Ashburner) if there
 * are any bugs. Many thanks to Philippe Thevenaz for advice with the code.
 *
 * See:
 *	M. Unser, A. Aldroubi and M. Eden.
 *	"B-Spline Signal Processing: Part I-Theory,"
 *	IEEE Transactions on Signal Processing 41(2):821-832 (1993).
 *
 *	M. Unser, A. Aldroubi and M. Eden.
 *	"B-Spline Signal Processing: Part II-Efficient Design and Applications,"
 *	IEEE Transactions on Signal Processing 41(2):834-848 (1993).
 *
 *	M. Unser.
 *	"Splines: A Perfect Fit for Signal and Image Processing,"
 *	IEEE Signal Processing Magazine 16(6):22-38 (1999).
 *
*/

#include <math.h>
#include <mex.h>
#include "spm_sys_deps.h"
#include "spm_mapping.h"
#include "spm_datatypes.h"


/***************************************************************************************
Starting periodic boundary condition based on Eq. 2.6 of Unser's 2nd 1993 paper.
	c - vector of unfiltered data
	m - length of c
	p - pole (root of polynomial)
	function returns value that c[0] should initially take

The expression for the first pass of the recursive convolution is:
	for (i=1; i<m; i++) c[i] += p*c[i-1];

If m==4, then:
	c0 = c0 + p*c3;
	c1 = c1 + p*c0;
	c2 = c2 + p*c1
	c3 = c3 + p*c2;
	c0 = c0 + p*c3;
	etc...
After recursive substitution, c0 becomes:
	(1  +p^4+p^8 +p^12 ...)*c0 +
	(p  +p^5+p^9 +p^13 ...)*c3 +
	(p^2+p^6+p^10+p^14 ...)*c2 +
	(p^3+p^7+p^11+p^15 ...)*c1

Using maple...
	sum('p^(k*m+n)','k'=0..infinity)
These series converge to...
	(p^n)/(1-p^m)

So c0 becomes:
	(c0 + c3*p + c2*p^2 + c1*p^3)/(1-p^4)
*/
static double cc_wrap(double c[], int m, double p)
{
	double s, pi;
	int    i, m1;

	m1 = ceil(-30/log(fabs(p)));
	if (m1>m) m1=m;

	pi   = p;
	s    = c[0];
	for (i=1; i<m1; i++)
	{
		s   += pi*c[m-i];
		pi  *= p;
	}
	return(s/(1.0-pi));
}

/***************************************************************************************
Starting mirrored boundary condition based on Eq. 2.6 of Unser's 2nd 1993 paper.
	c - vector of unfiltered data
	m - length of c
	p - pole (root of polynomial)
	function returns value that c[0] should initially take
*/
static double cc_mirror(double c[], int m, double p)
{
	double s, pi, p2i, ip;
	int    i, m1;

	m1 = ceil(-30/log(fabs(p)));
	if (m1 < m)
	{
		pi = p;
		s  = c[0];
		for (i=1; i<m1; i++)
		{
			s  += pi * c[i];
			pi *= p;
		}
		return(s);
	}
	else
	{
		pi   = p;
		ip   = 1.0/p;
		p2i  = pow(p,m-1.0);
		s    = c[0] + p2i*c[m-1];
		p2i *= p2i * ip;
		for (i=1; i<m-1; i++)
		{
			s   += (pi+p2i)*c[i];
			pi  *= p;
			p2i *= ip;
		}
		return(s/(1.0-pi*pi));
	}
}

/***************************************************************************************
End periodic boundary condition
	c - first pass filtered data
	m - length of filtered data (must be > 1)
	p - pole
	function returns value for c[m-1] before 2nd filter pass

The expression for the second pass of the recursive convolution is:
	for (i=m-2; i>=0; i--) c[i] = p*(c[i+1]-c[i]);
If m==4, then:
	c3 = p*(c0-c3);
	c2 = p*(c3-c2);
	c1 = p*(c2-c1);
	c0 = p*(c1-c0);
	c3 = p*(c0-c3);
	etc...

After recursive substitution, c0 becomes:
	-(p  +p^5+p^9  ...)*c3
	-(p^2+p^6+p^10 ...)*c0
	-(p^3+p^7+p^11 ...)*c1
	-(p^4+p^8+p^12 ...)*c2

These series converge to...
	(p^n)/(p^m-1)

So c0 becomes:
	(c3*p + c0*p^2 + c1*p^3 + c2*p^4)/(p^4-1)
*/
static double icc_wrap(double c[],int m, double p)
{
	double s, pi;
	int    i, m1;

	m1 = ceil(-30/log(fabs(p)));
	if (m1>m) m1=m;

	pi = p;
	s  = pi*c[m-1];
	for (i=0; i<m1-1; i++)
	{
		pi  *= p;
		s   += pi*c[i];
	}
	return(s/(pi-1.0));
}

/***************************************************************************************
End mirrored boundary condition
	c - first pass filtered data
	m - length of filtered data (must be > 1)
	p - pole
	function returns value for c[m-1] before 2nd filter pass
*/
static double icc_mirror(double c[],int m, double p)
{
	return((p/(p*p-1.0))*(p*c[m-2]+c[m-1]));
}

/***************************************************************************************
Compute gains required for zero-pole representation - see tf2zp.m in Matlab's
 Signal Processing Toolbox.
	p - poles
	np - number of poles
	function returns the gain of the system
*/
static double gain(double p[], int np)
{
	int j;
	double lambda = 1.0;
	for (j = 0; j < np; j++)
		lambda = lambda*(1.0-p[j])*(1.0-1.0/p[j]);
	return(lambda);
}

/***************************************************************************************
One dimensional recursive filtering - assuming wrapped boundaries
See Eq. 2.5 of Unsers 2nd 1993 paper.
	c - original vector on input, coefficients on output
	m - length of vector
	p - poles (polynomial roots)
	np - number of poles
*/
static void splinc_wrap(double c[], int m, double p[], int np)
{
	double lambda = 1.0;
	int i, k;

	if (m == 1) return;

	/* compute gain and apply it */
	lambda = gain(p,np);
	for (i = 0; i < m; i++)
		c[i] *= lambda;

	/* loop over poles */
	for (k = 0; k < np; k++)
	{
		double pp = p[k];
		c[0] = cc_wrap(c, m, pp);

		for (i=1; i<m; i++)
			c[i] += pp*c[i-1];

		c[m-1] = icc_wrap(c, m, pp);
		for (i=m-2; i>=0; i--)
			c[i] = pp*(c[i+1]-c[i]);
	}
}

/***************************************************************************************
One dimensional recursive filtering - assuming mirror boundaries
See Eq. 2.5 of Unsers 2nd 1993 paper.
	c - original vector on input, coefficients on output
	m - length of vector
	p - poles (polynomial roots)
	np - number of poles
*/
static void splinc_mirror(double c[], int m, double p[], int np)
{
	double lambda = 1.0;
	int i, k;

	if (m == 1) return;

	/* compute gain and apply it */
	lambda = gain(p,np);
	for (i = 0; i < m; i++)
		c[i] *= lambda;

	/* loop over poles */
	for (k = 0; k < np; k++)
	{
		double pp = p[k];
		c[0] = cc_mirror(c, m, pp);

		for (i=1; i<m; i++)
			c[i] += pp*c[i-1];

		c[m-1] = icc_mirror(c, m, pp);
		for (i=m-2; i>=0; i--)
			c[i] = pp*(c[i+1]-c[i]);
	}
}

/***************************************************************************************
Return roots of B-spline kernels.
	 d - degree of B-spline
	 np - number of roots of magnitude less than one
	 p - roots.
*/
static int get_poles(int d, int *np, double p[])
{
	/* Return polynomial roots that are less than one. */
	switch (d) {
		case 0:
			*np = 0;
			break;
		case 1:
			*np = 0;
			break;
		case 2: /* roots([1 6 1]) */
			*np = 1;
			p[0] = sqrt(8.0) - 3.0;
			break;
		case 3: /* roots([1 4 1]) */
			*np = 1;
			p[0] = sqrt(3.0) - 2.0;
			break;
		case 4: /* roots([1 76 230 76 1]) */
			*np = 2;
			p[0] = sqrt(664.0 - sqrt(438976.0)) + sqrt(304.0) - 19.0;
			p[1] = sqrt(664.0 + sqrt(438976.0)) - sqrt(304.0) - 19.0;
			break;
		case 5: /* roots([1 26 66 26 1]) */
			*np   = 2;
			p[0] = sqrt(67.5 - sqrt(4436.25)) + sqrt(26.25) - 6.5;
			p[1] = sqrt(67.5 + sqrt(4436.25)) - sqrt(26.25) - 6.5;
			break;
		case 6: /* roots([1 722 10543 23548 10543 722 1]) */
			*np   = 3;
			p[0] = -0.488294589303044755130118038883789062112279161239377608394;
			p[1] = -0.081679271076237512597937765737059080653379610398148178525368;
			p[2] = -0.00141415180832581775108724397655859252786416905534669851652709;
			break;
		case 7: /* roots([1 120 1191 2416 1191 120 1]) */
			*np   = 3;
			p[0] = -0.5352804307964381655424037816816460718339231523426924148812;
			p[1] = -0.122554615192326690515272264359357343605486549427295558490763;
			p[2] = -0.0091486948096082769285930216516478534156925639545994482648003;
			break;
		default:
			return(1);
	}
	return(0);
}


/***************************************************************************************
Deconvolve the B-spline basis functions from the image volume
	vol - a handle for the volume to deconvolve
	c - the coefficients (arising from the deconvolution)
	d - the spline degree
	splinc0, splinc1, splinc2	- functions for 1D deconvolutions
*/
static int vol_coeffs(MAPTYPE *vol, double c[], int d[], void (*splinc[])())
{
	double	p[4], *cp;
	int	np;
	int	i, j, k, n;
	double f[10240];

	/* Check that dimensions don't exceed size of f */
	if (vol->dim[1]>10240 ||vol->dim[2]>10240)
		return(1);

	/* Do a straight copy */
	cp = c;
	for(k=0; k<vol->dim[2]; k++)
	{
		double dk = k+1;
		for(j=0; j<vol->dim[1]; j++)
		{
			double dj = j+1;
			for(i=0;i<vol->dim[0];i++, cp++)
			{
				double di = i+1;
				resample(1,vol,cp,&di,&dj,&dk,0, 0.0);

				/* Not sure how best to handle NaNs */
				if (!finite(*cp)) *cp = 0.0;
			}
		}
	}

	/* Deconvolve along the fastest dimension (X) */
	if (d[0]>1 && vol->dim[0]>1)
	{
		if (get_poles(d[0], &np, p)) return(1);
		for(k=0; k<vol->dim[2]; k++)
		{
			double dk = k+1;
			for(j=0; j<vol->dim[1]; j++)
			{
				cp = &c[vol->dim[0]*(j+vol->dim[1]*k)];
				splinc[0](cp, vol->dim[0], p, np);
			}
		}
	}

	/* Deconvolve along the middle dimension (Y) */
	if (d[1]>1 && vol->dim[1]>1)
	{
		if (get_poles(d[1], &np, p)) return(1);
		n =vol->dim[0];
		for(k=0; k<vol->dim[2]; k++)
		{
			for(i=0;i<vol->dim[0];i++)
			{
				cp = &c[i+vol->dim[0]*vol->dim[1]*k];
				for(j=0; j<vol->dim[1]; j++, cp+=n)
					f[j] = *cp;
				splinc[1](f, vol->dim[1], p, np);
				cp = &c[i+vol->dim[0]*vol->dim[1]*k];
				for(j=0; j<vol->dim[1]; j++, cp+=n)
					*cp = f[j];
			}
		}
	}

	/* Deconvolve along the slowest dimension (Z) */
	if (d[2]>1 && vol->dim[2]>1)
	{
		if (get_poles(d[2], &np, p)) return(1);
		n = vol->dim[0]*vol->dim[1];
		for(j=0; j<vol->dim[1]; j++)
		{
			for(i=0;i<vol->dim[0];i++)
			{
				cp = &c[i+vol->dim[0]*j];
				for(k=0; k<vol->dim[2]; k++, cp+=n)
					f[k] = *cp;
				splinc[2](f, vol->dim[2], p, np);
				cp = &c[i+vol->dim[0]*j];
				for(k=0; k<vol->dim[2]; k++, cp+=n)
					*cp = f[k];
			}
		}
	}
	return(0);
}

/***************************************************************************************
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int k, d[3], sts;
	MAPTYPE *vol, *get_maps();
	double *c;
	void (*splinc[3])();

	if (nrhs < 2 || nlhs > 1)
		mexErrMsgTxt("Inappropriate usage.");
	if (mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) ||
		(mxGetM(prhs[1])*mxGetN(prhs[1]) != 3 && mxGetM(prhs[1])*mxGetN(prhs[1]) != 6))
		mexErrMsgTxt("Inappropriate usage.");

	for(k=0; k<3; k++)
	{
		d[k] = floor(mxGetPr(prhs[1])[k]+0.5);
		if (d[k]<0 || d[k]>7)
			mexErrMsgTxt("Bad spline degree.");
	}

	for(k=0; k<3; k++) splinc[k] = splinc_mirror;
	if (mxGetM(prhs[1])*mxGetN(prhs[1]) == 6)
	{
		for(k=0; k<3; k++)
			if (mxGetPr(prhs[1])[k+3])
				splinc[k] = splinc_wrap;
	}

	vol=get_maps(prhs[0], &k);
	if (k!=1)
	{
		free_maps(vol, k);
		mexErrMsgTxt("Too many images.");
	}

	plhs[0] = mxCreateNumericArray(3,vol->dim, mxDOUBLE_CLASS, mxREAL);
	c = mxGetPr(plhs[0]);

	sts = vol_coeffs(vol, c, d, splinc);

	if (sts)
	{
		free_maps(vol, k);
		mexErrMsgTxt("Problem with deconvolution.");
	}
	free_maps(vol, k);
}
