#ifndef lint
static char sccsid[] = "@(#)spm_brainwarp.c	2.12 (c) John Ashburner MRCCU/FIL (& Matthew Brett MRCCU/CBU) 01/03/14";
#endif /* lint */

#include <math.h>
#include "spm_sys_deps.h"
#include "spm_mapping.h"
extern int AbackslashB(double *, double *, double *);
extern void MtimesX(double *, double *, double *);


/*
INPUTS
T[3*nz*ny*nx + ni*4] - current transform
VF              - image to normalize
ni              - number of templates
vols2[ni]       - templates

nx              - number of basis functions in x
B0[dim1[0]*nx]  - basis functions in x
dB0[dim1[0]*nx] - derivatives of basis functions in x
ny              - number of basis functions in y 
B1[dim1[1]*ny]  - basis functions in y
dB1[dim1[1]*ny] - derivatives of basis functions in y
nz              - number of basis functions in z
B2[dim1[2]*nz]  - basis functions in z
dB2[dim1[2]*nz] - derivatives of basis functions in z

M[4*4]          - transformation matrix
samp[3]         - frequency of sampling template.


OUTPUTS
alpha[(3*nz*ny*nx+ni*4)^2] - A'A
 beta[(3*nz*ny*nx+ni*4)]   - A'b
pss[1*1]                   - sum of squares difference
pnsamp[1*1]                - number of voxels sampled
ss_deriv[3]                - sum of squares of derivatives of residuals
*/

static void matmul(int m, int l, int n, double B[], double C[], double A[])
{
	int i, j, k;
	double tmp;
	for(i=0; i<m; i++)
		for(j=0; j<n; j++)
		{
			tmp = 0.0;
			for(k=0; k<l; k++)
				tmp += B[i+m*k]*C[k+l*j];
			A[i+m*j] = tmp;
		}
}

static void transform_grads(double M[16], double x[3])
{
	double x0, x1;
	x0   = M[0+4*0]*x[0] + M[1+4*0]*x[1] + M[2+4*0]*x[2];
	x1   = M[0+4*1]*x[0] + M[1+4*1]*x[1] + M[2+4*1]*x[2];
	x[2] = M[0+4*2]*x[0] + M[1+4*2]*x[1] + M[2+4*2]*x[2];
	x[1] = x1;
	x[0] = x0;
}

static void mrqcof(double T[], double alpha[], double beta[], double pss[],
	MAPTYPE *VF, int ni, MAPTYPE *VG,
	int nx, double B0[], double dB0[],
	int ny, double B1[], double dB1[],
	int nz, double B2[], double dB2[],
	double M[16], int samp[], int edgeskip[],
	double *pnsamp, double ss_deriv[3],
	MAPTYPE *VWG, MAPTYPE *VWF )
{
	int i1,i2, s0[3], x1,x2, y1,y2, z1,z2, m1, m2, ni4;
	double *dvdt, s2[3], *ptr1, *ptr2, *Tz, *Ty, tmp,
		*betaxy, *betax, *alphaxy, *alphax, ss=0.0,  *scal;
	double *Jz[3][3], *Jy[3][3], J[3][3];
	double *bz3[3], *by3[3], *bx3[3];
	double nsamp = 0.0;
	int *dim1;
	double MW[16];

	if (VWF != (MAPTYPE *)0)
	{
		double MW1[16];
		if (AbackslashB(VWF->mat, VF->mat, MW1))
			mexErrMsgTxt("Can't invert matrix");
		matmul(4,4,4, M, MW1, MW);
	}

	dim1 = VG[0].dim;

	bx3[0] = dB0;	bx3[1] =  B0;	bx3[2] =  B0;
	by3[0] =  B1;	by3[1] = dB1;	by3[2] =  B1;
	bz3[0] =  B2;	bz3[1] =  B2;	bz3[2] = dB2;

	ni4 = ni*4;

	/* rate of change of voxel with respect to change in parameters */
	dvdt    = (double *)mxCalloc( 3*nx    + ni4                 , sizeof(double));

	/* Intermediate storage arrays */
	Tz      = (double *)mxCalloc( 3*nx*ny                       , sizeof(double));
	Ty      = (double *)mxCalloc( 3*nx                          , sizeof(double));
	betax   = (double *)mxCalloc( 3*nx    + ni4                 , sizeof(double));
	betaxy  = (double *)mxCalloc( 3*nx*ny + ni4                 , sizeof(double));
	alphax  = (double *)mxCalloc((3*nx    + ni4)*(3*nx    + ni4), sizeof(double));
	alphaxy = (double *)mxCalloc((3*nx*ny + ni4)*(3*nx*ny + ni4), sizeof(double));

	for (i1=0; i1<3; i1++)
	{
		for(i2=0; i2<3; i2++)
		{
			Jz[i1][i2] = (double *)mxCalloc(nx*ny, sizeof(double));
			Jy[i1][i2] = (double *)mxCalloc(nx   , sizeof(double));
		}
	}

	/* pointer to scales for each of the template images */
	scal = T + 3*nx*ny*nz;

	/* only zero half the matrix */
	m1 = 3*nx*ny*nz+ni4;
	for (x1=0;x1<m1;x1++)
	{
		for (x2=0;x2<=x1;x2++)
			alpha[m1*x1+x2] = 0.0;
		beta[x1]= 0.0;
	}
	ss=ss_deriv[0]=ss_deriv[1]=ss_deriv[2]=0;

	for(s0[2]=1; s0[2]<dim1[2]; s0[2]+=samp[2]) /* For each plane of the template images */
	{
		/* build up the deformation field (and derivatives) from it's seperable form */
		for(i1=0, ptr1=T; i1<3; i1++, ptr1 += nz*ny*nx)
			for(x1=0; x1<nx*ny; x1++)
			{
				/* intermediate step in computing nonlinear deformation field */
				tmp = 0.0;
				for(z1=0; z1<nz; z1++)
					tmp  += ptr1[x1+z1*ny*nx] * B2[dim1[2]*z1+s0[2]];
				Tz[ny*nx*i1 + x1] = tmp;

				/* intermediate step in computing Jacobian of nonlinear deformation field */
				for(i2=0; i2<3; i2++)
				{
					tmp = 0;
					for(z1=0; z1<nz; z1++)
						tmp += ptr1[x1+z1*ny*nx] * bz3[i2][dim1[2]*z1+s0[2]];
					Jz[i2][i1][x1] = tmp;
				}
			}

		/* only zero half the matrix */
		m1 = 3*nx*ny+ni4;
		for (x1=0;x1<m1;x1++)
		{
			for (x2=0;x2<=x1;x2++)
				alphaxy[m1*x1+x2] = 0.0;
			betaxy[x1]= 0.0;
		}

		for(s0[1]=1; s0[1]<dim1[1]; s0[1]+=samp[1]) /* For each row of the template images plane */
		{
			/* build up the deformation field (and derivatives) from it's seperable form */
			for(i1=0, ptr1=Tz; i1<3; i1++, ptr1+=ny*nx)
			{
				for(x1=0; x1<nx; x1++)
				{
					/* intermediate step in computing nonlinear deformation field */
					tmp = 0.0;
					for(y1=0; y1<ny; y1++)
						tmp  += ptr1[x1+y1*nx] *  B1[dim1[1]*y1+s0[1]];
					Ty[nx*i1 + x1] = tmp;

					/* intermediate step in computing Jacobian of nonlinear deformation field */
					for(i2=0; i2<3; i2++)
					{
						tmp = 0;
						for(y1=0; y1<ny; y1++)
							tmp += Jz[i2][i1][x1+y1*nx] * by3[i2][dim1[1]*y1+s0[1]];
						Jy[i2][i1][x1] = tmp;
					}
				}
			}

			/* only zero half the matrix */
			m1 = 3*nx+ni4;
			for(x1=0;x1<m1;x1++)
			{
				for (x2=0;x2<=x1;x2++)
					alphax[m1*x1+x2] = 0.0;
				betax[x1]= 0.0;
			}

			for(s0[0]=1; s0[0]<dim1[0]; s0[0]+=samp[0]) /* For each pixel in the row */
			{
				double trans[3];

				/* nonlinear deformation of the template space, followed by the affine transform */
				for(i1=0, ptr1 = Ty; i1<3; i1++, ptr1 += nx)
				{
					/* compute nonlinear deformation field */
					tmp = 0.0;
					for(x1=0; x1<nx; x1++)
						tmp  += ptr1[x1] * B0[dim1[0]*x1+s0[0]];
					trans[i1] = tmp + s0[i1];

					/* compute Jacobian of nonlinear deformation field */
					for(i2=0; i2<3; i2++)
					{
						if (i1 == i2) tmp = 1.0; else tmp = 0;
						for(x1=0; x1<nx; x1++)
							tmp += Jy[i2][i1][x1] * bx3[i2][dim1[0]*x1+s0[0]];
						J[i2][i1] = tmp;
					}
				}

				/* Affine component */
				MtimesX(M, trans, s2);

				/* is the transformed position in range? */
				if (	s2[0]>=1+edgeskip[0] && s2[0]<VF->dim[0]-edgeskip[0] &&
					s2[1]>=1+edgeskip[1] && s2[1]<VF->dim[1]-edgeskip[1] &&
					s2[2]>=1+edgeskip[2] && s2[2]<VF->dim[2]-edgeskip[2] )
				{
					double f, df[3], dv, dvds0[3];
					double wtf, wtg, wt;
					double s0d[3];
					s0d[0]=s0[0];s0d[1]=s0[1];s0d[2]=s0[2];

					resample_d(1,VF,&f,&df[0],&df[1],&df[2],s2,s2+1,s2+2, 1, 0.0);
					transform_grads(M, df);

					if (VWG != (MAPTYPE *)0) resample(1,VWG,&wtg,s0d,s0d+1,s0d+2, 1, 0.0);
					else wtg = 1.0;

					if (VWF != (MAPTYPE *)0)
					{
						double s3[3];
						MtimesX(MW, trans, s3);
						resample(1,VWF,&wtf,s3,s3+1,s3+2, 1, 0.0);
					}
					else wtf = 1.0;

					if (wtf && wtg) wt = sqrt(1.0 /(1.0/wtf + 1.0/wtg));
					else wt = 0.0;

					/* nonlinear transform the gradients to the same space as the template */
					dvds0[0] = J[0][0]*df[0] + J[0][1]*df[1] + J[0][2]*df[2];
					dvds0[1] = J[1][0]*df[0] + J[1][1]*df[1] + J[1][2]*df[2];
					dvds0[2] = J[2][0]*df[0] + J[2][1]*df[1] + J[2][2]*df[2];

					dv = f;
					for(i1=0; i1<ni; i1++)
					{
						double g, dg[3], tmp;
						resample_d(1,&VG[i1],&g,dg,dg+1,dg+2,s0d,s0d+1,s0d+2, 1, 0.0);

						/* linear combination of image and image modulated by constant
						   gradients in x, y and z */
						dvdt[i1*4  +3*nx] = wt*g;
						dvdt[i1*4+1+3*nx] = wt*g*s2[0];
						dvdt[i1*4+2+3*nx] = wt*g*s2[1];
						dvdt[i1*4+3+3*nx] = wt*g*s2[2];

						tmp = scal[i1*4] + s2[0]*scal[i1*4+1] +
							s2[1]*scal[i1*4+2] + s2[2]*scal[i1*4+3];

						dv       -= tmp*g;
						dvds0[0] -= tmp*dg[0];
						dvds0[1] -= tmp*dg[1];
						dvds0[2] -= tmp*dg[2];
					}

					for(i1=0; i1<3; i1++)
					{
						double tmp = -wt*df[i1];
						for(x1=0; x1<nx; x1++)
							dvdt[i1*nx+x1] = tmp * B0[dim1[0]*x1+s0[0]];
					}

					/* cf Numerical Recipies "mrqcof.c" routine */
					m1 = 3*nx+ni4;
					for(x1=0; x1<m1; x1++)
					{
						for (x2=0;x2<=x1;x2++)
							alphax[m1*x1+x2] += dvdt[x1]*dvdt[x2];
						betax[x1] += dvdt[x1]*dv*wt;
					}

					/* sum of squares */
					wt          *= wt;
					nsamp       += wt;
					ss          += wt*dv*dv;
					ss_deriv[0] += wt*dvds0[0]*dvds0[0];
					ss_deriv[1] += wt*dvds0[1]*dvds0[1];
					ss_deriv[2] += wt*dvds0[2]*dvds0[2];
				}
			}

			m1 = 3*nx*ny+ni4;
			m2 = 3*nx+ni4;

			/* Kronecker tensor products */
			for(y1=0; y1<ny; y1++)
			{
				double wt1 = B1[dim1[1]*y1+s0[1]];

				for(i1=0; i1<3; i1++)	/* loop over deformations in x, y and z */
				{
					/* spatial-spatial covariances */
					for(i2=0; i2<=i1; i2++)	/* symmetric matrixes - so only work on half */
					{
						for(y2=0; y2<=y1; y2++)
						{
							/* Kronecker tensor products with B1'*B1 */
							double wt2 = wt1 * B1[dim1[1]*y2+s0[1]];

							ptr1 = alphaxy + nx*(m1*(ny*i1 + y1) + ny*i2 + y2);
							ptr2 = alphax  + nx*(m2*i1 + i2);

							for(x1=0; x1<nx; x1++)
							{
								for (x2=0;x2<=x1;x2++)
									ptr1[m1*x1+x2] += wt2 * ptr2[m2*x1+x2];
							}
						}
					}

					/* spatial-intensity covariances */
					ptr1 = alphaxy + nx*(m1*ny*3 + ny*i1 + y1);
					ptr2 = alphax  + nx*(m2*3 + i1);
					for(x1=0; x1<ni4; x1++)
					{
						for (x2=0;x2<nx;x2++)
							ptr1[m1*x1+x2] += wt1 * ptr2[m2*x1+x2];
					}

					/* spatial component of beta */
					for(x1=0; x1<nx; x1++)
						betaxy[x1+nx*(ny*i1 + y1)] += wt1 * betax[x1 + nx*i1];
				}
			}
			ptr1 = alphaxy + nx*(m1*ny*3 + ny*3);
			ptr2 = alphax  + nx*(m2*3 + 3);
			for(x1=0; x1<ni4; x1++)
			{
				/* intensity-intensity covariances  */
				for (x2=0; x2<=x1; x2++)
					ptr1[m1*x1 + x2] += ptr2[m2*x1 + x2];

				/* intensity component of beta */
				betaxy[nx*ny*3 + x1] += betax[nx*3 + x1];
			}
		}

		m1 = 3*nx*ny*nz+ni4;
		m2 = 3*nx*ny+ni4;

		/* Kronecker tensor products */
		for(z1=0; z1<nz; z1++)
		{
			double wt1 = B2[dim1[2]*z1+s0[2]];

			for(i1=0; i1<3; i1++)	/* loop over deformations in x, y and z */
			{
				/* spatial-spatial covariances */
				for(i2=0; i2<=i1; i2++)	/* symmetric matrixes - so only work on half */
				{
					for(z2=0; z2<=z1; z2++)
					{
						/* Kronecker tensor products with B2'*B2 */
						double wt2 = wt1 * B2[dim1[2]*z2+s0[2]];

						ptr1 = alpha   + nx*ny*(m1*(nz*i1 + z1) + nz*i2 + z2);
						ptr2 = alphaxy + nx*ny*(m2*i1 + i2);
						for(y1=0; y1<ny*nx; y1++)
						{
							for (y2=0;y2<=y1;y2++)
								ptr1[m1*y1+y2] += wt2 * ptr2[m2*y1+y2];
						}
					}
				}

				/* spatial-intensity covariances */
				ptr1 = alpha   + nx*ny*(m1*nz*3 + nz*i1 + z1);
				ptr2 = alphaxy + nx*ny*(m2*3 + i1);
				for(y1=0; y1<ni4; y1++)
				{
					for (y2=0;y2<ny*nx;y2++)
						ptr1[m1*y1+y2] += wt1 * ptr2[m2*y1+y2];
				}

				/* spatial component of beta */
				for(y1=0; y1<ny*nx; y1++)
					beta[y1 + nx*ny*(nz*i1 + z1)] += wt1 * betaxy[y1 + nx*ny*i1];
			}
		}

		ptr1 = alpha   + nx*ny*(m1*nz*3 + nz*3);
		ptr2 = alphaxy + nx*ny*(m2*3 + 3);
		for(y1=0; y1<ni4; y1++)
		{
			/* intensity-intensity covariances */
			for(y2=0;y2<=y1;y2++)
				ptr1[m1*y1 + y2] += ptr2[m2*y1 + y2];

			/* intensity component of beta */
			beta[nx*ny*nz*3 + y1] += betaxy[nx*ny*3 + y1];
		}
	}


	/* Fill in the symmetric bits
	   - OK I know some bits are done more than once. */

	m1 = 3*nx*ny*nz+ni4;
	for(i1=0; i1<3; i1++)	
	{
		double *ptrz, *ptry, *ptrx;
		for(i2=0; i2<=i1; i2++)
		{
			ptrz = alpha + nx*ny*nz*(m1*i1 + i2);
			for(z1=0; z1<nz; z1++)
				for(z2=0; z2<=z1; z2++)
				{
					ptry = ptrz + nx*ny*(m1*z1 + z2);
					for(y1=0; y1<ny; y1++)
						for (y2=0;y2<=y1;y2++)
						{
							ptrx = ptry + nx*(m1*y1 + y2);
							for(x1=0; x1<nx; x1++)
								for(x2=0; x2<x1; x2++)
									ptrx[m1*x2+x1] = ptrx[m1*x1+x2];
						}
					for(x1=0; x1<nx*ny; x1++)
						for (x2=0; x2<x1; x2++)
							ptry[m1*x2+x1] = ptry[m1*x1+x2];
				}
			for(x1=0; x1<nx*ny*nz; x1++)
				for (x2=0; x2<x1; x2++)
					ptrz[m1*x2+x1] = ptrz[m1*x1+x2];
		}
	}
	for(x1=0; x1<nx*ny*nz*3+ni4; x1++)
		for (x2=0; x2<x1; x2++)
			alpha[m1*x2+x1] = alpha[m1*x1+x2];

	*pss    = ss;
	*pnsamp = nsamp;

	mxFree((char *)dvdt);
	mxFree((char *)Tz);
	mxFree((char *)Ty);
	mxFree((char *)betax);
	mxFree((char *)betaxy);
	mxFree((char *)alphax);
	mxFree((char *)alphaxy);

	for (i1=0; i1<3; i1++)
	{
		for(i2=0; i2<3; i2++)
		{
			mxFree((char *)Jz[i1][i2]);
			mxFree((char *)Jy[i1][i2]);
		}
	}
}

static void scale(int m, double dat[], double s)
{
	int i;
	for(i=0; i<m; i++)
		dat[i]*=s;
}

#define MYMAX(a,b) (((a)>(b)) ? (b) : (a))

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	MAPTYPE *map1, *map2, *mapw, /** object */ *mapw2;
	int i, nx,ny,nz,ni=1, samp[3], edgeskip[3];
	double *M, *B0, *B1, *B2, *dB0, *dB1, *dB2, *T, fwhm, fwhm2, fwhm3, df, chi2=0.0, ss_deriv[3];
	double pixdim[3], nsamp;

	/* also accept 13th argument - object volume weighting */
	int iW;
        if (((nrhs != 11) && (nrhs != 12) && (nrhs != 13)) || (nlhs > 4))
        {
                mexErrMsgTxt("Inappropriate usage. ([A,B,var,fwhm]=f(V1,V2,M,B0,B1,B2,dB0,dB1,dB2,T,fwhm);)");
        }

        map1 = get_maps(prhs[0], &ni);
        map2 = get_maps(prhs[1], &i);
	if (i!=1)
	{
		free_maps(map1, ni);
		free_maps(map2,  i);
		mexErrMsgTxt("Inappropriate usage.");
	}

	for (i=2; i<=10; i++)
		if (!mxIsNumeric(prhs[i]) || mxIsComplex(prhs[i]) ||
			mxIsSparse(prhs[i]) || !mxIsDouble(prhs[i]))
		{
			free_maps(map1, ni);
			free_maps(map2,  1);
			mexErrMsgTxt("Inputs must be numeric, real, full and double.");
		}

	if (mxGetM(prhs[2]) != 4 || mxGetN(prhs[2]) != 4)
	{
		free_maps(map1, ni);
		free_maps(map2,  1);
		mexErrMsgTxt("Transformation matrix must be 4x4.");
	}
	M = mxGetPr(prhs[2]);


	if ( mxGetM(prhs[3]) != map1[0].dim[0])
	{
		free_maps(map1, ni);
		free_maps(map2,  1);
		mexErrMsgTxt("Wrong sized X basis functions.");
	}

	nx = mxGetN(prhs[3]);
	B0 = mxGetPr(prhs[3]);
	if ( mxGetM(prhs[6]) != map1[0].dim[0] || mxGetN(prhs[6]) != nx)
	{
		free_maps(map1, ni);
		free_maps(map2,  1);
		mexErrMsgTxt("Wrong sized X basis function derivatives.");
	}

	dB0 = mxGetPr(prhs[6]);

	if ( mxGetM(prhs[4]) != map1[0].dim[1])
	{
		free_maps(map1, ni);
		free_maps(map2,  1);
		mexErrMsgTxt("Wrong sized Y basis functions.");
	}

	ny = mxGetN(prhs[4]);
	B1 = mxGetPr(prhs[4]);
	if ( mxGetM(prhs[7]) != map1[0].dim[1] || mxGetN(prhs[7]) != ny)
	{
		free_maps(map1, ni);
		free_maps(map2,  1);
		mexErrMsgTxt("Wrong sized Y basis function derivatives.");
	}

	dB1 = mxGetPr(prhs[7]);

	if ( mxGetM(prhs[5]) != map1[0].dim[2])
	{
		free_maps(map1, ni);
		free_maps(map2,  1);
		mexErrMsgTxt("Wrong sized Z basis functions.");
	}
	nz = mxGetN(prhs[5]);
	B2 = mxGetPr(prhs[5]);
	if ( mxGetM(prhs[8]) != map1[0].dim[2] || mxGetN(prhs[8]) != nz)
	{
		free_maps(map1, ni);
		free_maps(map2,  1);
		mexErrMsgTxt("Wrong sized Z basis function derivatives.");
	}
	dB2 = mxGetPr(prhs[8]);

	T = mxGetPr(prhs[9]);
	if (mxGetM(prhs[9])*mxGetN(prhs[9]) != 3*nx*ny*nz+ni*4)
	{
		free_maps(map1, ni);
		free_maps(map2,  1);
		mexErrMsgTxt("Transform is wrong size.");
	}

	if (mxGetM(prhs[10])*mxGetN(prhs[10]) == 1)
	{
		fwhm  = mxGetPr(prhs[10])[0];
		fwhm2 = mxGetPr(prhs[10])[0];
	}
	else if (mxGetM(prhs[10])*mxGetN(prhs[10]) == 2)
	{
		fwhm  = mxGetPr(prhs[10])[0];
		fwhm2 = mxGetPr(prhs[10])[1];
	}
	else
		mexErrMsgTxt("FWHM should contain one or two values.");

	voxdim(&map2[0],pixdim);
	/* Because of edge effects from the smoothing, ignore voxels that are too close */
	edgeskip[0]   = rint(fwhm/pixdim[0]); edgeskip[0] = ((edgeskip[0]<1) ? 0 : edgeskip[0]);
	edgeskip[1]   = rint(fwhm/pixdim[1]); edgeskip[1] = ((edgeskip[1]<1) ? 0 : edgeskip[1]);
	edgeskip[2]   = rint(fwhm/pixdim[2]); edgeskip[2] = ((edgeskip[2]<1) ? 0 : edgeskip[2]);

	voxdim(&map1[0],pixdim);

	/* sample about every fwhm/2 */
	samp[0]   = rint(fwhm/2.0/pixdim[0]); samp[0] = ((samp[0]<1) ? 1 : samp[0]);
	samp[1]   = rint(fwhm/2.0/pixdim[1]); samp[1] = ((samp[1]<1) ? 1 : samp[1]);
	samp[2]   = rint(fwhm/2.0/pixdim[2]); samp[2] = ((samp[2]<1) ? 1 : samp[2]);

	for(i=0; i<ni; i++)
	{
		if (map1[0].dim[0] != map1[i].dim[0] || map1[0].dim[1] != map1[i].dim[1] || map1[0].dim[2] != map1[i].dim[2])
			mexErrMsgTxt("Volumes must have same dimensions.");
	}

	plhs[0] = mxCreateDoubleMatrix(3*nx*ny*nz+ni*4,3*nx*ny*nz+ni*4,mxREAL);
	plhs[1] = mxCreateDoubleMatrix(3*nx*ny*nz+ni*4,1,mxREAL);
	plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);
	plhs[3] = mxCreateDoubleMatrix(1,1,mxREAL);

	/* initialise weighting images */
	mapw = mapw2 = (MAPTYPE *)0;

	if (nrhs>=12)
	{
		/* No call to map routine if null */
		if (!mxIsEmpty(prhs[11]))
		{
			mapw = get_maps(prhs[11], &i);
			if (i!=1)
			{
				free_maps(map1, ni);
				free_maps(map2,  1);
				free_maps(mapw,  i);
				mexErrMsgTxt("Inappropriate usage.");
			}

		}
		/* Check for object weighting image */
		if (nrhs>=13 && !mxIsEmpty(prhs[12]))
		{
			mapw2 = get_maps(prhs[12], &iW);
			if (iW!=1)
			{
				free_maps(map1, ni);
				free_maps(map2,  1);
				free_maps(mapw2,  iW);
				free_maps(mapw,  i);
				mexErrMsgTxt("Inappropriate usage.");
			}
		}
	}

	mrqcof(T, mxGetPr(plhs[0]), mxGetPr(plhs[1]), &chi2,
		map2, ni,map1,
		nx,B0,dB0, ny,B1,dB1, nz,B2,dB2, M, samp, edgeskip,
		&nsamp, ss_deriv, mapw, mapw2);

	fwhm3 = ((pixdim[0]/sqrt(2.0*ss_deriv[0]/chi2))*sqrt(8.0*log(2.0)) +
	         (pixdim[1]/sqrt(2.0*ss_deriv[1]/chi2))*sqrt(8.0*log(2.0)) + 
	         (pixdim[2]/sqrt(2.0*ss_deriv[2]/chi2))*sqrt(8.0*log(2.0)))/3.0;

	*mxGetPr(plhs[3]) = fwhm3;

	if (fwhm3<fwhm2)
		fwhm2 = fwhm3;
	if (fwhm2<fwhm)
		fwhm2 = fwhm;

	/* W = fwhm/sqrt(8*log(2))
	   W*sqrt(2*pi) = fwhm*1.0645 */
	df = (MYMAX((pixdim[0]*samp[0])/(fwhm2*1.0645),1.0) *
	      MYMAX((pixdim[1]*samp[1])/(fwhm2*1.0645),1.0) *
	      MYMAX((pixdim[2]*samp[2])/(fwhm2*1.0645),1.0)) * (nsamp - (3*nx*ny*nz + ni*4));

	chi2 /= df;
	mxGetPr(plhs[2])[0] = chi2;

	scale((3*nx*ny*nz+ni*4)*(3*nx*ny*nz+ni*4), mxGetPr(plhs[0]), 1.0/chi2);
	scale((3*nx*ny*nz+ni*4)                  , mxGetPr(plhs[1]), 1.0/chi2);

	free_maps(map1, ni);
	free_maps(map2,  1);
	if (mapw2 != (MAPTYPE *)0)
		free_maps(mapw2,  iW);
	if (mapw != (MAPTYPE *)0)
		free_maps(mapw,  i);
}

