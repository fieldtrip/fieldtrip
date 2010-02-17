/*
En entree : ind_tr, XYZv, weight, defl
En Sortie : Cii_cog

weight=(sig[0]-sig[1])/((sig[0]+sig[1])*2.*pi) ;

*/
#include <math.h>
#include "mex.h"

void Cii_cog(double C[], double ind_tr[], int Ntri, double XYZv[], int Nvert,
	double weight, double defl)
{
	double dO ;
	double si[3], su[3], sv[3], sw[3] ;
	double r1[3], r2[3], r3[3], det, den, d1, d2, d3 ;
	int i, j, u, v, w ;

	printf("Ntri = %d\n",Ntri);
	printf("weight = %f , defl =  %f \n",weight,defl) ;

	for (i=0;i<Ntri;++i) {
		u=ind_tr[i] - 1 ;
		v=ind_tr[i+Ntri] - 1 ;
		w=ind_tr[i+2*Ntri] - 1 ;
		si[0]=(XYZv[u] + XYZv[v] + XYZv[w]) / 3 ;
		si[1]=(XYZv[u+Nvert] + XYZv[v+Nvert] + XYZv[w+Nvert]) / 3 ;
		si[2]=(XYZv[u+2*Nvert] + XYZv[v+2*Nvert] +
			 XYZv[w+2*Nvert]) / 3 ;

		for (j=0;j<Ntri;++j) {
			u=ind_tr[j] - 1 ;
			v=ind_tr[j+Ntri] - 1 ;
			w=ind_tr[j+2*Ntri] - 1 ;

			r1[0]=XYZv[u]-si[0]; r1[1]=XYZv[u+Nvert]-si[1];
				r1[2]=XYZv[u+2*Nvert]-si[2];
			r2[0]=XYZv[v]-si[0]; r2[1]=XYZv[v+Nvert]-si[1];
				r2[2]=XYZv[v+2*Nvert]-si[2];
			r3[0]=XYZv[w]-si[0]; r3[1]=XYZv[w+Nvert]-si[1];
				r3[2]=XYZv[w+2*Nvert]-si[2];

			det=r1[0]*r2[1]*r3[2]+r2[0]*r3[1]*r1[2]
				+r3[0]*r1[1]*r2[2]
				-(r1[2]*r2[1]*r3[0]+r1[1]*r2[0]*r3[2]
				+r1[0]*r2[2]*r3[1]) ;
			d1=sqrt(r1[0]*r1[0]+r1[1]*r1[1]+r1[2]*r1[2]) ;
			d2=sqrt(r2[0]*r2[0]+r2[1]*r2[1]+r2[2]*r2[2]) ;
			d3=sqrt(r3[0]*r3[0]+r3[1]*r3[1]+r3[2]*r3[2]) ;
	
			den=d1*d2*d3+(r1[0]*r2[0]+r1[1]*r2[1]+r1[2]*r2[2])*d3
				+(r1[0]*r3[0]+r1[1]*r3[1]+r1[2]*r3[2])*d2
				+(r2[0]*r3[0]+r2[1]*r3[1]+r2[2]*r3[2])*d1 ;
			dO=2.*atan2(det,den) ;

			C[i+j*Ntri] = dO * weight ;
		}
		/* negative identity matrix on the diagonal*/
		C[i+i*Ntri] = -1 ;
	}
	/* Deflation */
	for (i=0;i<(Ntri*Ntri);++i) {
		C[i] = C[i] - defl ;
	}
}

/***************************************************************************/

void mexFunction(
	int nlhs, mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	double *Cii ;
	double *ind_tr, *XYZv, weight, defl ;
	unsigned int mind, nind, mXYZv, nXYZv ;

	if (nrhs!=4) {
		mexErrMsgTxt("4 Inputs : ind_tri, XYZv, weight, defl") ;
	} else if (nlhs>1) {
		mexErrMsgTxt("1 Output : Cii_cog") ;
	}

	mind=mxGetM(prhs[0]) ;
	nind=mxGetN(prhs[0]) ;
	if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) 
		|| mxIsSparse(prhs[0]) || !(nind==3)) {
		mexErrMsgTxt("tri must be a matrix : Ntri x 3") ;
	}

	mXYZv=mxGetM(prhs[1]) ;
	nXYZv=mxGetN(prhs[1]) ;
	if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) 
		|| mxIsSparse(prhs[1]) || (nXYZv!=3) ) {
		mexErrMsgTxt("X must be a matrix : Nvert x 3") ;
	}

	if (!mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2])
		|| mxIsSparse(prhs[2]) || !mxIsDouble(prhs[2])
		|| mxGetN(prhs[2])*mxGetM(prhs[2])!=1 ) {
		mexErrMsgTxt("weight must be a scalar") ;
	}

	if (!mxIsNumeric(prhs[3]) || mxIsComplex(prhs[3])
		|| mxIsSparse(prhs[3]) || !mxIsDouble(prhs[3])
		|| mxGetN(prhs[3])*mxGetM(prhs[3])!=1 ) {
		mexErrMsgTxt("defl must be a scalar") ;
	}

	plhs[0]=mxCreateDoubleMatrix(mind,mind,mxREAL) ;

	Cii=mxGetPr(plhs[0]) ;
	ind_tr=mxGetPr(prhs[0]) ;
	XYZv=mxGetPr(prhs[1]) ;
	weight = mxGetScalar(prhs[2]) ;
	defl = mxGetScalar(prhs[3]) ;

	Cii_cog(Cii, ind_tr, mind, XYZv, mXYZv, weight, defl) ;

	mxSetPr(plhs[0],Cii) ;
}

