/*
Direct potential for a vertex approach

IN : 	XYZ_v, XYZ_d, weight

OUT : Gi *transposed* !!!!


weight = 1/(2*pi*(siga+sigb))
*/
#include <math.h>
#include "mex.h"

void G_vert(double Gi[], double XYZ_v[], int Nvert, double XYZ_d[], 
	int Ndip, double weight)
{
	double ndiff, w2 ;
	double ri[3], diff[3] ;
	int i, j ;


	printf("Nvert = %d\n",Nvert);
	printf("Ndip  = %d\n",Ndip);

	for (i=0;i<Nvert;++i) { 
		/*printf("j= %d\n",j);*/
		ri[0]=XYZ_v[i] ;
		ri[1]=XYZ_v[i+Nvert] ;
		ri[2]=XYZ_v[i+2*Nvert] ;

		for (j=0;j<Ndip;++j) { 
			diff[0]=ri[0]-XYZ_d[j] ;
			diff[1]=ri[1]-XYZ_d[j+Ndip] ;
			diff[2]=ri[2]-XYZ_d[j+2*Ndip] ;

			ndiff=sqrt(diff[0]*diff[0]+diff[1]*diff[1]
				+diff[2]*diff[2]) ;
			w2=weight/(ndiff*ndiff*ndiff) ;

			Gi[i*3*Ndip+3*j]=w2*diff[0] ;
			Gi[i*3*Ndip+3*j+1]=w2*diff[1] ;
			Gi[i*3*Ndip+3*j+2]=w2*diff[2] ;
		}
	}
}

void mexFunction(
	int nlhs, mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	double *Gi ;
	double *XYZ_v, *XYZ_d, weight ;
	unsigned int mXYZv, nXYZv, mXYZd, nXYZd ;

	if (nrhs!=3) {
		mexErrMsgTxt("3 inputs : XYZ_v, XYZ_d, weight") ;
	} else if (nlhs>1) {
		mexErrMsgTxt("1 output : Gi transposed") ;
	}

	mXYZv=mxGetM(prhs[0]) ;
	nXYZv=mxGetN(prhs[0]) ;
	if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) 
		|| mxIsSparse(prhs[0]) || !(nXYZv==3)) {
		mexErrMsgTxt("XYZ_v must be a matrix : Nvert x 3") ;
	}

	mXYZd=mxGetM(prhs[1]) ;
	nXYZd=mxGetN(prhs[1]) ;
	if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) 
		|| mxIsSparse(prhs[1]) || !(nXYZd==3)) {
		mexErrMsgTxt("X_d must be a matrix : Ndip x 3") ;
	}


	if (!mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2])
		|| mxIsSparse(prhs[2]) || !mxIsDouble(prhs[2])
		|| mxGetN(prhs[2])*mxGetM(prhs[2])!=1 ) {
		mexErrMsgTxt("weight must be a scalar") ;
	}


	plhs[0]=mxCreateDoubleMatrix(3*mXYZd,mXYZv,mxREAL) ;

	Gi=mxGetPr(plhs[0]) ;
	XYZ_v=mxGetPr(prhs[0]) ;
	XYZ_d=mxGetPr(prhs[1]) ;
	weight = mxGetScalar(prhs[2]) ;

	G_vert(Gi, XYZ_v, mXYZv, XYZ_d, mXYZd, weight) ;

	mxSetPr(plhs[0],Gi) ;
}
