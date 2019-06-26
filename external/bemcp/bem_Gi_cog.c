/*
  Direct potential for a COG approach

  Input : XYZv, ind_tr, XYZd, weight
  Output : Gi *transposed* !!!!

  weight = 1./(2.*pi*(siga[0]+sigb[0]));
*/

#include <math.h>
#include "mex.h"

void Gi_cog(double Gi[],
	double XYZv[], int Nvert, double ind_tr[], int Ntri,
	double XYZd[], int Ndip, double weight)
{
	double ndiff, w2 ;
	double ri[3], diff[3] ;
	int i, j, u, v, w ;

	/*printf("Ndip = %d\n",Ndip);*/
	/*printf("Ntri = %d\n",Ntri);*/

	for (i=0;i<Ntri;++i) { 
		u=ind_tr[i] - 1 ;
		v=ind_tr[i+Ntri] - 1 ;
		w=ind_tr[i+2*Ntri] - 1 ;
		ri[0]=(XYZv[u] + XYZv[v] + XYZv[w]) / 3 ;
		ri[1]=(XYZv[u+Nvert] + XYZv[v+Nvert] + XYZv[w+Nvert]) / 3 ;
		ri[2]=(XYZv[u+2*Nvert] + XYZv[v+2*Nvert] +
			 XYZv[w+2*Nvert]) / 3 ;

		for (j=0;j<Ndip;++j) { 
			diff[0]=ri[0]-XYZd[j] ;
			diff[1]=ri[1]-XYZd[j+Ndip] ;
			diff[2]=ri[2]-XYZd[j+2*Ndip] ;

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
	double *XYZv, *ind_tr, *XYZd, weight ;
	unsigned int mXYZv, nXYZv, mind, nind, mXYZd, nXYZd ;

	if (nrhs!=4) {
		mexErrMsgTxt("4 Inputs : XYZv, ind_tr, XYZd, weight") ;
	} else if (nlhs>1) {
		mexErrMsgTxt("1 Output : Gi transposed") ;
	}

	mXYZv=mxGetM(prhs[0]) ;
	nXYZv=mxGetN(prhs[0]) ;

	if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) 
		|| mxIsSparse(prhs[0]) || (nXYZv!=3)) {
		mexErrMsgTxt("XYZv must be a matrix : Nvert x 3") ;
	}

	mind=mxGetM(prhs[1]) ;
	nind=mxGetN(prhs[1]) ;
	if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) 
		|| mxIsSparse(prhs[1]) || (nind!=3) ) {
		mexErrMsgTxt("ind_tr must be a matrix : Ntri x 3") ;
	}
 
	mXYZd=mxGetM(prhs[2]) ;
	nXYZd=mxGetN(prhs[2]) ;
	if (!mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2]) 
		|| mxIsSparse(prhs[2]) || (nXYZd!=3)) {
		mexErrMsgTxt("XYZd  must be a matrix : Ndip x 3") ;
	}

	if (!mxIsNumeric(prhs[3]) || mxIsComplex(prhs[3])
		|| mxIsSparse(prhs[3]) || !mxIsDouble(prhs[3])
		|| mxGetN(prhs[3])*mxGetM(prhs[3])!=1 ) {
		mexErrMsgTxt("weight must be a scalar") ;
	}


	plhs[0]=mxCreateDoubleMatrix(3*mXYZd,mind,mxREAL) ;

	Gi=mxGetPr(plhs[0]) ;
	XYZv=mxGetPr(prhs[0]) ;
	ind_tr=mxGetPr(prhs[1]) ;
	XYZd=mxGetPr(prhs[2]) ;
	weight = mxGetScalar(prhs[3]) ;

	Gi_cog(Gi,XYZv,mXYZv,ind_tr,mind,XYZd,mXYZd,weight) ;

	mxSetPr(plhs[0],Gi) ;
}
