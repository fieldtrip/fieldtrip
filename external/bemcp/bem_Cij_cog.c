/*
  Input : XYZva, ind_tra, XYZvb, ind_trb, weight, defl
  Output : Cij_cog

  weight = (sig[0]-sig[1])/((sig[2]+sig[3])*2.*pi) ;
*/

#include <math.h>
#include "mex.h"

void sol_angle(double *dO, double *r1, double *r2, double *r3)
{
	double det, d1, d2, d3, den ;

	det=r1[0]*r2[1]*r3[2]+r2[0]*r3[1]*r1[2]+r3[0]*r1[1]*r2[2]
		-(r1[2]*r2[1]*r3[0]+r1[1]*r2[0]*r3[2]+r1[0]*r2[2]*r3[1]) ;
	d1=sqrt(r1[0]*r1[0]+r1[1]*r1[1]+r1[2]*r1[2]) ;
	d2=sqrt(r2[0]*r2[0]+r2[1]*r2[1]+r2[2]*r2[2]) ;
	d3=sqrt(r3[0]*r3[0]+r3[1]*r3[1]+r3[2]*r3[2]) ;
	
	den=d1*d2*d3+(r1[0]*r2[0]+r1[1]*r2[1]+r1[2]*r2[2])*d3
		+(r1[0]*r3[0]+r1[1]*r3[1]+r1[2]*r3[2])*d2
		+(r2[0]*r3[0]+r2[1]*r3[1]+r2[2]*r3[2])*d1 ;
	*dO=2.*atan2(det,den) ;
}

void Cij(double C[], double XYZva[], int Nverta,
	double ind_tra[], int Ntria, double XYZvb[], int Nvertb,
	double ind_trb[], int Ntrib, double weight, double defl)
{
	double dO ;
	double si[3] ;
	double r1[3], r2[3], r3[3], det, den, d1, d2, d3 ;
	int i, j, u, v, w ;


	/*printf("Ntria = %d\n",Ntria);*/
	/*printf("Ntrib = %d\n",Ntrib);*/


	for (i=0;i<Ntria;++i) {
		u=ind_tra[i] - 1 ;
		v=ind_tra[i+Ntria] - 1 ;
		w=ind_tra[i+2*Ntria] - 1 ;
		si[0]=(XYZva[u] + XYZva[v] + XYZva[w]) / 3 ;
		si[1]=(XYZva[u+Nverta] + XYZva[v+Nverta] + XYZva[w+Nverta])/3 ;
		si[2]=(XYZva[u+2*Nverta] + XYZva[v+2*Nverta] +
			 XYZva[w+2*Nverta]) / 3 ;

		for (j=0;j<Ntrib;++j) {
			u=ind_trb[j] - 1 ;
			v=ind_trb[j+Ntrib] - 1 ;
			w=ind_trb[j+2*Ntrib] - 1 ;

			r1[0]=XYZvb[u]-si[0]; r1[1]=XYZvb[u+Nvertb]-si[1];
				r1[2]=XYZvb[u+2*Nvertb]-si[2];
			r2[0]=XYZvb[v]-si[0]; r2[1]=XYZvb[v+Nvertb]-si[1];
				r2[2]=XYZvb[v+2*Nvertb]-si[2];
			r3[0]=XYZvb[w]-si[0]; r3[1]=XYZvb[w+Nvertb]-si[1];
				r3[2]=XYZvb[w+2*Nvertb]-si[2];

			sol_angle(&dO,r1,r2,r3) ;

			C[i+j*Ntria] = dO * weight - defl ;
 			C[i+j*Ntria] = dO * weight - defl ;
 			C[i+j*Ntria] = dO * weight - defl ;
		}
	}

}

/****************************************************************************/

void mexFunction(
	int nlhs, mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	double *C ;
	double *XYZva, *ind_tra, *XYZvb, *ind_trb, weight, defl ;
	unsigned int mXYZva, nXYZva, minda, ninda ;
	unsigned int mXYZvb, nXYZvb, mindb, nindb ;

	if (nrhs!=6) {
		mexErrMsgTxt("6 Inputs : XYZva, ind_tra, XYZvb, ind_trb, weight, defl") ;
	} else if (nlhs>1) {
		mexErrMsgTxt("1 Outputs : Cij_cog") ;
	}

	mXYZva=mxGetM(prhs[0]) ;
	nXYZva=mxGetN(prhs[0]) ;
	if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) 
		|| mxIsSparse(prhs[0]) || (nXYZva!=3)) {
		mexErrMsgTxt("XYZva must be a matrix : Nverta x 3") ;
	}

	minda=mxGetM(prhs[1]) ;
	ninda=mxGetN(prhs[1]) ;
	if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) 
		|| mxIsSparse(prhs[1]) || !(ninda==3)) {
		mexErrMsgTxt("ind_tra must be a matrix : Ntria x 3") ;
	}

	mXYZvb=mxGetM(prhs[2]) ;
	nXYZvb=mxGetN(prhs[2]) ;
	if (!mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2]) 
		|| mxIsSparse(prhs[2]) || (nXYZvb!=3)) {
		mexErrMsgTxt("XYZvb must be a matrix : Nvertb") ;
	}

	mindb=mxGetM(prhs[3]) ;
	nindb=mxGetN(prhs[3]) ;
	if (!mxIsNumeric(prhs[3]) || mxIsComplex(prhs[3]) 
		|| mxIsSparse(prhs[3]) || !(nindb==3)) {
		mexErrMsgTxt("ind_trb must be a matrix : Ntrib x 3") ;
	}
	
	if (!mxIsNumeric(prhs[4]) || mxIsComplex(prhs[4])
		|| mxIsSparse(prhs[4]) || !mxIsDouble(prhs[4])
		|| mxGetN(prhs[4])*mxGetM(prhs[4])!=1 ) {
		mexErrMsgTxt("weight must be a scalar") ;
	}

	if (!mxIsNumeric(prhs[5]) || mxIsComplex(prhs[5])
		|| mxIsSparse(prhs[5]) || !mxIsDouble(prhs[5])
		|| mxGetN(prhs[5])*mxGetM(prhs[5])!=1 ) {
		mexErrMsgTxt("defl must be a scalar") ;
	}

	plhs[0]=mxCreateDoubleMatrix(minda,mindb,mxREAL) ;

	C=mxGetPr(plhs[0]) ;
	XYZva=mxGetPr(prhs[0]) ;
	ind_tra=mxGetPr(prhs[1]) ;
	XYZvb=mxGetPr(prhs[2]) ;
	ind_trb=mxGetPr(prhs[3]) ;
	weight = mxGetScalar(prhs[4]) ;
	defl = mxGetScalar(prhs[5]) ;

	Cij(C,XYZva,mXYZva,ind_tra,minda, 
		XYZvb,mXYZvb,ind_trb,mindb,weight,defl);

	mxSetPr(plhs[0],C) ;
}
