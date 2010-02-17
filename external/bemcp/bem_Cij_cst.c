/*
En entree :  XYZva, XYZvb, trib, weight, defl
En sortie : Cii


siga==sig[0], sigb==sig[1], sigc==sig[2], sigd==sig[3]
weight=(sig[0]-sig[1])/((sig[2]+sig[3])*2.*pi) ;

*/
#include <math.h>
#include "mex.h"

void Cij_cst(double C[], double XYZva[], int Nverta, double XYZvb[], int Nvertb,
	double trib[], int Ntrib, double weight, double defl)
{
	double dO, dO1, dO2, dO3 ;
	double si[3], r1[3], r2[3], r3[3], det, den, d1, d2, d3 ;
	double r21[3], r32[3], r13[3], l21, l32, l13 ;
	double g1, g2 ,g3, Ome[3], z1[3], z2[3], z3[3] ;
	double tmp[3], n[3], nn, beta, ideuxA ;
	int i, j, u, v, w ;

	printf("Nverta = %d\n",Nverta);
	printf("Ntrib = %d\n",Ntrib);

	for (i=0;i<Nverta;++i) { /* Run through all vert of surf a */ 
		si[0] = XYZva[i] ;
		si[1] = XYZva[i+Nverta] ;
		si[2] = XYZva[i+2*Nverta] ;

		for (j=0;j<Ntrib;++j) { /* Run through all tri of surf b */
			u=trib[j] - 1 ;
			v=trib[j+Ntrib] - 1 ;
			w=trib[j+2*Ntrib] - 1 ;

			r1[0]=XYZvb[u]-si[0];
				r1[1]=XYZvb[u+Nvertb]-si[1];
				r1[2]=XYZvb[u+2*Nvertb]-si[2];
			r2[0]=XYZvb[v]-si[0];
				r2[1]=XYZvb[v+Nvertb]-si[1];
				r2[2]=XYZvb[v+2*Nvertb]-si[2];
			r3[0]=XYZvb[w]-si[0];
				r3[1]=XYZvb[w+Nvertb]-si[1];
				r3[2]=XYZvb[w+2*Nvertb]-si[2];
			d1=sqrt(r1[0]*r1[0]+r1[1]*r1[1]+r1[2]*r1[2]) ;
			d2=sqrt(r2[0]*r2[0]+r2[1]*r2[1]+r2[2]*r2[2]) ;
			d3=sqrt(r3[0]*r3[0]+r3[1]*r3[1]+r3[2]*r3[2]) ;

			det=r1[0]*r2[1]*r3[2]+r2[0]*r3[1]*r1[2]
					+r3[0]*r1[1]*r2[2]
				-(r1[2]*r2[1]*r3[0]+r1[1]*r2[0]*r3[2]
					+r1[0]*r2[2]*r3[1]) ;
	
			den=d1*d2*d3+(r1[0]*r2[0]+r1[1]*r2[1]+r1[2]*r2[2])*d3
				+(r1[0]*r3[0]+r1[1]*r3[1]+r1[2]*r3[2])*d2
				+(r2[0]*r3[0]+r2[1]*r3[1]+r2[2]*r3[2])*d1 ;
			dO=2.*atan2(det,den) ;
			/* Solid angle */

			C[i+u*Nverta] = C[i+u*Nverta] + dO/3*weight ;
 			C[i+v*Nverta] = C[i+v*Nverta] + dO/3*weight ;
 			C[i+w*Nverta] = C[i+w*Nverta] + dO/3*weight ;
		}
	}

	/* Deflation */
	for (i=0;i<(Nverta*Nvertb);++i) {
		C[i] = C[i] - defl ;
	}
}

/****************************************************************************/

void mexFunction(
	int nlhs, mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	double *C ;
	double *XYZva, *XYZvb, *trib, weight, defl ;
	unsigned int mXYZva, nXYZva, mXYZvb, nXYZvb, mtrib, ntrib ;

	if (nrhs!=5) {
		mexErrMsgTxt("5 Inputs : XYZva, XYZvb, trib, weight, Npt.") ;
	} else if (nlhs>1) {
		mexErrMsgTxt("1 Output : Cii") ;
	}

	mXYZva=mxGetM(prhs[0]) ;
	nXYZva=mxGetN(prhs[0]) ;
	if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) 
		|| mxIsSparse(prhs[0]) || (nXYZva!=3)) {
		mexErrMsgTxt("XYZva must be a matrix : Nverta x 3") ;
	}

	mXYZvb=mxGetM(prhs[1]) ;
	nXYZvb=mxGetN(prhs[1]) ;
	if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) 
		|| mxIsSparse(prhs[1]) || (nXYZvb!=3)) {
		mexErrMsgTxt("XYZvb must be a matrix : Nvertb x 3") ;
	}

	mtrib=mxGetM(prhs[2]) ;
	ntrib=mxGetN(prhs[2]) ;
	if (!mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2]) 
		|| mxIsSparse(prhs[2]) || !(ntrib==3)) {
		mexErrMsgTxt("trib must be a matrix : Ntrb x 3") ;
	}
	
	if (!mxIsNumeric(prhs[3]) || mxIsComplex(prhs[3])
		|| mxIsSparse(prhs[3]) || !mxIsDouble(prhs[3])
		|| mxGetN(prhs[3])*mxGetM(prhs[3])!=1 ) {
		mexErrMsgTxt("weight must be a scalar") ;
	}

	if (!mxIsNumeric(prhs[4]) || mxIsComplex(prhs[4])
		|| mxIsSparse(prhs[4]) || !mxIsDouble(prhs[4])
		|| mxGetN(prhs[4])*mxGetM(prhs[4])!=1 ) {
		mexErrMsgTxt("defl must be a scalar") ;
	}

	plhs[0]=mxCreateDoubleMatrix(mXYZva,mXYZvb,mxREAL) ;

	C=mxGetPr(plhs[0]) ;
	XYZva=mxGetPr(prhs[0]) ;
	XYZvb=mxGetPr(prhs[1]) ;
	trib=mxGetPr(prhs[2]) ;
	weight = mxGetScalar(prhs[3]) ;
	defl = mxGetScalar(prhs[4]) ;

	Cij_cst(C, XYZva, mXYZva, XYZvb, mXYZvb, trib, mtrib, weight, defl) ;
	mxSetPr(plhs[0],C) ;
}






