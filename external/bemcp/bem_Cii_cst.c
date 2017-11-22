/*
  Input : ind_tr, XYZv, weight, defl, tri_pt4
  Output : Cii_lin

  weight = (sig1-sig2)/((sig1+sig2)*2*pi) ;
*/

#include <math.h>
#include "mex.h"

#define NR_NEIGH_TRI 10 /* Max # of triangle neighbours */

void Cii_constant(double C[], double ind_tr[], int Ntri, double XYZv[],
	int Nvert, double weight, double defl, double tri_pt4[])
{

	double pi, dO, dO1, dO2, dO3 ;
	double si[3], r1[3], r2[3], r3[3], det, den ;
	double d1, d2, d3, sumOm, di2pi ;
	double r21[3], r32[3], r13[3], l21, l32, l13 ;
	double g1, g2 ,g3, Ome[3], z1[3], z2[3], z3[3] ;
	double tmp[3], n[3], nn, beta, ideuxA, epsi ;
	double a[3][3], b[3], deter, rc[3], r_b[3], R ; 
		/* r0 = si, for the calculuus of the ASA */
	double psi1, psi2, phi12, Ome_n[NR_NEIGH_TRI], Ome_ro ;
	int i, j, u, v, w, nze, lze[NR_NEIGH_TRI] ;
		/* 10 neighbouring triangles max for the ASA */

	/*printf("Nvert = %d\n",Nvert) ;*/
	/*printf("Ntri = %d\n",Ntri) ;*/
	/*printf("weight = %f , defl =  %f \n",weight,defl) ;*/

	epsi = .00000000000000001 ;
	pi = 4.*atan(1.) ;

	for (i=0;i<Nvert;++i) { 
	/* Run through all the points */
	/* printf("i = %d\n",i) ;*/

		si[0] = XYZv[i] ;
		si[1] = XYZv[i+Nvert] ;
		si[2] = XYZv[i+2*Nvert] ;
		/* printf("si = %f  %f  %f\n",si[0],si[1],si[2]) ; */

		nze = 0 ; /* # of neighb for ASA */
		sumOm = 0 ; /* total SA for this vert */

		for (j=0;j<Ntri;++j) {
		/* Run through all the triangles */
		/* printf("j = %d\n",j) ; */
			u=ind_tr[j] - 1 ;
			v=ind_tr[j+Ntri] - 1 ;
			w=ind_tr[j+2*Ntri] - 1 ;
			/* triangle vert. indices */
			/* printf("[u v w] = %d %d %d\n",u,v,w) ; */
			if ((i!=u) & (i!=v) & (i!=w)) {
				/* if not an ASA */
				r1[0]=XYZv[u]-si[0];
					r1[1]=XYZv[u+Nvert]-si[1];
					r1[2]=XYZv[u+2*Nvert]-si[2];
				r2[0]=XYZv[v]-si[0];
					r2[1]=XYZv[v+Nvert]-si[1];
					r2[2]=XYZv[v+2*Nvert]-si[2];
				r3[0]=XYZv[w]-si[0];
					r3[1]=XYZv[w+Nvert]-si[1];
					r3[2]=XYZv[w+2*Nvert]-si[2];
				det=r1[0]*r2[1]*r3[2]+r2[0]*r3[1]*r1[2]+
					r3[0]*r1[1]*r2[2]
					-(r1[2]*r2[1]*r3[0]+r1[1]*r2[0]*r3[2]+
					r1[0]*r2[2]*r3[1]) ;
				d1=sqrt(r1[0]*r1[0]+r1[1]*r1[1]+r1[2]*r1[2]) ;
				d2=sqrt(r2[0]*r2[0]+r2[1]*r2[1]+r2[2]*r2[2]) ;
				d3=sqrt(r3[0]*r3[0]+r3[1]*r3[1]+r3[2]*r3[2]) ;
				den=d1*d2*d3+(r1[0]*r2[0]+r1[1]*r2[1]+
					r1[2]*r2[2])*d3
					+(r1[0]*r3[0]+r1[1]*r3[1]+
					r1[2]*r3[2])*d2
					+(r2[0]*r3[0]+r2[1]*r3[1]+
					r2[2]*r3[2])*d1 ;
				dO=2.*atan2(det,den) ;
				/* SA sustended by triangle j and vert i */
				/* printf("dO = %f\n",dO) ; */

				C[i+u*Nvert] = C[i+u*Nvert] + dO/3 * weight ;
 				C[i+v*Nvert] = C[i+v*Nvert] + dO/3 * weight ;
 				C[i+w*Nvert] = C[i+w*Nvert] + dO/3 * weight ;
				/* printf("C(i,[u v w]) = %f  %f  %f\n"
				,C[i+u*Nvert],C[i+v*Nvert],C[i+w*Nvert]) ;
				/* Distribution of the SA on the 3 vert */
				sumOm = sumOm + dO ;
				/* printf("sumOm = %f \n",sumOm) ;*/
				/* Total SA */
			}
			else {
			/* list and # of triangle for autosolid angle */
				lze[nze] = j ;
				nze++ ;
			}
		}
	/*******************************/
	/* Autosolid angle calculation */
		di2pi = (2.*pi-sumOm) ; /* Missing SA */
		/*printf("Missing angle = %f\n",di2pi) ;*/
		Ome_ro = 0 ;
		for (j=0;j<nze;++j) {
		/* run through the neighb tri of vert i 
		take care of the distribution of ASA between triangles */
			u=ind_tr[lze[j]] - 1 ;
			v=ind_tr[lze[j]+Ntri] - 1 ;
			w=ind_tr[lze[j]+2*Ntri] - 1 ;

			if (i==w) {
				r1[0] = XYZv[u] ;
					r1[1] = XYZv[u+Nvert] ;
					r1[2] = XYZv[u+2*Nvert] ; 
				r2[0] = XYZv[v] ;
					r2[1] = XYZv[v+Nvert] ;
					r2[2] = XYZv[v+2*Nvert] ; 
			}
			if (i==v) {
				r1[0] = XYZv[w] ;
					r1[1] = XYZv[w+Nvert] ;
					r1[2] = XYZv[w+2*Nvert] ; 
				r2[0] = XYZv[u] ;
					r2[1] = XYZv[u+Nvert] ;
					r2[2] = XYZv[u+2*Nvert] ; 
			}
			if (i==u) {
				r1[0] = XYZv[v] ;
					r1[1] = XYZv[v+Nvert] ;
					r1[2] = XYZv[v+2*Nvert] ; 
				r2[0] = XYZv[w] ;
					r2[1] = XYZv[w+Nvert] ;
					r2[2] = XYZv[w+2*Nvert] ; 
			}
			/* i is one of the 3 indices of the triangle  */
			/* Take care of the other of the vertices !!! */

			r3[0] = tri_pt4[lze[j]] ;
				r3[1] = tri_pt4[lze[j]+Ntri] ;
				r3[2] = tri_pt4[lze[j]+2*Ntri] ;
			/* triangle 4th point for sphere fitting */

			/*printf("r1 = %f %f %f ; r2 = %f %f %f ; 
			r3 = %f %f %f\n",r1[0],r1[1],r1[2], r2[0],r2[1],r2[2], 
			r3[0],r3[1],r3[2]) ;*/

			a[0][0] = si[0]+r1[0] ;
				a[0][1] = si[1]+r1[1] ; a[0][2] = si[2]+r1[2] ;
			a[1][0] = si[0]+r2[0] ;
				a[1][1] = si[1]+r2[1] ; a[1][2] = si[2]+r2[2] ;
			a[2][0] = si[0]+r3[0] ;
				a[2][1] = si[1]+r3[1] ; a[2][2] = si[2]+r3[2] ;

			b[0] = (si[0]*si[0]+si[1]*si[1]+si[2]*si[2]-
				(r1[0]*r1[0]+r1[1]*r1[1]+r1[2]*r1[2]))/2 ;
			b[1] = (si[0]*si[0]+si[1]*si[1]+si[2]*si[2]-
				(r2[0]*r2[0]+r2[1]*r2[1]+r2[2]*r2[2]))/2 ;
			b[2] = (si[0]*si[0]+si[1]*si[1]+si[2]*si[2]-
				(r3[0]*r3[0]+r3[1]*r3[1]+r3[2]*r3[2]))/2 ;

			deter = (a[0][0]*a[1][1]*a[2][2]+
				a[0][1]*a[1][2]*a[2][0]+a[0][2]*a[1][0]*a[2][1])
				- (a[2][0]*a[1][1]*a[0][2]+
				a[1][0]*a[0][1]*a[2][2]+a[0][0]*a[2][1]*a[1][2]) ;
			rc[0] = ((b[0]*a[1][1]*a[2][2]+
				a[0][1]*a[1][2]*b[2]+a[0][2]*b[1]*a[2][1])
				- (b[2]*a[1][1]*a[0][2]+b[1]*a[0][1]*a[2][2]+
				b[0]*a[2][1]*a[1][2]))/deter ;
			rc[1] = ((a[0][0]*b[1]*a[2][2]+
				b[0]*a[1][2]*a[2][0]+a[0][2]*a[1][0]*b[2])
				- (a[2][0]*b[1]*a[0][2]+a[1][0]*b[0]*a[2][2]+
				a[0][0]*b[2]*a[1][2]))/deter ;
			rc[2] = ((a[0][0]*a[1][1]*b[2]+
				a[0][1]*b[1]*a[2][0]+b[0]*a[1][0]*a[2][1])
				- (a[2][0]*a[1][1]*b[0]+a[1][0]*a[0][1]*b[2]+
				a[0][0]*a[2][1]*b[1]))/deter ;
			R = sqrt(pow(si[0]-rc[0],2) + pow(si[1]-rc[1],2) +
					pow(si[2]-rc[2],2)) ;
			/* rc and R = center and radius of the fitted sphere */

			/*printf("rc = %f %f %f ; R = %f \n", 
						rc[0],rc[1],rc[2],R) ;*/

			psi1 = 2*asin(sqrt(pow(si[0]-r1[0],2) +
				pow(si[1]-r1[1],2) + pow(si[2]-r1[2],2))/(2*R)) ;
			psi2 = 2*asin(sqrt(pow(si[0]-r2[0],2) +
				pow(si[1]-r2[1],2) + pow(si[2]-r2[2],2))/(2*R)) ;

			/*printf("psi1 = %f , psi2 = %f \n",psi1,psi2) ;*/

			r_b[0] = rc[0]+(sin(psi1-psi2)*(si[0]-rc[0])+
				sin(psi2)*(r1[0]-rc[0]))/sin(psi1) ;
			r_b[1] = rc[1]+(sin(psi1-psi2)*(si[1]-rc[1])+
				sin(psi2)*(r1[1]-rc[1]))/sin(psi1) ;
			r_b[2] = rc[2]+(sin(psi1-psi2)*(si[2]-rc[2])+
				sin(psi2)*(r1[2]-rc[2]))/sin(psi1) ;

			/*printf("rb = %f %f %f \n",r_b[0],r_b[1],r_b[2]) ;*/

			phi12 = 2*asin(sqrt(pow(r_b[0]-r2[0],2)+
					pow(r_b[1]-r2[1],2)+pow(r_b[2]-r2[2],2))
					/(2*R*sin(psi2))) ;

			/*printf("phi12 = %f \n",phi12) ;*/
			Ome_n[j] = (psi1+psi2)/4*phi12 ;
			Ome_ro = Ome_ro + Ome_n[j] ;
			/* sum of SA of the neighb triangles */

			/*printf("Ome_ro = %f \n",Ome_ro) ;*/
		}

		for (j=0;j<nze;++j) { 
		/* run through the neighb tri of vert i 
		distribution of ASA between the vert of tri */
			u=ind_tr[lze[j]] - 1 ;
			v=ind_tr[lze[j]+Ntri] - 1 ;
			w=ind_tr[lze[j]+2*Ntri] - 1 ;
			C[i+u*Nvert] = C[i+u*Nvert] 
						+ Ome_n[j]/Ome_ro * di2pi / 3 * weight ;
			C[i+v*Nvert] = C[i+v*Nvert] 
						+ Ome_n[j]/Ome_ro * di2pi / 3 * weight ;
			C[i+w*Nvert] = C[i+w*Nvert] 
						+ Ome_n[j]/Ome_ro * di2pi / 3 * weight ;

		}



		/* remove the identity matrix */
		C[i+i*Nvert] = C[i+i*Nvert] - 1 ;
	}
	/* Deflation */
	for (i=0;i<(Nvert*Nvert);++i) {
		C[i] = C[i] - defl ;
	}
}

/***************************************************************************/

void mexFunction(
	int nlhs, mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	double *Cii_c ;
	double *ind_tr, *XYZv, weight, defl, *tri_pt4 ;
	unsigned int mind, nind, mXYZv, nXYZv, mtri_pt4, ntri_pt4 ;

	if (nrhs!=5) {
		mexErrMsgTxt("5 inputs required: ind_tr, XYZv, weight, defl, tri_pt4") ;
	} else if (nlhs>1) {
		mexErrMsgTxt("Only ONE output : Cii_lin") ;
	}

	mind=mxGetM(prhs[0]) ;
	nind=mxGetN(prhs[0]) ;
	if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) 
		|| mxIsSparse(prhs[0]) || (nind!=3)) {
		mexErrMsgTxt("tri must be a matrix : Ntri x 3") ;
	}

	mXYZv=mxGetM(prhs[1]) ;
	nXYZv=mxGetN(prhs[1]) ;
	if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) 
		|| mxIsSparse(prhs[1]) || (nXYZv!=3) ) {
		mexErrMsgTxt("XYZv must be a matrix : Nvert x 3") ;
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


	mtri_pt4=mxGetM(prhs[4]) ;
	ntri_pt4=mxGetN(prhs[4]) ;
	if (!mxIsNumeric(prhs[4]) || mxIsComplex(prhs[4]) 
		|| mxIsSparse(prhs[4]) || (ntri_pt4!=3) || (mtri_pt4!=mind)) {
		mexErrMsgTxt("tri_pt4 must be a matrix : Ntri x 3") ;
	}

	plhs[0]=mxCreateDoubleMatrix(mXYZv,mXYZv,mxREAL) ; 

	Cii_c = mxGetPr(plhs[0]) ;
	ind_tr = mxGetPr(prhs[0]) ;
	XYZv = mxGetPr(prhs[1]) ;
	weight = mxGetScalar(prhs[2]) ;
	defl = mxGetScalar(prhs[3]) ;
	tri_pt4 = mxGetPr(prhs[4]) ;

	Cii_constant(Cii_c, ind_tr, mind, XYZv, mXYZv, weight, defl, tri_pt4) ;

	mxSetPr(plhs[0],Cii_c) ;
}
