/* SPINV    Evaluate the sparse inverse matrix
 *
 * z = spinv(A)  returns the elements of inv(A)_ij, for which A_ij
 *      is different from zero. 
 *
 * z = spinv(LD, 1)  returns the elements of inv(A)_ij, for which A_ij
 *      is different from zero, and where LD is the LDL cholesky
 *      decomposition of A. LD has to be in the form returned by ldlchol.
 *
 *
 *   Note! If z = spinv(LD, 1) is used LD must not be modified in Matlab
 *   after ldlchol. Matlab destroys the symbolic sparsity structure in
 *   the Cholesky decomposition, which is needed in the spinv
 *   algorithm. If LD is modified the worst scenario is memory corruption. 
 *   
 *
 *     For details, see:
 *     Jarno Vanhatalo and Aki Vehtari (2008). Modelling local and
 *     global phenomena with sparse Gaussian processes. Proceedings of
 *     the 24th Conference on Uncertainty in Artificial Intelligence
 *
 *
 */

 /* -----------------------------------------------------------------------------
 * The function requires SuiteSparse package by Timothy A. Davis to work. 
 *
 * Part of the code in this file is copied from CHOLMOD/MATLAB/ldlchol.c 
 * in SuiteSparse version 3.2.0. 
 * -------------------------------------------------------------------------- */

/* -----------------------------------------------------------------------------
 * Copyright (C) 2005-2006 Timothy A. Davis
 * Copyright (c) 2008-2010 Jarno Vanhatalo
 *
 * This software is distributed under the GNU General Public
 * License (version 3 or later); please refer to the file
 * License.txt, included with the software, for details.
 * -----------------------------------------------------------------------------
 */

#include <stdio.h>
#include "cholmod_matlab.h"
#include "mex.h"
#include "matrix.h"
#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))

void cumsum2 (mwIndex *p, mwIndex *c, mwSize n);

void mexFunction
(
    int	nargout,
    mxArray *pargout [ ],
    int	nargin,
    const mxArray *pargin [ ]
    )
{
  double dummy = 0, beta [2], *px, *C, *Ct, *C2, *fil, *Zt, *zt, done=1.0, *zz, dzero=0.0;
  cholmod_sparse Amatrix, *A, *Lsparse ;
  cholmod_factor *L ;
  cholmod_common Common, *cm ;
  Int minor, *It2, *Jt2 ;
  mwIndex l, k2, h, k, i, j, ik, *I, *J, *Jt, *It, *I2, *J2, lfi, *w, *w2, *r;
  mwSize nnz, nnzlow, m, n;
  int nz = 0;
  mwSignedIndex one=1, lfi_si;
  mxArray *Am, *Bm;
  char *uplo="L", *trans="N";
  

  /* ---------------------------------------------------------------------- */
  /* Only one input. We have to find first the Cholesky factorization.      */ 
  /* start CHOLMOD and set parameters */ 
  /* ---------------------------------------------------------------------- */

  if (nargin == 1) {
    cm = &Common ;
    cholmod_l_start (cm) ;
    sputil_config (SPUMONI, cm) ;
    
    /* convert to packed LDL' when done */
    cm->final_asis = FALSE ;
    cm->final_super = FALSE ;
    cm->final_ll = FALSE ;
    cm->final_pack = TRUE ;
    cm->final_monotonic = TRUE ;

    /* since numerically zero entries are NOT dropped from the symbolic
     * pattern, we DO need to drop entries that result from supernodal
     * amalgamation. */
    cm->final_resymbol = TRUE ;

    cm->quick_return_if_not_posdef = (nargout < 2) ;
  }

  /* This will disable the supernodal LL', which will be slow. */
  /* cm->supernodal = CHOLMOD_SIMPLICIAL ; */
  
  /* ---------------------------------------------------------------------- */
  /* get inputs */
  /* ---------------------------------------------------------------------- */
  
  if (nargin > 3)
    {
      mexErrMsgTxt ("usage: Z = sinv(A), or Z = sinv(LD, 1)") ;
    }
  
  n = mxGetM (pargin [0]) ;
  m = mxGetM (pargin [0]) ;
  
  if (!mxIsSparse (pargin [0]))
    {
      mexErrMsgTxt ("A must be sparse") ;
    }
  if (n != mxGetN (pargin [0]))
    {
      mexErrMsgTxt ("A must be square") ;
    }

  /* Only one input. We have to find first the Cholesky factorization.      */
  if (nargin == 1) {
    /* get sparse matrix A, use tril(A)  */
    A = sputil_get_sparse (pargin [0], &Amatrix, &dummy, -1) ; 
    
    A->stype = -1 ;	    /* use lower part of A */
    beta [0] = 0 ;
    beta [1] = 0 ;
    
    /* ---------------------------------------------------------------------- */
    /* analyze and factorize */
    /* ---------------------------------------------------------------------- */
    
    L = cholmod_l_analyze (A, cm) ;
    cholmod_l_factorize_p (A, beta, NULL, 0, L, cm) ;
    
    if (cm->status != CHOLMOD_OK)
      {
	mexErrMsgTxt ("matrix is not positive definite") ;
      }
    
    /* ---------------------------------------------------------------------- */
    /* convert L to a sparse matrix */
    /* ---------------------------------------------------------------------- */

    Lsparse = cholmod_l_factor_to_sparse (L, cm) ;
    if (Lsparse->xtype == CHOLMOD_COMPLEX)
      {
	mexErrMsgTxt ("matrix is complex") ;
      }
    
    /* ---------------------------------------------------------------------- */
    /* Set the sparse Cholesky factorization in Matlab format */
    /* ---------------------------------------------------------------------- */
    /*Am = sputil_put_sparse (&Lsparse, cm) ;
      I = mxGetIr(Am);
      J = mxGetJc(Am);
      C = mxGetPr(Am);
      nnz = mxGetNzmax(Am); */

    It2 = Lsparse->i;
    Jt2 = Lsparse->p;
    Ct = Lsparse->x;
    nnz = (mwSize) Lsparse->nzmax;

    Am = mxCreateSparse(m, m, nnz, mxREAL) ;
    I = mxGetIr(Am);
    J = mxGetJc(Am);
    C = mxGetPr(Am);
    for (j = 0 ;  j < n+1 ; j++)  J[j] = (mwIndex) Jt2[j];
    for ( i = 0 ; i < nnz ; i++) {
	I[i] = (mwIndex) It2[i];
	C[i] = Ct[i];
    }
    
    cholmod_l_free_sparse (&Lsparse, cm) ;

    /*FILE *out1 = fopen( "output1.txt", "w" );
    if( out1 != NULL )
      fprintf( out1, "Hello %d\n", nnz );
      fclose (out1);*/
    
  } else {
    /* The cholesky factorization is given as an input.      */
    /* We have to copy it into workspace                     */
    It = mxGetIr(pargin [0]);
    Jt = mxGetJc(pargin [0]);
    Ct = mxGetPr(pargin [0]);
    nnz = mxGetNzmax(pargin [0]);
    
    Am = mxCreateSparse(m, m, nnz, mxREAL) ;
    I = mxGetIr(Am);
    J = mxGetJc(Am);
    C = mxGetPr(Am);
    for (j = 0 ;  j < n+1 ; j++)  J[j] = Jt[j];
    for ( i = 0 ; i < nnz ; i++) {
	I[i] = It[i];
	C[i] = Ct[i];
    }    
  }

  /* Evaluate the sparse inverse */
  C[nnz-1] = 1.0/C[J[m-1]];               /* set the last element of sparse inverse */
  fil = mxCalloc((mwSize)1,sizeof(double));
  zt = mxCalloc((mwSize)1,sizeof(double));
  Zt = mxCalloc((mwSize)1,sizeof(double));
  zz = mxCalloc((mwSize)1,sizeof(double));
  for (j=m-2;j!=-1;j--){
    lfi = J[j+1]-(J[j]+1);
    
    /* if (lfi > 0) */
    if ( J[j+1] > (J[j]+1) )
      {
	/*	printf("lfi = %u \n ", lfi);
	printf("lfi*double = %u \n", (mwSize)lfi*sizeof(double));
	printf("lfi*lfi*double = %u \n", (mwSize)lfi*(mwSize)lfi*sizeof(double));
	printf("\n \n");
	*/
	
	fil = mxRealloc(fil,(mwSize)lfi*sizeof(double));
	for (i=0;i<lfi;i++) fil[i] = C[J[j]+i+1];                   /* take the j'th lower triangular column of the Cholesky */
	
	zt = mxRealloc(zt,(mwSize)lfi*sizeof(double));              /* memory for the sparse inverse elements to be evaluated */
	Zt = mxRealloc(Zt,(mwSize)lfi*(mwSize)lfi*sizeof(double));  /* memory for the needed sparse inverse elements */
	
	/* Set the lower triangular for Zt */
	k2 = 0;
	for (k=J[j]+1;k<J[j+1];k++){
	  ik = I[k];
	  h = k2;
	  for (l=J[ik];l<=J[ik+1];l++){
	    if (I[l] == I[ J[j]+h+1 ]){
	      Zt[h+lfi*k2] = C[l];
	      h++;
	    }
	  }
	  k2++;
	}
	
	
	/* evaluate zt = fil*Zt */
	lfi_si = (mwSignedIndex) lfi;
	dsymv_(uplo, &lfi_si, &done, Zt, &lfi_si, fil, &one, &dzero, zt, &one);
	
	/* Set the evaluated sparse inverse elements, zt, into C */
	k=lfi-1;
	for (i = J[j+1]-1; i!=J[j] ; i--){
	  C[i] = -zt[k];
	  k--;
	}
	/* evaluate the j'th diagonal of sparse inverse */
	dgemv_(trans, &one, &lfi_si, &done, fil, &one, zt, &one, &dzero, zz, &one); 
	C[J[j]] = 1.0/C[J[j]] + zz[0];
      }
    else
      {
	/* evaluate the j'th diagonal of sparse inverse */
	C[J[j]] = 1.0/C[J[j]];	
      }
  }
    
  /* Free the temporary variables */
  mxFree(fil);
  mxFree(zt);
  mxFree(Zt);
  mxFree(zz);

  /* ---------------------------------------------------------------------- */
  /* Permute the elements according to r(q) = 1:n                           */
  /* Done only if the Cholesky was evaluated here                           */
  /* ---------------------------------------------------------------------- */
  if (nargin == 1) {
   
    Bm = mxCreateSparse(m, m, nnz, mxREAL) ;     
    It = mxGetIr(Bm);
    Jt = mxGetJc(Bm);
    Ct = mxGetPr(Bm);                            /* Ct = C(r,r) */ 
    
    r = (mwIndex *) L->Perm;                         /* fill reducing ordering */
    w = mxCalloc(m,sizeof(mwIndex));                 /* column counts of Am */
    
    /* count entries in each column of Bm */
    for (j=0; j<m; j++){
      k = r ? r[j] : j ;       /* column j of Bm is column k of Am */
      for (l=J[j] ; l<J[j+1] ; l++){
	i = I[l];
	ik = r ? r[i] : i ;    /* row i of Bm is row ik of Am */
	w[ max(ik,k) ]++;
      }
    }
    cumsum2(Jt, w, m);
    for (j=0; j<m; j++){
      k = r ? r[j] : j ;             /* column j of Bm is column k of Am */
      for (l=J[j] ; l<J[j+1] ; l++){
	i= I[l];
	ik = r ? r[i] : i ;          /* row i of Bm is row ik of Am */
	It [k2 = w[max(ik,k)]++ ] = min(ik,k);
	Ct[k2] = C[l];
      }
    }
    mxFree(w);
    
    /* ---------------------------------------------------------------------- */
    /* Transpose the permuted (upper triangular) matrix Bm into Am */
    /* (this way we get sorted columns)                            */
    /* ---------------------------------------------------------------------- */
    w = mxCalloc(m,sizeof(mwIndex));                 
    for (i=0 ; i<Jt[m] ; i++) w[It[i]]++;        /* row counts of Bm */
    cumsum2(J, w, m);                            /* row pointers */
    for (j=0 ; j<m ; j++){
      for (i=Jt[j] ; i<Jt[j+1] ; i++){
	I[ l=w[ It[i] ]++ ] = j;
	C[l] = Ct[i];
      }
    }
    mxFree(w);
    mxDestroyArray(Bm);
  }
  
  /* ---------------------------------------------------------------------- */
  /* Fill the upper triangle of the sparse inverse */
  /* ---------------------------------------------------------------------- */
  
  w = mxCalloc(m,sizeof(mwIndex));        /* workspace */
  w2 = mxCalloc(m,sizeof(mwIndex));       /* workspace */
  for (k=0;k<J[m];k++) w[I[k]]++;     /* row counts of the lower triangular */
  for (k=0;k<m;k++) w2[k] = w[k] + J[k+1] - J[k] - 1;   /* column counts of the sparse inverse */
  
  nnz = (mwSize)2*nnz - m;                       /* The number of nonzeros in Z */
  pargout[0] = mxCreateSparse(m,m,nnz,mxREAL);   /* The sparse matrix */
  It = mxGetIr(pargout[0]);
  Jt = mxGetJc(pargout[0]);
  Ct = mxGetPr(pargout[0]);
  
  cumsum2(Jt, w2, m);               /* column starting points */
  for (j = 0 ; j < m ; j++){           /* fill the upper triangular */
    for (k = J[j] ; k < J[j+1] ; k++){
      It[l = w2[ I[k]]++] = j ;	 /* place C(i,j) as entry Ct(j,i) */
      if (Ct) Ct[l] = C[k] ;
    }
  }
  for (j = 0 ; j < m ; j++){           /* fill the lower triangular */
    for (k = J[j]+1 ; k < J[j+1] ; k++){
      It[l = w2[j]++] = I[k] ;         /* place C(j,i) as entry Ct(j,i) */
      if (Ct) Ct[l] = C[k] ;
    }
  }
  
  mxFree(w2);
  mxFree(w);
  
  /* ---------------------------------------------------------------------- */
  /* return to MATLAB */
  /* ---------------------------------------------------------------------- */
  
  /* ---------------------------------------------------------------------- */
  /* free workspace and the CHOLMOD L, except for what is copied to MATLAB */
  /* ---------------------------------------------------------------------- */
  if (nargin == 1) {
    cholmod_l_free_factor (&L, cm) ;
    cholmod_l_finish (cm) ;
    cholmod_l_print_common (" ", cm) ;
  }
  mxDestroyArray(Am);
  
}

void cumsum2 (mwIndex *p, mwIndex *c, mwSize n)
{
  mwIndex i;
  mwIndex nz = 0;
  if(!p || !c) return;
  for (i=0;i<n;i++){
    p[i]=nz;
    nz+=c[i];
    c[i]=p[i];
  }
  p[n]=nz;
}
