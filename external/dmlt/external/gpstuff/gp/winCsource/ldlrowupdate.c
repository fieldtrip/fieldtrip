/* LDLROWUPDATE a function to conduct row-update for LDL Cholesky factorization
 *
 * Usage:
 *
 *	LD = ldlrowupdate (i1,LD,C,'+')	update an LDL' factorization after inclusion of row
 *	LD = ldlrowupdate (i1,LD,C,'-')	downdate an LDL' factorization after removal of row
 *
 * Here i1 is to row that is modified in the matrix K = LDL. LD is the
 * LDL Cholesky factorisation of K and C is the added/removed row of K.
 *
 * See Davis and Hager 2005 (Row Modification of a sparse Cholesky 
 * factorization) section 4 for details of the algorithm.
 */

 /* -----------------------------------------------------------------------------
 * The function requires SuiteSparse package by Timothy A. Davis to work. 
 *
 * Part of the code in this file is copied from CHOLMOD/MATLAB/ldlupdate.c 
 * in SuiteSparse version 3.2.0. 
 * -------------------------------------------------------------------------- */

/* -----------------------------------------------------------------------------
 * Copyright (C) 2005-2006 Timothy A. Davis
 * Copyright (c) 2009-2010      Jarno Vanhatalo
 *
 * This software is distributed under the GNU General Public
 * License (version 3 or later); please refer to the file
 * License.txt, included with the software, for details.
 * -----------------------------------------------------------------------------
 */


#include "cholmod_matlab.h"

void mexFunction
(
    int nargout,
    mxArray *pargout [ ],
    int nargin,
    const mxArray *pargin [ ]
)
{
    int ki;
    double dummy = 0 ;
    double *Lx, *Lx2 ;
    Int *Li, *Lp, *Li2, *Lp2, *Lnz2, *ColCount ;
    cholmod_sparse Cmatrix, *R, *Lsparse ;
    cholmod_factor *L ;
    cholmod_common Common, *cm ;
    Int j, k, s, update, n, lnz ;
    char buf [LEN] ;

    /* ---------------------------------------------------------------------- */
    /* start CHOLMOD and set parameters */ 
    /* ---------------------------------------------------------------------- */

    cm = &Common ;
    cholmod_l_start (cm) ;
    sputil_config (SPUMONI, cm) ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    if (nargout > 1 || nargin < 3 || nargin > 4)
    {
	mexErrMsgTxt ("Usage: L = ldlrowupdate (k, L, R, '+')") ; 
    }

    n = mxGetN (pargin [1]) ;
    k = mxGetN (pargin [2]) ;

    if (!mxIsSparse (pargin [1]) || !mxIsSparse (pargin [2])
	    || n != mxGetM (pargin [1]) || n != mxGetM (pargin [2])
	    || mxIsComplex (pargin [1]) || mxIsComplex (pargin [2]))
    {
      k = mxGetM (pargin [2]);
      j = mxGetM (pargin [1]);
      printf("n=%d  L=%d  R=%d \n", n, j, k);
	mexErrMsgTxt ("ldlrowupdate: R and/or L not sparse, complex, or wrong"
		" dimensions") ;
    }

    /* ---------------------------------------------------------------------- */
    /* determine if we're doing an update or downdate */
    /* ---------------------------------------------------------------------- */

    update = TRUE ;
    if (nargin > 3 && mxIsChar (pargin [3]))
    {
	mxGetString (pargin [3], buf, LEN) ;
	if (buf [0] == '-')
	{
	    update = FALSE ;
	}
	else if (buf [0] != '+')
	{
	    mexErrMsgTxt ("ldlrowupdate: update string must be '+' or '-'") ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* get ki: column integer of update */
    /* ---------------------------------------------------------------------- */
    ki = (int) *mxGetPr(pargin[0]);
    ki = ki-1;

    /* ---------------------------------------------------------------------- */
    /* get R: sparse matrix of incoming/outgoing columns */
    /* ---------------------------------------------------------------------- */

    R = sputil_get_sparse (pargin [2], &Cmatrix, &dummy, 0) ;

    /* ---------------------------------------------------------------------- */
    /* construct a copy of the input sparse matrix L */
    /* ---------------------------------------------------------------------- */

    /* get the MATLAB L */
    Lp = (Int *) mxGetJc (pargin [1]) ;
    Li = (Int *) mxGetIr (pargin [1]) ;
    Lx = mxGetPr (pargin [1]) ;

    /* allocate the CHOLMOD symbolic L */
    L = cholmod_l_allocate_factor (n, cm) ;
    L->ordering = CHOLMOD_NATURAL ;
    ColCount = L->ColCount ;
    for (j = 0 ; j < n ; j++)
    {
	ColCount [j] = Lp [j+1] - Lp [j] ;
    }

    /* allocate space for a CHOLMOD LDL' packed factor */
    cholmod_l_change_factor (CHOLMOD_REAL, FALSE, FALSE, TRUE, TRUE, L, cm) ;

    /* copy MATLAB L into CHOLMOD L */
    Lp2 = L->p ;
    Li2 = L->i ;
    Lx2 = L->x ;
    Lnz2 = L->nz ;
    lnz = L->nzmax ;
    for (j = 0 ; j <= n ; j++)
    {
	Lp2 [j] = Lp [j] ;
    }
    for (j = 0 ; j < n ; j++)
    {
	Lnz2 [j] = Lp [j+1] - Lp [j] ;
    }
    for (s = 0 ; s < lnz ; s++)
    {
	Li2 [s] = Li [s] ;
    }
    for (s = 0 ; s < lnz ; s++)
    {
	Lx2 [s] = Lx [s] ;
    }

    /* ---------------------------------------------------------------------- */
    /* update/downdate the LDL' factorization */
    /* ---------------------------------------------------------------------- */
    /* add row */
    if (update){
      if (!cholmod_l_rowadd (ki, R, L, cm))
	{
	  mexErrMsgTxt ("rowadd failed\n") ;
	}
    }
    /* delete row */
    else {
      if (!cholmod_l_rowdel (ki, NULL, L, cm))
	{
	  mexErrMsgTxt ("rowdel failed\n") ;
	}
    }
      

    /* ---------------------------------------------------------------------- */
    /* copy the results back to MATLAB */
    /* ---------------------------------------------------------------------- */

    /* change L back to packed LDL' (it may have become unpacked if the
     * sparsity pattern changed).  This change takes O(n) time if the pattern
     * of L wasn't updated. */
    Lsparse = cholmod_l_factor_to_sparse (L, cm) ;

    /* return L as a sparse matrix */
    pargout [0] = sputil_put_sparse (&Lsparse, cm) ;

    /* ---------------------------------------------------------------------- */
    /* free workspace and the CHOLMOD L, except for what is copied to MATLAB */
    /* ---------------------------------------------------------------------- */

    cholmod_l_free_factor (&L, cm) ;
    cholmod_l_finish (cm) ;
    cholmod_l_print_common (" ", cm) ;
    /*
    if (cm->malloc_count != 3 + mxIsComplex (pargout[0])) mexErrMsgTxt ("!") ;
    */
}
