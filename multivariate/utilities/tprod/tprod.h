#ifndef tprodH
#define tprodH
#include "mxInfo.h"

/*
  
  Tensor multiplication routine headers.

  $Id: tprod.h,v 1.1.1.1 2008/02/27 14:42:55 roboos Exp $

 Copyright 2006-     by Jason Farquhar (jdrf@zepler.org)
 Permission is granted for anyone to copy, use, or modify this
 software and accompanying documents for any uncommercial purposes,
 provided this copyright notice is retained, and note is made of
 any changes that have been made. This software and documents are
 distributed without any warranty, express or implied.


 */

/* use MATLAB memory management routines if in matlab mode */
#ifdef MATLAB_MEX_FILE
#include "mex.h"
#include "matrix.h"
#define MALLOC mxMalloc
#define CALLOC mxCalloc
#define FREE   mxFree
#define BOOLTYPE mxLogical
#else
#include <stdlib.h>
#define MALLOC malloc
#define CALLOC calloc
#define FREE   free
#define BOOLTYPE char
#endif

typedef enum { RXRY=0, CXRY=1, RXCY=2, CXCY=3} RealComplexInputTypes; /* enum for the bits of input cases */

#ifndef MAX
#define MAX(A,B)  ((A) > (B) ? (A) : (B))
#define MIN(A,B)  ((A) < (B) ? (A) : (B))
#endif

/* check the compilier state to use the appropriate inline directive */
#ifdef __GNUC__ /* use the GNUC special form */
#define INLINE __inline__
#elif defined(__STDC__) && __STDC_VERSION__ >= 199901L /*C99 compat compilier*/
#define INLINE static inline
#else /* fall back on C89 version, i.e. *no inlines* */
#define INLINE
#endif

typedef enum { /* types of error we could encounter in repop */
  OK=0, 
  ZTYPEMISMATCH, 
  XTYPEMISMATCH, 
  YTYPEMISMATCH, 
  INTYPEMISMATCH, 
  OTHERERROR  
} TprodErrorCode;


/* compute the output info from xinf, yinf and list of accumulated dims */
MxInfo initzmxInfo(int znd, const MxInfo xinfo, const MxInfo yinfo,
						 const int *x2yIdx, int xnidx, int ynidx);

/* compute the accumulated and left over sub-matrices info */
void initrestmaccmxInfo(int znd,
							  const MxInfo xinfo, const MxInfo yinfo, 
							  const int x2yIdx[], int xnidx, int ynidx,
							  MxInfo *xrestinfo,  MxInfo *yrestinfo,
							  MxInfo *xmaccinfo,  MxInfo *ymaccinfo);

void squeezemxInfoPair(MxInfo *xinf, MxInfo *yinf);

void optimisetprodQuery(MxInfo *zrest, MxInfo *xrest, MxInfo *yrest, 
								MxInfo *xmacc, MxInfo *ymacc);

TprodErrorCode tprod(const MxInfo zinfo, const MxInfo xrest, const MxInfo yrest,
							const MxInfo xmacc, const MxInfo ymacc,
							int blksz);
#endif
