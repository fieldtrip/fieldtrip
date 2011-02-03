#ifndef tprodH
#define tprodH
#include "mxInfo.h"

/*
  
  Tensor multiplication routine headers.

  $Id$

 Copyright 2006-     by Jason Farquhar (jdrf@zepler.org)
 Permission is granted for anyone to copy, use, or modify this
 software and accompanying documents for any uncommercial purposes,
 provided this copyright notice is retained, and note is made of
 any changes that have been made. This software and documents are
 distributed without any warranty, express or implied.


 */

/* use MATLAB memory management routines if in matlab mode */
#ifdef MATLAB_MEX_FILE
#define BOOLTYPE mxLogical
#else
#define BOOLTYPE char
#endif

#include <stdlib.h>

/* Library functions we use to provide memory management and error handling */
/* N.B. *YOU* must provide implementations of these prototypes */
void *CALLOC(size_t nmemb, size_t size);
void *MALLOC(size_t size);
void FREE(void *ptr);
void ERROR(const char *msg);
void WARNING(const char *msg);

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
  UNSUPPORTEDINPUTS,
  LESSXIND,
  LESSYIND,
  DUPXLABEL,
  DUPYLABEL,
  NONSINGLENTONIGNORED,
  NOMATCHINGXY,
  NOMATCHINGXYACC,
  MATCHDIMMISMATCH,
  XYTOOSMALL,
  ACCDIMMISMATCH,
  NONCONTIGUOUSZ,
  OTHERERROR
} TprodErrorCode;

/* compute the x2yIdx version of the input dimspecs */
TprodErrorCode compx2yIdx(MxInfo xinfo, int xnidx, int *xidx,
					MxInfo yinfo, int ynidx, int *yidx,
									int**x2yidxp, int *znd, int *maccnd, int *seqnd);
TprodErrorCode compx2yIdx_dd(MxInfo xinfo, int xnidx, double *xidx,
						 MxInfo yinfo, int ynidx, double *yidx,
									  int **x2yidxp, int *znd, int *maccnd, int *seqnd);

/* compute the output info from xinf, yinf and list of accumulated dims */
MxInfo initzmxInfo(int znd, const MxInfo xinfo, const MxInfo yinfo,
						 const int *x2yIdx, int xnidx, int ynidx);

/* compute the accumulated and left over sub-matrices info */
TprodErrorCode  initrestmaccmxInfo(int znd,
											  const MxInfo xinfo, const MxInfo yinfo, 
											  const int x2yIdx[], int xnidx, int ynidx,
											  MxInfo *xrestinfo,  MxInfo *yrestinfo,
											  MxInfo *xmaccinfo,  MxInfo *ymaccinfo);

void squeezemxInfoPair(MxInfo *xinf, MxInfo *yinf);

TprodErrorCode optimisetprodQuery(MxInfo *zrest, MxInfo *xrest, MxInfo *yrest,
											  MxInfo *xmacc, MxInfo *ymacc);

TprodErrorCode tprod(const MxInfo zinfo, const MxInfo xrest, const MxInfo yrest,
							const MxInfo xmacc, const MxInfo ymacc,
							int blksz);
#endif
