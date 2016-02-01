#ifndef REPOPH_
#define REPOPH_
/*

  Generic header file for the replicating operators code

  $Id$

 Copyright 2006-     by Jason D.R. Farquhar (jdrf@zepler.org)
 Permission is granted for anyone to copy, use, or modify this
 software and accompanying documents for any uncommercial
 purposes, provided this copyright notice is retained, and note is
 made of any changes that have been made. This software and
 documents are distributed without any warranty, express or
 implied


 */



#ifdef MATLAB_MEX_FILE
#include "mex.h"
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

/* check the compilier state to use the appropriate inline directive */
#ifdef __GNUC__ /* use the GNUC special form */
#define INLINE __inline__
#elif defined(__STDC__) && __STDC_VERSION__ >= 199901L /*C99 compat compilier*/
#define INLINE static inline
#else /* fall back on C89 version, i.e. *no inlines* */
#define INLINE
#endif

/* N.B. use defines rather than enum so pre-processor can also test them */
#define PLUS    1
#define MINUS   2
#define TIMES   3
#define RDIVIDE 4 
#define LDIVIDE 5 
#define POWER   6 
#define EQ      7 
#define NE      8
#define LT      9
#define GT      10
#define LE      11 
#define GE      12
#define MINOP   13
#define MAXOP   14
#define MOD     15

/* enum { PLUS=1, MINUS, TIMES, RDIVIDE, LDIVIDE, POWER, EQ, NE, LT, GT, LE, GE }; */
typedef enum { RXRY=0, CXRY=1, RXCY=2, CXCY=3} RealComplexInputTypes; /* enum for the types input cases */

#ifndef MAX
#define min(A,B)  ((A) < (B) ? (A) : (B))
#define max(A,B)  ((A) > (B) ? (A) : (B))
#endif

/* macro to automatically construct the function name by concatenating:
   thename itself, the type of x/y inputs.
   N.B. we seem to need a double nested call to fix a (bug?) in gcc*/
#ifndef CAT
#define CAT2(A,B) A##B
#define CAT(A,B) CAT2(A,B)
#endif

char anycomplexCoeff(MxInfo *inf);
char removeZeroImag(MxInfo *zinf);
int getOpid(const char*);
MxInfo initzinfo(const MxInfo xinfo, const MxInfo yinfo, int repNonUnitDim);
void repopqueryOptimise(MxInfo *zinfo, MxInfo *xinfo, MxInfo *yinfo, 
								int *repxy);

typedef enum { /* types of error we could encounter in repop */
  OK=0, 
  ZTYPEMISMATCH, 
  XTYPEMISMATCH, 
  YTYPEMISMATCH, 
  INTYPEMISMATCH, 
  UNDEFOPERATOR,
  OTHERERROR  
} RepopErrorCode;

#if !defined(OPNM) && ! defined(OPID) 
/* only if we aren't in the 1-at-a-time mode */
RepopErrorCode repop(MxInfo zinfo, const MxInfo xinfo, const MxInfo yinfo, int opid);
#endif

#endif
