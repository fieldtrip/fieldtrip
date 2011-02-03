/*
 * Copyright (C) 2010, Stefan Klanke
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
*/
#ifndef __buffer_mxutils_h
#define __buffer_mxutils_h

#include "mex.h"
#include "matrix.h"
#include "buffer.h"


/* This is a fix to allow compilation on very old Matlab versions (e.g., 7.1), 
   where mwSize is not available yet and int's are used instead.
*/
#ifndef mwSize
#define mwSize int
#endif


/* Like the normal fieldtrip append, but using mxMalloc and mxRealloc. This gives automatic
   error reporting and automated cleanup after a call to mxErrMsgTxt elsewhere. 
   If you free the memory yourself pointed to by *buf yourself, use mxFree !
*/
unsigned int ft_mx_append(void **buf, unsigned int size, const void *addbuf, unsigned int addsize);

/* Convert FieldTrip datatype to Matlab datatype (class ID) */
mxClassID class_id_from_ft_type(UINT32_T type);

/* Convert Matlab datatype (class ID) to FieldTrip datatype */
UINT32_T ft_type_from_class_id(mxClassID cid);

/* Convenience wrapper: Get FieldTrip data type from that of a Matlab array */
#define ft_type_from_array(A)  ft_type_from_class_id(mxGetClassID(A))

/* Creates a 2D Matlab array from a given Fieldtrip type, rows x columns, and data pointer */
mxArray *matrix_from_ft_type_data(UINT32_T type, UINT32_T rows, UINT32_T cols, const void *data);

#endif
