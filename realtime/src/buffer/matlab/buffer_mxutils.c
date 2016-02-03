/*
 * Copyright (C) 2010, Stefan Klanke
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
*/

#include "mex.h"
#include "matrix.h"
#include "buffer.h"
#include "buffer_mxutils.h"

/* Like the normal fieldtrip append, but using mxMalloc and mxRealloc. This gives automatic
   cleanup after mexErrMsgTxt calls etc. If you free the memory yourself, use mxFree !
*/
unsigned int ft_mx_append(void **buf, unsigned int size, const void *addbuf, unsigned int addsize) {
	if (((  *buf != NULL) && (   size == 0)) || 
		((  *buf == NULL) && (   size != 0)) ||
		((addbuf == NULL) && (addsize != 0))) {
		mexErrMsgTxt("Internal error: bad call to ft_mx_append. Please report as a bug.");
	}
	if (addsize == 0) return size;
	
	/* these will yield an error message in Matlab if memory is not available */
	if (*buf == NULL) {
		*buf = mxMalloc(addsize);
	} else {
		*buf = mxRealloc(*buf, size + addsize);
	}
	
	memcpy((char *)(*buf) + size, addbuf, addsize);
	
	return size + addsize;
}

mxClassID class_id_from_ft_type(UINT32_T type) {
	switch(type) {
        case DATATYPE_UINT8:
			return mxUINT8_CLASS;
        case DATATYPE_UINT16:
			return mxUINT16_CLASS;
        case DATATYPE_UINT32:
			return mxUINT32_CLASS;
        case DATATYPE_UINT64:
			return mxUINT64_CLASS;
        case DATATYPE_INT8:
			return mxINT8_CLASS;
        case DATATYPE_INT16:
			return mxINT16_CLASS;
        case DATATYPE_INT32:
			return mxINT32_CLASS;
        case DATATYPE_INT64:
			return mxINT64_CLASS;
		case DATATYPE_FLOAT32:
			return mxSINGLE_CLASS;
		case DATATYPE_FLOAT64:
			return mxDOUBLE_CLASS;
		case DATATYPE_CHAR:
			return mxCHAR_CLASS;
	}
	return mxUNKNOWN_CLASS;
}

UINT32_T ft_type_from_class_id(mxClassID cid) {
	switch(cid) {
		case mxCHAR_CLASS:
			return DATATYPE_CHAR;
		case mxDOUBLE_CLASS:
			return DATATYPE_FLOAT64;
		case mxSINGLE_CLASS:
			return DATATYPE_FLOAT32;
		case mxINT8_CLASS:
			return DATATYPE_INT8;
		case mxINT16_CLASS:
			return DATATYPE_INT16;
		case mxINT32_CLASS:
			return DATATYPE_INT32;
		case mxINT64_CLASS:
			return DATATYPE_INT64;
		case mxUINT8_CLASS:
			return DATATYPE_UINT8;
		case mxUINT16_CLASS:
			return DATATYPE_UINT16;
		case mxUINT32_CLASS:
			return DATATYPE_UINT32;
		case mxUINT64_CLASS:
			return DATATYPE_UINT64;			
	}
	return DATATYPE_UNKNOWN;
}

mxArray *matrix_from_ft_type_data(UINT32_T type, UINT32_T rows, UINT32_T cols, const void *data) {
	mxArray *A;

	switch (type) {
		case DATATYPE_CHAR:
			{
				/* TODO: is there a proper opposite of mxArrayToString that we could use here??? 
					DATATYPE_CHAR in the FieldTrip buffer is defined to be 8 bits, but 
					mxChar elements are 16-bit unicode (I think).
				*/
				mxChar *dest;
				mwSize dim[2];
				const char *src = (const char *) data;
				UINT32_T n;
				
				dim[0] = rows;
				dim[1] = cols;
				
				A = mxCreateCharArray(2, dim);
				dest = mxGetChars(A);
				for (n=0;n<rows*cols;n++) {
					*dest++ = *src++;
				}
			}
			break;
        case DATATYPE_UINT8:
			A = mxCreateNumericMatrix(rows, cols, mxUINT8_CLASS, mxREAL);
			memcpy(mxGetData(A), data, rows*cols*WORDSIZE_UINT8);
			break;
        case DATATYPE_UINT16:
			A = mxCreateNumericMatrix(rows, cols, mxUINT16_CLASS, mxREAL);
			memcpy(mxGetData(A), data, rows*cols*WORDSIZE_UINT16);
			break;
        case DATATYPE_UINT32:
			A = mxCreateNumericMatrix(rows, cols, mxUINT32_CLASS, mxREAL);
			memcpy(mxGetData(A), data, rows*cols*WORDSIZE_UINT32);
			break;
        case DATATYPE_UINT64:
			A = mxCreateNumericMatrix(rows, cols, mxUINT64_CLASS, mxREAL);
			memcpy(mxGetData(A), data, rows*cols*WORDSIZE_UINT64);
			break;
			
        case DATATYPE_INT8:
			A = mxCreateNumericMatrix(rows, cols, mxINT8_CLASS, mxREAL);
			memcpy(mxGetData(A), data, rows*cols*WORDSIZE_INT8);
			break;
        case DATATYPE_INT16:
			A = mxCreateNumericMatrix(rows, cols, mxINT16_CLASS, mxREAL);
			memcpy(mxGetData(A), data, rows*cols*WORDSIZE_INT16);
			break;
        case DATATYPE_INT32:
			A = mxCreateNumericMatrix(rows, cols, mxINT32_CLASS, mxREAL);
			memcpy(mxGetData(A), data, rows*cols*WORDSIZE_INT32);
			break;
        case DATATYPE_INT64:
			A = mxCreateNumericMatrix(rows, cols, mxINT64_CLASS, mxREAL);
			memcpy(mxGetData(A), data, rows*cols*WORDSIZE_INT64);
			break;
			
        case DATATYPE_FLOAT32:
			A = mxCreateNumericMatrix(rows, cols, mxSINGLE_CLASS, mxREAL);
			memcpy(mxGetData(A), data, rows*cols*WORDSIZE_FLOAT32);
			break;
        case DATATYPE_FLOAT64:
			A = mxCreateNumericMatrix(rows, cols, mxDOUBLE_CLASS, mxREAL);
			memcpy(mxGetData(A), data, rows*cols*WORDSIZE_FLOAT64);
			break;
			
        default:
			A = mxCreateDoubleMatrix(0,0, mxREAL);
			/* TODO: should we rather return NULL here, or maybe
				A = mxCreateString("--invalid data type--");
			*/
	}
	return A;
}


