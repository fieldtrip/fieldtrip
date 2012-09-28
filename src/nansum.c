#include <matrix.h>
#include <mex.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "compiler.h"

/* Microsoft doesn't define isnan, but has an _isnan instead. We simply do: */
#define isnan(x) !(x == x)

/* Shortcuts: */
#define size(x) mxGetDimensions(x)
#define ndim(x) mxGetNumberOfDimensions(x)

/* Helper function to find linearized distance between elements on a given
 * dimension of a multidimensional array. Lower indices are close in memory
 * (column-major or FORTRAN/MATLAB order: */
mwSize stride(int dim, int shape_len, const mwSize *shape)
{
  int i, s = 1;
  for (i = 0; i < dim; i++)
    s *= shape[i];
  return s;
}

/* Convert a offset in a linearized array to a n-dimensional array index using
 * stride(). */
void offset_to_index(mwSize offset, int ndim, const mwSize *shape, 
  mwSize *index)
{
  int i, s;
  for (i = ndim; i-- > 0;) {
    s = stride(i, ndim, shape);
    index[i] = offset / s;
    offset -= index[i] * s;
  }
}

/* Convert an n-dimensional index for an multi-dimensional array to a
 * linearized offset using stride(). */
mwSize index_to_offset(int ndim, mwSize *index, const mwSize *shape)
{
  int i = 0;
  mwSize offset = 0;
  for (i = 0; i < ndim; ++i)
    offset += stride(i, ndim, shape) * index[i];
  return offset;
}

/* What is C a horrible language. To hoops we have to get through to make this
 * work for different data types...
 * Since overloading does not work, we generate type-specific functions: */
#define fname(name, suffix) name ## _ ## suffix
#define nansum_template(TYPE, INTERMEDIATE_TYPE)\
double fname(nansum, TYPE)(int n, TYPE *x0, mwSize stride) \
{\
  int i; INTERMEDIATE_TYPE result = 0;\
  for (i = 0; i < n; ++i) {\
    if (!isnan(x0[i * stride]))\
      result += x0[i * stride];\
  }\
  return result;\
}

nansum_template(float, float); /* Note that the calculations are performed with
                                  limited precision as well to be fully
                                  compatible with MATLABs nansum. */

nansum_template(double, double);
nansum_template(int8_T, double); nansum_template(uint8_T, double);
nansum_template(int16_T, double); nansum_template(uint16_T, double);
nansum_template(int32_T, double); nansum_template(uint32_T, double);
nansum_template(int64_T, double); nansum_template(uint64_T, double);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  const mxArray *X = prhs[0], *Y = plhs[0];
  mwSize *index, *size_X, *size_Y, stride_x;
  mwIndex j;
  int i, n, squash_dim, squash_len;
  mxClassID classid;

  /* Check and handle input. */
  switch(nrhs) {
    case 0: 
      mexErrMsgTxt("Too few input arguments!");
      break;

    case 1: /* No dimension is given. Find first non-singleton dimension: */
      squash_dim = ndim(X);
      for(i = ndim(X); i-- > 0;) {
        if(size(X)[i] > 1)
          squash_dim = i;
      }
      break;

    case 2: 
      squash_dim = (int) mxGetScalar(prhs[1]);
      if(squash_dim <= 0) {
        /* Note that squash_dim > ndim(X) is allowed! */
        mexErrMsgTxt("Invalid value for dimension (argument 2)!");
      }
      squash_dim -= 1; /* convert form MATLAB to C indexing */
      break;

    default:
      mexErrMsgTxt("Too many input arguments!");
  }
  printf("Found squash_dim = %d.\n", squash_dim);

  /* Create output array Y: */
  size_Y = mxMalloc(ndim(X) * sizeof(size_Y));
  memcpy(size_Y, size(X), ndim(X) * sizeof(size_Y));
  size_Y[squash_dim] = 1; 

  classid = mxGetClassID(X);  /* copy class ID from X */
  
  if (mxGetClassID(X) == mxSINGLE_CLASS) {
    /* Return single-precision if input is single precision. */
    Y = plhs[0] = mxCreateNumericArray(
      ndim(X), size_Y, mxSINGLE_CLASS, mxGetImagData(X) != NULL);
  }
  else {
    Y = plhs[0] = mxCreateNumericArray(
      ndim(X), size_Y, mxDOUBLE_CLASS, mxGetImagData(X) != NULL);
  }

  /*
  for(i = ndim(Y); i-- > 0;) {
    printf("i=%d, size(Y)[i] = %d\n", i, size(Y)[i]);
  }
  */

  /* Store calls to stat function with offset and stride in Y: */
  index = (mwSize *) mxMalloc(ndim(X) * sizeof(mwSize));
  stride_x = stride(squash_dim, ndim(X), size(X));
      
  /* Find the number of elements in our squashed dimension: */
  if(squash_dim >= ndim(X)) {
    /* MATLAB's nansum supports out-of-range dims to operate on. */
    squash_len = 1; 
  }
  else {
    squash_len = size(X)[squash_dim];
  }

  {
    void *src = mxGetData(X), *src_imag = mxGetImagData(X);
    double *dest = mxGetData(Y), *dest_imag = mxGetImagData(Y);

    printf("src_i = %p, dest_i = %p\n", src_imag, dest_imag);

    for (i = 0; i < mxGetNumberOfElements(Y); ++i) {
      /* For each element in the output array, do: */

      /* Transform the element number in an index in Y. */
      offset_to_index(i, ndim(Y), size(Y), index);
      /* And map this index back to an offset in X: */
      j = index_to_offset(ndim(X), index, size(X));
      
      switch (classid) {
        /* Oh, the fucking horror of dumb, statically typed languages... */
        case mxINT8_CLASS:
          dest[i] = nansum_int8_T(squash_len, (int8_T *) src + j, stride_x);
          break;

        case mxUINT8_CLASS:
          dest[i] = nansum_uint8_T(squash_len, (uint8_T *) src + j, stride_x);
          break;

        case mxINT16_CLASS:
          dest[i] = nansum_int16_T(squash_len, (int16_T *) src + j, stride_x);
          break;

        case mxUINT16_CLASS:
          dest[i] = nansum_uint16_T(squash_len, (uint16_T *) src + j, stride_x);
          break;

        case mxINT32_CLASS:
          dest[i] = nansum_int32_T(squash_len, (int32_T *) src + j, stride_x);
          break;

        case mxUINT32_CLASS:
          dest[i] = nansum_uint32_T(squash_len, (uint32_T *) src + j, stride_x);
          break;

        case mxINT64_CLASS:
          dest[i] = nansum_int64_T(squash_len, (int64_T *) src + j, stride_x);
          break;

        case mxUINT64_CLASS:
          dest[i] = nansum_uint64_T(squash_len, (uint64_T *) src + j, stride_x);
          break;

        case mxSINGLE_CLASS:
          {
            float *d = (float *) dest;
            float *d_i = (float *) dest_imag;
            float v = 0;

            d[i] = nansum_float(squash_len, (float *) src + j, stride_x);
            if (src_imag) {
              d_i[i] = (float) nansum_float(
                squash_len, (float *) src_imag + j, stride_x);
            }
          }
          break;

        case mxDOUBLE_CLASS:
          dest[i] = nansum_double(squash_len, (double *) src + j, stride_x);
          if (src_imag) {
            dest_imag[i] = nansum_double(squash_len, (double *) src_imag + j,                 stride_x);
          }
          break;

        case mxLOGICAL_CLASS:
          mexWarnMsgTxt(
            "The byte layout of logical arrays is not documented. Guessing.");
          dest[i] = nansum_uint8_T(squash_len, (uint8_T *) src + j, stride_x);
          break;

        case mxCHAR_CLASS:
          mexWarnMsgTxt(
            "The byte layout of char arrays is not documented. Guessing.");
          dest[i] = nansum_uint16_T(squash_len, (uint16_T *) src + j, stride_x);
          break;

        default:
          mexErrMsgTxt("Unsupported type!");
      }
    }
  }

  mxFree(index);
  mxFree(size_Y);
  /* TODO: aren't we forgetting anything? */
}
