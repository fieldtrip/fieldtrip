#include "nanaccum.h"

/* Helper function to find linearized distance between elements on a given
 * dimension of a multidimensional array. Lower indices are close in memory
 * (column-major or FORTRAN/MATLAB order: */
mwSize stride(int dim, int shape_len, const mwSize *shape)
{
  int i, s = 1;
  for (i = 0; i < dim; i++) {
    s *= shape[i];
  }
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
      squash_dim -= 1; /* convert from MATLAB to C indexing */
      break;

    default:
      mexErrMsgTxt("Too many input arguments!");
  }

  /* Create output array Y: */
  size_Y = mxMalloc(ndim(X) * sizeof(size_Y));
  memcpy(size_Y, size(X), ndim(X) * sizeof(size_Y));
  size_Y[squash_dim] = 1; 

  classid = mxGetClassID(X);  /* Copy class ID from X. */
  
  if (classid == mxSINGLE_CLASS) {
    /* Return single-precision if input is single precision. */
    Y = plhs[0] = mxCreateNumericArray(
      ndim(X), size_Y, mxSINGLE_CLASS, mxGetImagData(X) != NULL);
  }
  else {
    /* By default, use double precision. */
    Y = plhs[0] = mxCreateNumericArray(
      ndim(X), size_Y, mxDOUBLE_CLASS, mxGetImagData(X) != NULL);
  }

  /* Prepare variables for linear indexing: */
  index = (mwSize *) mxMalloc(ndim(X) * sizeof(mwSize));
  stride_x = stride(squash_dim, ndim(X), size(X));

  /* MATLAB's nansum supports out-of-range dims to operate on. */
  /* Find the number of elements in our squashed dimension: */
  if(squash_dim >= ndim(X))
    squash_len = 1; 
  else
    squash_len = size(X)[squash_dim];

  {
    void *src = mxGetData(X), *src_imag = mxGetImagData(X);
    double *dest = mxGetData(Y), *dest_imag = mxGetImagData(Y);

    for (i = 0; i < mxGetNumberOfElements(Y); ++i) {
      /* For each element in the output array, do: */

      /* Transform the element number in an index in Y. Please not that size()
       * uses mxGetNumberOfDimensions(), which *ignores trailing singleton
       * dimensions*. In that case, this elegant conversion does not work
       * anymore. But, since we know that ndim(Y) == ndim(X), we use that
       * instead:*/
      offset_to_index(i, ndim(X), size(Y), index);

      /* And map this index back to an offset in X: */
      j = index_to_offset(ndim(X), index, size(X));
      
      switch (classid) {
        /* Oh, the fucking horror of dumb, statically typed languages... */
        case mxINT8_CLASS:
          dest[i] = nanstat_int8_T(squash_len, (int8_T *) src + j, stride_x);
          break;

        case mxUINT8_CLASS:
          dest[i] = nanstat_uint8_T(squash_len, (uint8_T *) src + j, stride_x);
          break;

        case mxINT16_CLASS:
          dest[i] = nanstat_int16_T(squash_len, (int16_T *) src + j, stride_x);
          break;

        case mxUINT16_CLASS:
          dest[i] = nanstat_uint16_T(squash_len, (uint16_T *) src + j, 
            stride_x);
          break;

        case mxINT32_CLASS:
          dest[i] = nanstat_int32_T(squash_len, (int32_T *) src + j, stride_x);
          break;

        case mxUINT32_CLASS:
          dest[i] = nanstat_uint32_T(squash_len, (uint32_T *) src + j, 
            stride_x);
          break;

        case mxINT64_CLASS:
          dest[i] = nanstat_int64_T(squash_len, (int64_T *) src + j, stride_x);
          break;

        case mxUINT64_CLASS:
          dest[i] = nanstat_uint64_T(squash_len, (uint64_T *) src + j, 
            stride_x);
          break;

        case mxSINGLE_CLASS:
          ((float *) dest)[i] = (float) nanstat_float(
            squash_len, (float *) src + j, stride_x);
          if (src_imag) {
            ((float *) dest_imag)[i] = (float) nanstat_float(
              squash_len, (float *) src_imag + j, stride_x);
          }
          break;

        case mxDOUBLE_CLASS:
          dest[i] = nanstat_double(squash_len, (double *) src + j, stride_x);
          if (src_imag) {
            dest_imag[i] = nanstat_double(squash_len, (double *) src_imag + j,                 stride_x);
          }
          break;

        case mxLOGICAL_CLASS:
          mexWarnMsgTxt(
            "The byte layout of logical arrays is not documented. Guessing.");
          dest[i] = nanstat_uint8_T(squash_len, (uint8_T *) src + j, stride_x);
          break;

        case mxCHAR_CLASS:
          mexWarnMsgTxt(
            "The byte layout of char arrays is not documented. Guessing.");
          dest[i] = nanstat_uint16_T(squash_len, (uint16_T *) src + j, 
            stride_x);
          break;

        default:
          mexErrMsgTxt("Unsupported type!");
      }
    }
  }

  if(index) mxFree(index);
  if(size_Y) mxFree(size_Y);
}
