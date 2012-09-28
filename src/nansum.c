#include <matrix.h>
#include <mex.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "compiler.h"

#if defined (COMPILER_MSVC)
#include <math.h>
#define isnan _isnan
#define INFINITY (HUGE_VAL+HUGE_VAL)
#define NAN (INFINITY - INFINITY)
#elif defined(COMPILER_LCC)
#include <math.h>
#define INFINITY (DBL_MAX+DBL_MAX)
#define NAN (INFINITY - INFINITY)
#else
#include <math.h>
#endif

/* Shortcuts: */
#define size(x) mxGetDimensions(x)
#define ndim(x) mxGetNumberOfDimensions(x)

mwSize stride(int dim, int shape_len, const mwSize *shape)
{
  if (dim == shape_len - 1) return 1;
  if (dim < shape_len - 1)
    return shape[shape_len - 1] * stride(dim, shape_len - 1, shape);
  return -1; /* something odd happened */
}

void offset_to_index(mwSize offset, int ndim, const mwSize *shape, 
  mwSize *index)
{
  int i = 0;
  for (i = 0; i < ndim; ++i) {
    index[i] = offset % stride(i, ndim, shape);
    offset -= offset / stride(i, ndim, shape);
  }
}

mwSize index_to_offset(int ndim, mwSize *index, const mwSize *shape)
{
  int i = 0;
  mwSize offset = 0;
  for (i = 0; i < ndim; ++i)
    offset += stride(i, ndim, shape) * index[i];
  return offset;
}

double sum(int n, double *x0, mwSize stride) 
{
  int i;
  double result = 0;

  printf("in sum with n=%d, x0=%p, stride=%d. ", n, x0, stride);
  for (i = 0; i < n; ++i) {
    printf("%f ", x0[i * stride]);
    if (!isnan(x0[i * stride]))
      result += x0[i * stride];
  }
  printf("result=%f\n", result);
  return result;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  const mxArray *X = prhs[0], *Y = plhs[0];
  const mwSize *dims;
  mwSize *index, *size_X, *size_Y, stride_x;
  mwIndex j;
  mwIndex indx;
  int i, numdims, squash_dim;
  int numelin, numelout, x0, x1, y1;
  mxClassID classid;

  /* Check and handle input */
  switch(nrhs) {
    case 0: 
      mexErrMsgTxt("Too few input arguments!");
      break;

    case 1: /* No dimension is given. Find first non-singleton dimension: */
      squash_dim = 0;
      for(i = ndim(X); i-- > 0;) {
        if(size(X)[i] > 1)
          squash_dim = i;
      }
      break;

    case 2: 
      squash_dim = mxGetScalar(prhs[1]);
      if(squash_dim <= 0 || squash_dim > ndim(X))
        mexErrMsgTxt("Invalid value for dimension (argument 2)!");
      squash_dim -= 1; /* convert form MATLAB to C indexing */
      break;

    default:
      mexErrMsgTxt("Too many input arguments!");
  }

  /* Create output array Y: */
  size_Y = mxMalloc(ndim(X) * sizeof(size_Y));
  memcpy(size_Y, size(X), ndim(X) * sizeof(size_Y));
  size_Y[squash_dim] = 1; 

  classid = mxGetClassID(X);  /* copy class ID from X */
  
  Y = plhs[0] = mxCreateNumericArray(
    ndim(X), size_Y, classid, mxGetImagData(X) != NULL);
  assert(Y == plhs[0]); 

  /* Store calls to stat function with offset and stride in Y: */
  index = (mwSize *) mxMalloc(ndim(X) * sizeof(mwSize));
  stride_x = stride(squash_dim, ndim(X), size(X));

  {
    double *dest = mxGetData(Y);
    double *src = mxGetData(X);
    printf("&X = %p.\n", src);

    /* FIXME: something is wrong in calculation of offset & stride, probably
     * C/Fortran confusion. MATLAB uses FORTRAN order. */

    for (i = 0; i < mxGetNumberOfElements(Y); ++i) {

      offset_to_index(i, ndim(Y), size(Y), index);

      printf("i = %d. index=[", i);
      for (j = 0; j < ndim(Y); ++j) printf("%d ", index[j]);
      printf("]", i);

      j = index_to_offset(ndim(X), index, size(X));

      printf(", j = %d. ", j);

      dest[i] = sum(size(X)[squash_dim], src + j, stride_x);
    }
  }

  mxFree(index);
  mxFree(size_Y);
}
