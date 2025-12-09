#include "mex.h"

// define EXTERN_C
#ifndef EXTERN_C
#ifdef __cplusplus
#define EXTERN_C extern "C"
#else
#define EXTERN_C
#endif
#endif

// existing MATRIX_DLL_EXPORT_SYM namespace hack (optional, keeps original code)
#if defined(__cplusplus) && defined(MATRIX_DLL_EXPORT_SYM)
  namespace matrix{ namespace detail{ namespace noninlined{ namespace mx_array_api{
#endif

EXTERN_C mxArray* mxSerialize(mxArray const *);
EXTERN_C mxArray* mxDeserialize(const void *, size_t);

#if defined(__cplusplus) && defined(MATRIX_DLL_EXPORT_SYM)
  }}}}
  using namespace matrix::detail::noninlined::mx_array_api;
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nlhs && nrhs) {
        plhs[0] = (mxArray *) mxSerialize(prhs[0]);
    }
}
