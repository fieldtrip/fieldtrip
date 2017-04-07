#include <mex.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#ifdef _MSC_VER
  /* Microsoft doesn't define isnan, but has an _isnan instead. Since it
   * appears to be buggy, we simply do: */
  #define isnan(x) !(x == x)
  #define NAN sqrt(-1)
#endif


/* Shortcuts: */
#define size(x) mxGetDimensions(x)
#define ndim(x) mxGetNumberOfDimensions(x)
