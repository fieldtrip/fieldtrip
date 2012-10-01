#include "nanaccum.h"

/* What is C a horrible language. To hoops we have to get through to make this
 * work for different data types... Since overloading does not work, we
 * generate type-specific functions: */
#define fname(name, suffix) name ## _ ## suffix
#define nanstat_template(TYPE, INTERMEDIATE_TYPE)\
double fname(nanstat, TYPE)(int n, TYPE *x0, mwSize stride) \
{\
  int i; INTERMEDIATE_TYPE result = 0, c=0;\
  for (i = 0; i < n; ++i) {\
    if (!isnan(x0[i * stride])){\
      result += x0[i * stride];\
      c += 1;\
    }\
  }\
  return result / c;\
}

nanstat_template(float, float); /* Note that the calculations are performed with
                                  limited precision as well to be fully
                                  compatible with MATLABs nanmean. */

nanstat_template(double, double);
nanstat_template(int8_T, double); nanstat_template(uint8_T, double);
nanstat_template(int16_T, double); nanstat_template(uint16_T, double);
nanstat_template(int32_T, double); nanstat_template(uint32_T, double);
nanstat_template(int64_T, double); nanstat_template(uint64_T, double);

#include "nanaccum.c" /* takes care of generic nan-handling of arrays :D. */
