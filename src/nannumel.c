#include "nanaccum.h"

/* What is C a horrible language. To hoops we have to get through to make this
 * work for different data types... Since overloading does not work, we
 * generate type-specific functions: */
#define fname(name, suffix) name ## _ ## suffix
#define nanstat_template(TYPE)\
double fname(nanstat, TYPE)(int n, TYPE *x0, mwSize stride) \
{\
  int i; int count = 0;\
  for (i = 0; i < n; ++i) {\
    if (!isnan(x0[i * stride]))\
      count += 1;\
  }\
  return count;\
}

nanstat_template(float);
nanstat_template(double);
nanstat_template(int8_T); nanstat_template(uint8_T);
nanstat_template(int16_T); nanstat_template(uint16_T);
nanstat_template(int32_T); nanstat_template(uint32_T);
nanstat_template(int64_T); nanstat_template(uint64_T);

#include "nanaccum.c" /* takes care of generic nan-handling of arrays :D. */
