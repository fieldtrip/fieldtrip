#include <time.h>

#ifdef __cplusplus
extern "C" {
#endif

#define CLOCK_REALTIME 0

/* OS X does not have clock_gettime, make a drop-in replacement that uses clock_get_time */
/* the first argument to this function is CLOCK_REALTIME, which gets ignored */
int clock_gettime(int ignore, struct timespec *tp);

#ifdef __cplusplus
}
#endif


