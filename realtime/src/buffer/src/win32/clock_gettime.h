#ifndef _CLOCK_GETTIME
#define _CLOCK_GETTIME
#include <time.h>

/* MinGW does not have clock_gettime, make a drop-in replacement that uses clock_get_time */
/* the first argument to this function is CLOCK_REALTIME, which gets ignored */

int clock_gettime(int ignore, struct timeval *tp);

#define CLOCK_REALTIME 0

#endif
