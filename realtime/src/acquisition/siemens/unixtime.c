#include <unixtime.h>

#ifdef WIN32

#define WIN32_LEAN_AND_MEAN
#include <windows.h>

/* WIN32 epoch starts in 1601, Unix in 1970 */

#define EPOCHFILETIME (116444736000000000LL)

double getUnixTimeAsDouble() {
    FILETIME        ft;
    LARGE_INTEGER   li;
    __int64         t;

    GetSystemTimeAsFileTime(&ft);
	
    li.LowPart  = ft.dwLowDateTime;
    li.HighPart = ft.dwHighDateTime;
	
    t  = li.QuadPart;       /* In 100-nanosecond intervals */
    t -= EPOCHFILETIME;     /* Offset to the Epoch time */
	
    t /= 10;                /* In microseconds */
	return (double) t * 1e-6;
}

#else

#include <sys/time.h>

double getUnixTimeAsDouble() {
	struct timeval tv;
	
	gettimeofday(&tv,0);
	return (double) tv.tv_sec + (double) tv.tv_usec*1e-6;
}

#endif

