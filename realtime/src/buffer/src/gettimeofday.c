#include "platform.h"
#include "compiler.h"
#include "platform_includes.h"

/*
 * timeval.h    1.0 01/12/19
 *
 * Defines gettimeofday, timeval, etc. for Win32
 *
 * By Wu Yongwei
 *
 */


#ifdef PLATFORM_WINDOWS
#ifndef COMPILER_MINGW

#define EPOCHFILETIME ((INT64_T) 116444736000000000LL)

#ifdef COMPILER_LCC
VOID STDCALL GetSystemTimeAsFileTime(LPFILETIME);
#endif

int gettimeofday(struct timeval *tv, struct timezone *tz)
{
    FILETIME        ft;
    LARGE_INTEGER   li;
    INT64_T         t;
    static int      tzflag;

    if (tv) {
        GetSystemTimeAsFileTime(&ft);
        li.LowPart  = ft.dwLowDateTime;
        li.HighPart = ft.dwHighDateTime;
        t  = li.QuadPart;       /* In 100-nanosecond intervals */
        t -= EPOCHFILETIME;     /* Offset to the Epoch time */
        t /= 10;                /* In microseconds */
        tv->tv_sec  = (long)(t / 1000000);
        tv->tv_usec = (long)(t % 1000000);
    }

	#ifndef COMPILER_LCC
	/* LCC that comes with MATLAB has problems with _timezone and _daylight,
		and we don't need it anyway */
    if (tz) {
        if (!tzflag) {
            _tzset();
            tzflag++;
        }
        tz->tz_minuteswest = _timezone / 60;
        tz->tz_dsttime = _daylight;
    }
	#endif

    return 0;
}

#endif
#endif

