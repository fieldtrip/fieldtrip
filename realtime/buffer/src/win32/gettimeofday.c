/*
 * timeval.h    1.0 01/12/19
 *
 * Defines gettimeofday, timeval, etc. for Win32
 *
 * By Wu Yongwei
 *
 */

#ifndef _TIMEVAL_H
#define _TIMEVAL_H

#if TIME_WITH_SYS_TIME
#include <sys/time.h>
#include <time.h>
#elif HAVE_SYS_TIME_H
#include <sys/time.h>
#else
#include <time.h>
#endif

#if ! HAVE_GETTIMEOFDAY || MINGW

#if __WIN32__
#include <windows.h>
#else
#include <errno.h>
#endif

#ifndef __GNUC__
#define EPOCHFILETIME (116444736000000000i64)
#else
#define EPOCHFILETIME (116444736000000000LL)
#endif

struct timezone {
    int tz_minuteswest; /* minutes W of Greenwich */
    int tz_dsttime;     /* type of dst correction */
};

/*!
    \brief A Windows gettimeofday implementation.
 */

int gettimeofday(struct timeval *tv, struct timezone *tz)
{
    #if __WIN32__
    FILETIME        ft;
    LARGE_INTEGER   li;
    __int64         t;
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
    
    if (tz) {
        if (!tzflag) {
            _tzset();
            tzflag++;
        }
        tz->tz_minuteswest = _timezone / 60;
        tz->tz_dsttime = _daylight;
    }
    
    return 0;
    #else
    errno = ENOSYS;
    return -1;
    #endif
}
#endif
#endif


