#include "compiler.h"
#include "platform.h"
#include "platform_includes.h"

#ifdef PLATFORM_OSX
#include <Availability.h>
#if __MAC_OS_X_VERSION_MAX_ALLOWED < 101300

/*
 * OS X did not have clock_gettime for a long time.
 *
 * In OS X El Captain and the accompanying XCode with MacOSX10.12.sdk it is not available.
 * In macOS High Sierra and the accompanying XCode with MacOSX10.14.sdk it is available.
 *
 * This is a drop-in replacement based on https://gist.github.com/jbenet/1087739
 * that uses clock_get_time. The first argument to this function is CLOCK_REALTIME,
 * which gets ignored.
 *
 * See also https://github.com/zeromq/libzmq/issues/2175 where this is discussed.
 */

#include <mach/clock.h>
#include <mach/mach.h>

int clock_gettime(int ignore, struct timespec *ts) {
  clock_serv_t cclock;
  mach_timespec_t mts;
  host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
  clock_get_time(cclock, &mts);
  mach_port_deallocate(mach_task_self(), cclock);
  ts->tv_sec  = mts.tv_sec;
  ts->tv_nsec = mts.tv_nsec;
  return 0;
}
#endif /* __MAC_OS_X_VERSION_MAX_ALLOWED */

#elif defined (COMPILER_MINGW_ORG)

/* MinGW does not have clock_gettime, make a drop-in replacement that uses clock_get_time */
/* the first argument to this function is CLOCK_REALTIME, which gets ignored */

LARGE_INTEGER
getFILETIMEoffset()
{
    SYSTEMTIME s;
    FILETIME f;
    LARGE_INTEGER t;

    s.wYear = 1970;
    s.wMonth = 1;
    s.wDay = 1;
    s.wHour = 0;
    s.wMinute = 0;
    s.wSecond = 0;
    s.wMilliseconds = 0;
    SystemTimeToFileTime(&s, &f);
    t.QuadPart = f.dwHighDateTime;
    t.QuadPart <<= 32;
    t.QuadPart |= f.dwLowDateTime;
    return (t);
}

int
clock_gettime(int ignore, struct timeval *tv)
{
    LARGE_INTEGER           t;
    FILETIME                f;
    double                  microseconds;
    static LARGE_INTEGER    offset;
    static double           frequencyToMicroseconds;
    static int              initialized = 0;
    static BOOL             usePerformanceCounter = 0;

    if (!initialized) {
        LARGE_INTEGER performanceFrequency;
        initialized = 1;
        usePerformanceCounter = QueryPerformanceFrequency(&performanceFrequency);
        if (usePerformanceCounter) {
            QueryPerformanceCounter(&offset);
            frequencyToMicroseconds = (double)performanceFrequency.QuadPart / 1000000.;
        } else {
            offset = getFILETIMEoffset();
            frequencyToMicroseconds = 10.;
        }
    }
    if (usePerformanceCounter) QueryPerformanceCounter(&t);
    else {
        GetSystemTimeAsFileTime(&f);
        t.QuadPart = f.dwHighDateTime;
        t.QuadPart <<= 32;
        t.QuadPart |= f.dwLowDateTime;
    }

    t.QuadPart -= offset.QuadPart;
    microseconds = (double)t.QuadPart / frequencyToMicroseconds;
    t.QuadPart = microseconds;
    tv->tv_sec = t.QuadPart / 1000000;
    tv->tv_usec = t.QuadPart % 1000000;
    return (0);
}

#endif

/*
#include <stdio.h>

int main(int argc, char **argv) {

  struct timespec ts;
  clock_gettime(0, &ts);
  printf("s:  %lu\n", ts.tv_sec);
  printf("ns: %lu\n", ts.tv_nsec);
  return 0;
}
*/
