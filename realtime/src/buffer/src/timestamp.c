/*
 * This is derived from
 *  http://stackoverflow.com/questions/21665641/ns-precision-monotonic-clock-in-c-on-linux-and-os-x/21665642#21665642
 *  https://gist.github.com/jbenet/1087739
 *
 * OS X
 *   compile with gcc get_monotonic_time.c
 *   from http://stackoverflow.com/questions/11680461/monotonic-clock-on-osx
 *   SYSTEM_CLOCK returns the time since boot time
 *   CALENDAR_CLOCK returns the UTC time since 1970-01-01
 *
 * Linux
 *   compile with gcc get_monotonic_time.c -lrt
 *   from http://man7.org/linux/man-pages/man2/clock_gettime.2.html
 *   CLOCK_REALTIME represents seconds and nanoseconds since the Epoch
 *   CLOCK_MONOTONIC represents monotonic time since some unspecified starting point
 *   CLOCK_BOOTTIME identical to CLOCK_MONOTONIC, except it also includes any time that the system is suspended
 *
 * Windows ??
 */

#include "platform.h"
#include "compiler.h"

#if defined(PLATFORM_OSX)
#include <mach/clock.h>
#include <mach/mach.h>
#include <unistd.h>

#elif defined(PLATFORM_LINUX)
#include <time.h>
#include <sys/time.h>

#ifndef CLOCK_BOOTTIME
/* should be present since Linux 2.6.39, but is missing on Raspbian */
#define CLOCK_BOOTTIME CLOCK_MONOTONIC
#endif

#elif defined(PLATFORM_WINDOWS)
#include <Windows.h>
#include <Winbase.h>

#define CLOCK_BOOTTIME 1
#define CLOCK_REALTIME 2

struct timespec { long tv_sec; long tv_nsec; };    //header part

int clock_gettime(int unused, struct timespec *spec)      //C-file part
{
  __int64 wintime;
  GetSystemTimeAsFileTime((FILETIME*)&wintime);
  wintime      -=116444736000000000;           //1jan1601 to 1jan1970
  spec->tv_sec  =wintime / 10000000;           //seconds
  spec->tv_nsec =wintime % 10000000 *100;      //nano-seconds
  return 0;
}

#endif

#include <stdio.h>
#include "timestamp.h"

void get_monotonic_time(struct timespec *ts, int timestamp_ref){
#if defined(PLATFORM_OSX)
  clock_serv_t cclock;
  mach_timespec_t mts;
  switch (timestamp_ref) {
    case TIMESTAMP_REF_BOOT:
      host_get_clock_service(mach_host_self(), SYSTEM_CLOCK, &cclock);
      break;
    case TIMESTAMP_REF_EPOCH:
      host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
      break;
  }
  clock_get_time(cclock, &mts);
  mach_port_deallocate(mach_task_self(), cclock);
  ts->tv_sec = mts.tv_sec;
  ts->tv_nsec = mts.tv_nsec;
#elif defined(PLATFORM_LINUX)
  switch (timestamp_ref) {
    case TIMESTAMP_REF_BOOT:
      clock_gettime(CLOCK_BOOTTIME, ts);
      break;
    case TIMESTAMP_REF_EPOCH:
      clock_gettime(CLOCK_REALTIME, ts);
      break;
  }
#elif defined(PLATFORM_WINDOWS) && defined(COMPILER_MINGW)
  switch (timestamp_ref) {
    case TIMESTAMP_REF_BOOT:
      clock_gettime(CLOCK_BOOTTIME, ts);
      break;
    case TIMESTAMP_REF_EPOCH:
      clock_gettime(CLOCK_REALTIME, ts);
      break;
  }
#endif
}

double get_elapsed_time(struct timespec *before, struct timespec *after){
  double deltat_s  = after->tv_sec - before->tv_sec;
  double deltat_ns = after->tv_nsec - before->tv_nsec;
  return deltat_s + deltat_ns*1e-9;
}

/*
int main(){
  struct timespec before, after;
  get_monotonic_time(&before, TIMESTAMP_REF_EPOCH);
  printf("time since epoch = %10ld s, %10ld ns\n", before.tv_sec, before.tv_nsec);
  get_monotonic_time(&before, TIMESTAMP_REF_BOOT);
  printf("time since boot  = %10ld s, %10ld ns\n", before.tv_sec, before.tv_nsec);
  while (1) {
    get_monotonic_time(&after, TIMESTAMP_REF_BOOT);
    printf("deltaT = %f s\n", get_elapsed_time(&before,&after));
    usleep(100000);
  }
}
*/
