#include "compiler.h"
#include "platform.h"
#include "platform_includes.h"

/*
 * work around lack of clock_gettime in os x
 * https://gist.github.com/jbenet/1087739
 */

#ifdef PLATFORM_OSX

#include <mach/clock.h>
#include <mach/mach.h>

/* OS X does not have clock_gettime, make a drop-in replacement that uses clock_get_time */
/* the first argument to this function is CLOCK_REALTIME, which gets ignored */

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
