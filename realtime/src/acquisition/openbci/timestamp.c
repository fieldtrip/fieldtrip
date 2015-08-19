/* 

This is the snipper get_monotonic_time.c from
http://stackoverflow.com/questions/21665641/ns-precision-monotonic-clock-in-c-on-linux-and-os-x/21665642#21665642
which in turn is based on the snippet current_utc_time.c from
https://gist.github.com/jbenet/1087739

   OS X, compile with: gcc get_monotonic_time.c
  Linux, compile with: gcc get_monotonic_time.c -lrt
Windows, ??

*/

#include <time.h>
#include <sys/time.h>
#include <stdio.h>

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#include <unistd.h>
#endif

// Use clock_gettime in linux, clock_get_time in OS X.
void get_monotonic_time(struct timespec *ts){
#ifdef __MACH__
  clock_serv_t cclock;
  mach_timespec_t mts;
  host_get_clock_service(mach_host_self(), SYSTEM_CLOCK, &cclock);
  clock_get_time(cclock, &mts);
  mach_port_deallocate(mach_task_self(), cclock);
  ts->tv_sec = mts.tv_sec;
  ts->tv_nsec = mts.tv_nsec;
#else
  clock_gettime(CLOCK_MONOTONIC, ts);
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
  get_monotonic_time(&before);
  while (1) {
    usleep(100000);
    get_monotonic_time(&after);
    printf("deltaT = %e s\n", get_elapsed_time(&before,&after));
  }
}
*/
