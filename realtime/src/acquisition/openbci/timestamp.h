#ifndef _TIMESTAMP_H_
#define _TIMESTAMP_H_

void get_monotonic_time(struct timespec *ts);
double get_elapsed_time(struct timespec *before, struct timespec *after);

#endif
