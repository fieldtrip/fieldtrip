#ifndef _TIMESTAMP_H_
#define _TIMESTAMP_H_

#define TIMESTAMP_REF_BOOT  1
#define TIMESTAMP_REF_EPOCH 2

void get_monotonic_time(struct timespec *ts, int timestamp_ref);
double get_elapsed_time(struct timespec *before, struct timespec *after);

#endif
