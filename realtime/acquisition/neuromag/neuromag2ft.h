/*
 * Code by Gustavo Sudre (gsudre@cmu.edu)
 * based on the sample code by Lauri Parkkonen et al
 */

#include <fiff.h>

#ifndef FAIL
#define FAIL -1
#endif

#ifndef OK
#define OK 0
#endif

extern int fieldtrip_sock;
extern int process_tag(fiffTag tag);
extern void process_data(float **data, int nchan, int nsamp, fiffChInfo *ch_info);
extern void clean_up();
