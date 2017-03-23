#ifndef _USLEEP
#define _USLEEP

/*
 * work around lack of usleep on MinGW, see
 * https://searchcode.com/codesearch/view/20805233/
*/

#include <unistd.h>
#include <windows.h>
#include <sys/types.h>
#include <errno.h>

int usleep (useconds_t);

#endif
