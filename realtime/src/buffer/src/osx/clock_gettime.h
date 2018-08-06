#include <time.h>

#ifdef PLATFORM_OSX

/*
 * OS X did not have clock_gettime for a long time, whereas in the most recent
 * macOS High Sierra and the accompanying XCode with MacOSX10.12.sdk it is available.
 *
 * This is a drop-in replacement based on https://gist.github.com/jbenet/1087739
 * that uses clock_get_time. The first argument to this function is CLOCK_REALTIME,
 * which gets ignored.
 *
 * See also https://github.com/zeromq/libzmq/issues/2175 where this is discussed.
 */

#ifdef not(__CLOCK_AVAILABILITY)

#define CLOCK_REALTIME 0

#ifdef __cplusplus
extern "C" {
#endif

int clock_gettime(int ignore, struct timespec *tp);

#ifdef __cplusplus
}
#endif

#endif /* if __CLOCK_AVAILABILITY */
#endif /* if PLATFORM_OSX */
