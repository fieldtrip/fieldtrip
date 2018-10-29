#include <time.h>

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

#define CLOCK_REALTIME 0

#ifdef __cplusplus
extern "C" {
#endif

int clock_gettime(int ignore, struct timespec *tp);

#ifdef __cplusplus
}
#endif

#endif /* __MAC_OS_X_VERSION_MAX_ALLOWED */
#endif /* if PLATFORM_OSX */
