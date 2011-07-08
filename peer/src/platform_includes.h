#ifndef PLATFORM_INCLUDES_H
#define PLATFORM_INCLUDES_H

#include "platform.h"
#include "compiler.h"

/* these platforms will always use a similar gcc compiler */
#if defined (PLATFORM_LINUX) || defined (PLATFORM_OSX)
#include <sys/socket.h>
#include <sys/resource.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <netdb.h>
#include <ifaddrs.h>
#include <pwd.h>           /* for getpwuid and geteuid */

#include <sys/types.h>
#include <sys/socket.h>
#include <sys/un.h>

/* closesocket exists on windows and is used in the code */
#define closesocket(s) (close(s))

#elif defined (PLATFORM_WIN32) || defined(PLATFORM_WIN64)
/* there are various compiler options for windows */

#if defined (COMPILER_BORLAND)
#include <windows.h>
#include "win32/gettimeofday.h"
/* without the following, compilation with the Borland command line tools fails -- SK */
typedef __int8            int8_t;
typedef __int16           int16_t;
typedef __int32           int32_t;
typedef __int64           int64_t;
typedef unsigned __int8   uint8_t;
typedef unsigned __int16  uint16_t;
typedef unsigned __int32  uint32_t;
typedef unsigned __int64  uint64_t;
#define random()          (random(INT32_MAX))
#define bzero(b,len)      (memset((b), '\0', (len)), (void) 0)
#define usleep(x)         (Sleep((x)/1000))
#define strcasecmp(a,b)   (strcmpi(a,b))
#define strncasecmp(a,b,n)(strncmpi(a,b,n))

#elif defined (COMPILER_MSVC)
#define WIN32_MEAN_AND_LEAN /* otherwise multiple includes of windows.h fail */
#define __CLEANUP_CXX  /* for pthreads */
#include <winsock2.h>  /* this includes windows.h */
#include <ws2tcpip.h>
#include <stdio.h>
#include "win32/stdint.h"
#include "win32/gettimeofday.h"
#pragma comment (lib, "Ws2_32.lib")
#define bzero(b,len) (memset((b), '\0', (len)), (void) 0)
#define usleep(x)    (Sleep((x)/1000))
#define sleep(x)     (Sleep((x)*1000))
#define strcasecmp strcmpi
#define strncasecmp strnicmp

#elif defined (COMPILER_MINGW)
#include <ws2tcpip.h>
#include <windows.h>
#include <stdint.h>
//#include <win32compat.h>
#define bzero(b,len) (memset((b), '\0', (len)), (void) 0)
#define sleep(x)     (Sleep((x)*1000))
#define usleep(x)    (Sleep((x)/1000))

#elif defined (COMPILER_CYGWIN)
#include <windows.h>

#endif /* compiler */
#endif /* platform */

#if defined(PLATFORM_WIN32) || defined(PLATFORM_WIN64)
#if defined(COMPILER_MINGW)
	#include "stdint.h"
#else
	#include "win32/stdint.h"
#endif
#else
	#include <stdint.h>
	#include <syslog.h>
#endif

#endif /* PLATFORM_INCLUDES_H */
