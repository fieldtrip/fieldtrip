/* 
 * Copyright (C) 2008, Robert Oostenveld & Christian Hesse
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 */

#ifndef PLATFORM_INCLUDES_H
#define PLATFORM_INCLUDES_H

#include "platform.h"
#include "compiler.h"

#if defined (PLATFORM_OSX)
    /* Linux and OSX use a similar gcc compiler */
    #include <sys/types.h>
    #include <sys/socket.h>
    #include <sys/un.h>
    #include <netinet/in.h>
    #include <netinet/tcp.h>
    #include <arpa/inet.h>
    #include <netdb.h>
    #include <unistd.h>
    #include <strings.h>
    #include <stdint.h>

    #define closesocket(s) (close(s))

    /* OSX needs an alternative for clock_gettime */
    #include "osx/clock_gettime.h"

#elif defined (PLATFORM_LINUX)
    /* Linux and OSX use a similar gcc compiler */
    #include <sys/types.h>
    #include <sys/socket.h>
    #include <sys/un.h>
    #include <netinet/in.h>
    #include <netinet/tcp.h>
    #include <arpa/inet.h>
    #include <netdb.h>
    #include <unistd.h>
    #include <strings.h>
    #include <stdint.h>

    #define closesocket(s) (close(s))

#elif defined (PLATFORM_WIN64) && defined (COMPILER_MSVC)
    #include <winsock2.h>                       /* for timeval */
    #include "win32/gettimeofday.h"
    #include "win32/stdint.h"

    #define bzero(b,len)    (memset((b), '\0', (len)), (void) 0)
    #define usleep(x)       (Sleep((x)/1000)) // should be replaced by win32/usleep.h
    #define strcasecmp(a,b) (strcmpi(a,b))

#elif defined (PLATFORM_WIN64) && defined (COMPILER_MINGW_W64)
    #include <winsock2.h>
    #include <windows.h>
    #include <sys/time.h>
    #include <stdint.h>
    #include "win32/fsync.h"
    #include "win32/usleep.h"

    #define bzero(b,len)            (memset((b), '\0', (len)), (void) 0)
    #define bcopy(src, dest, len)   (memmove(dest, src, len))

#elif defined (PLATFORM_WIN32) && defined (COMPILER_BORLAND)
    #include <windows.h>
    #include "win32/gettimeofday.h"

    #define bzero(b,len)    (memset((b), '\0', (len)), (void) 0)
    #define usleep(x)       (Sleep((x)/1000))
    #define strcasecmp(a,b) (strcmpi(a,b))
    /* without the following, compilation with the Borland command line tools fails -- SK */
    typedef __int8            int8_t;
    typedef __int16           int16_t;
    typedef __int32           int32_t;
    typedef __int64           int64_t;
    typedef unsigned __int8   uint8_t;
    typedef unsigned __int16  uint16_t;
    typedef unsigned __int32  uint32_t;
    typedef unsigned __int64  uint64_t;

#elif defined (PLATFORM_WIN32) && defined (COMPILER_MSVC)
    #include <winsock2.h>                       /* for timeval */
    #include "win32/gettimeofday.h"
    #include "win32/stdint.h"

    #define bzero(b,len)    (memset((b), '\0', (len)), (void) 0)
    #define usleep(x)       (Sleep((x)/1000))  // should be replaced by win32/usleep.h
    #define strcasecmp(a,b) (strcmpi(a,b))

#elif defined (PLATFORM_WIN32) && defined (COMPILER_MINGW_W64)
    #include <winsock2.h>
    #include <windows.h>
    #include <sys/time.h>
    #include <stdint.h>
    #include "win32/fsync.h"
    #include "win32/usleep.h"

    #define bzero(b,len)            (memset((b), '\0', (len)), (void) 0)
    #define bcopy(src, dest, len)   (memmove(dest, src, len))

#elif defined (PLATFORM_WIN32) && defined (COMPILER_MINGW_ORG)
    #include <winsock2.h>
    #include <windows.h>
    #include <sys/time.h>
    #include <stdint.h>
    #include "win32/clock_gettime.h"

    #define bzero(b,len)    (memset((b), '\0', (len)), (void) 0)
    #define usleep(x)       (Sleep((x)/1000)) // should be replaced by win32/usleep.h
    #define strcasecmp(a,b) (strcmpi(a,b))

#elif defined (PLATFORM_CYGWIN) && defined (COMPILER_CYGWIN)
    //  #include <winsock2.h>
    //  #include <windows.h>
    #include <sys/socket.h>
    #include <sys/un.h>
    #include <netinet/in.h>
    #include <netinet/tcp.h>
    #include <netinet/ip.h>
    #include <unistd.h>  /* for close() */
    #include <netdb.h>

    #define closesocket(s) (close(s))

#elif defined (PLATFORM_WIN32) && defined (COMPILER_LCC)
    #include <winsock2.h>
    #include <windows.h>
    #include "win32/gettimeofday.h"

    #define strcasecmp(a,b) (strcmpi(a,b))

    #define bzero(b,len) (memset((b), '\0', (len)), (void) 0)
    #define usleep(x)    (Sleep((x)/1000))

    #ifndef UINT8_T
        #define UINT8_T   unsigned char
    #endif

    #ifndef INT8_T
        #define INT8_T    char
    #endif

    #ifndef UINT16_T
        #define UINT16_T  unsigned short
    #endif

    #ifndef INT16_T
        #define INT16_T   short
    #endif

    #ifndef UINT32_T
        #define UINT32_T  unsigned int
    #endif

    #ifndef INT32_T
        #define INT32_T   int
    #endif

    #ifndef UINT64_T
        #define UINT64_T  unsigned long long
    #endif

    #ifndef INT64_T
        #define INT64_T   long long
    #endif

/*  #define PTW32_STATIC_LIB
    #define __cdecl
    #define PTW_CDECL
  */

#else
    /* it cannot be determined at compile time */
    #error "Unknown combination of platform and compiler - please report this to http://bugzilla.fieldtriptoolbox.org"

#endif

#endif /* PLATFORM_INCLUDES_H */
