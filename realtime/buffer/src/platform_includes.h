#include "platform.h"
#include "compiler.h"

#if defined (COMPILER_BORLAND)

#endif

/* these platforms will always use a similar gcc compiler */
#if defined (PLATFORM_LINUX) || defined (PLATFORM_OSX)
	#include <sys/types.h>
	#include <sys/socket.h>
	#include <sys/un.h>
	#include <netinet/in.h>
	#include <arpa/inet.h>
	#include <netdb.h>
	#include <unistd.h>
	#include <strings.h>
	#include <stdint.h>

	#define closesocket(s) (close(s))
	
#elif defined (PLATFORM_WIN64)

#if defined (COMPILER_MSVC)
	#include <winsock2.h>						/* for timeval */
	#include "win32/gettimeofday.h"
	#include "win32/stdint.h"

	#define bzero(b,len) (memset((b), '\0', (len)), (void) 0)
	#define usleep(x)    (Sleep((x)/1000))
	#define strcasecmp(a,b) (strcmpi(a,b))
#elif defined (COMPILER_MINGW)
	#include <winsock2.h>
	#include <sys/time.h>
	#include <stdint.h>

	#define bzero(b,len) memset(b,0,len)
	#define usleep(x)    (Sleep((x)/1000))
		
	#ifndef strcasecmp
		#define strcasecmp(a,b) (strcmpi(a,b))
	#endif
#else
	#error "Unsupported compiler"
#endif

#elif defined (PLATFORM_WIN32)
/* there are various compiler options for windows */
#if defined (COMPILER_BORLAND)
	#include <windows.h>
	#include "win32/gettimeofday.h"

    #define bzero(b,len) (memset((b), '\0', (len)), (void) 0)
	#define usleep(x)    (Sleep((x)/1000))
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

#elif defined (COMPILER_MSVC)
	#include <winsock2.h>						/* for timeval */
	#include "win32/gettimeofday.h"
	#include "win32/stdint.h"

	#define bzero(b,len) (memset((b), '\0', (len)), (void) 0)
	#define usleep(x)    (Sleep((x)/1000))
	#define strcasecmp(a,b) (strcmpi(a,b))

#elif defined (COMPILER_MINGW)
	#include <winsock2.h>
	#include <sys/time.h>
	#include <stdint.h>

	#define bzero(b,len) (memset((b), '\0', (len)), (void) 0)
	#define usleep(x)    (Sleep((x)/1000))
	#define strcasecmp(a,b) (strcmpi(a,b))

#elif defined (COMPILER_CYGWIN)
	#include <winsock2.h>

	#define strcasecmp(a,b) (strcmpi(a,b))
#elif defined (COMPILER_LCC)
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
	
/*	#define PTW32_STATIC_LIB 
	#define __cdecl
	#define PTW_CDECL
*/
#endif /* compiler */

#endif /* platform */
