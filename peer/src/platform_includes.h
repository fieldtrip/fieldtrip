#include "platform.h"
#include "compiler.h"

/* these platforms will always use a similar gcc compiler */
#if defined (PLATFORM_LINUX) || defined (PLATFORM_OSX)
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <netdb.h>
#include <ifaddrs.h>
#include <pwd.h>           /* for getpwuid and geteuid */

/* closesocket exists on windows and is used in the code */
#define closesocket(s) (close(s))

#elif defined (PLATFORM_WIN32)

/* there are various compiler options for windows */
#if defined (COMPILER_BORLAND)
#include <windows.h>
#include "win32/gettimeofday.h"
#define random()     (random(INT32_MAX))
#define bzero(b,len) (memset((b), '\0', (len)), (void) 0)
#define usleep(x)    (Sleep((x)/1000))

#elif defined (COMPILER_MSVC)
#include <windows.h>
#define bzero(b,len) (memset((b), '\0', (len)), (void) 0)
#define usleep(x)    (Sleep((x)/1000))
#include "win32/gettimeofday.h"

#elif defined (COMPILER_MINGW)
#include <windows.h>
#include <win32compat.h>
#define bzero(b,len) (memset((b), '\0', (len)), (void) 0)
#define usleep(x)    (Sleep((x)/1000))

#elif defined (COMPILER_CYGWIN)
#include <windows.h>

#endif /* compiler */

#endif /* platform */
