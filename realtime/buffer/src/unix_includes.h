
#if defined (__BORLANDC__)
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
#elif defined (_MSC_VER)
  #include <winsock2.h>						/* for timeval */
  #include "win32/gettimeofday.h"
  #include "win32/stdint.h"
#elif defined (__MINGW32__)
  /* #include <windows.h> */
  #include <sys/time.h>
  #include <stdint.h>
#else
  #include <unistd.h>
  #include <strings.h>
  #include <stdint.h>
#endif

#if defined (__MINGW32__) || defined (_MSC_VER)
  #define bzero(b,len) (memset((b), '\0', (len)), (void) 0)
  #define usleep(x)    (Sleep((x)/1000))
#endif

#if defined (__BORLANDC__) || defined (__WIN32__)
  #define strcasecmp(a,b) (strcmpi(a,b))
#endif
