
#if defined (__BORLANDC__)
  #include <windows.h>
  #include "gettimeofday.h"

#elif defined (_MSC_VER)
  #include <windows.h>
  #include "gettimeofday.h"
  #include "stdint.h"

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
