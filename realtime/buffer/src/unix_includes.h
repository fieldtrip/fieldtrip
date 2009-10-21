#ifdef __BORLANDC__
  #include <windows.h>
  #include "gettimeofday.h"
#else
  #include <unistd.h>
  #include <strings.h>
#endif
#ifdef WIN32
  #define bzero(b,len) (memset((b), '\0', (len)), (void) 0)
  #define usleep(x)    (Sleep((x)/1000))
#endif
#if defined (__BORLANDC__) || defined (__WIN32__)
  #define strcasecmp(a,b) (strcmpi(a,b))
#endif
