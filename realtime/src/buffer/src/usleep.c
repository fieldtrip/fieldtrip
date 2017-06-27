
#ifdef PLATFORM_WINDOWS
/*
 * work around lack of usleep on MinGW, see
 * https://searchcode.com/codesearch/view/20805233/
*/

#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif

#include <unistd.h>
#include <windows.h>
#include <sys/types.h>
#include <errno.h>

int __cdecl usleep (useconds_t);

int __cdecl
usleep (useconds_t us)
{
  if (us >= 1000000)
    return EINVAL;
  else if (us != 0)
    Sleep (us / 1000);

  return 0;
}

#endif
