#if defined (__BORLANDC__)
    #include <windows.h>

#elif defined (_MSC_VER)
    #include <windows.h>

#elif defined (__CYGWIN32__)
    #include <windows.h>

#elif defined (__MINGW32__)
    #include <winsock2.h>
    #include <win32compat.h>
    #ifndef WIN32
      #define WIN32 1
    #endif

#else
    #include <sys/types.h>
    #include <sys/socket.h>
    #include <netinet/in.h>
    #include <arpa/inet.h>
    #include <netdb.h>
#endif

#ifndef WIN32
  #define closesocket(s) (close(s))
#endif

