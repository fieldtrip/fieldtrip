#ifdef __BORLANDC__
    #include <windows.h>
#else
  #ifdef __MINGW32__
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
#endif

#ifndef WIN32
  #define closesocket(s) (close(s))
#endif

