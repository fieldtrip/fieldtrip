
/* prevent double include */
#ifndef PLATFORM_H
#define PLATFORM_H

#if defined(linux) || defined(__linux) || defined(__linux__) || defined(__GNU__) || defined(__GLIBC__) 
/* linux, also other platforms (Hurd etc) that use GLIBC */
#define PLATFORM_LINUX

#elif defined(__FreeBSD__) || defined(__NetBSD__) || defined(__OpenBSD__) || defined(__DragonFly__)
/* BSD */
#define PLATFORM_BSD

#elif defined(sun) || defined(__sun)
/* Solaris */
#define PLATFORM_SUN

#elif defined(__sgi)
/* SGI Irix */
#define PLATFORM_SGI

#elif defined(__hpux)
/* hp unix */
#define PLATFORM_HP

#elif defined(__CYGWIN__)
/* cygwin is not win32 */
#define PLATFORM_CYGWIN

#elif defined(_WIN64) || defined(__WIN64__) || defined(WIN64)
/* win64 */
#define PLATFORM_WIN64
#define PLATFORM_WINDOWS

#elif defined(_WIN32) || defined(__WIN32__) || defined(WIN32)
/* win32 */
#define PLATFORM_WIN32
#define PLATFORM_WINDOWS

#elif defined(__BEOS__)
/* BeOS */
#define PLATFORM_BEOS

#elif defined (__APPLE__) && defined (__MACH__)
/* MacOSX */
#define PLATFORM_OSX

#elif defined(macintosh) || defined(__APPLE__) || defined(__APPLE_CC__)
/* MacOS classic */
#define PLATFORM_MAC

#elif defined(__IBMCPP__) || defined(_AIX)
/* IBM */
#define PLATFORM_AIX

#elif defined(__amigaos__)
/* AmigaOS */
#define PLATFORM_AMIGA

#elif defined(__QNXNTO__)
/* QNX */
#define PLATFORM_QNX

#elif defined(__VXWORKS__)
/* vxWorks */
#define PLATFORM_VXWORKS

#elif defined(unix) || defined(__unix) || defined(_XOPEN_SOURCE) || defined(_POSIX_SOURCE)
/* generic unix platform */
#define PLATFORM_UNIX

#else
/* the platform cannot be determined at compile time */
#error "Unknown platform - please report the platform details to http://www.fieldtriptoolbox.org"
#endif

#endif /* PLATFORM_H */

