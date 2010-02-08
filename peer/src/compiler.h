/* prevent double include */
#ifndef COMPILER_H
#define COMPILER_H

#if    defined (__BORLANDC__)
#define COMPILER_BORLAND

#elif defined (_MSC_VER)
#define COMPILER_MSVC

#elif defined (__CYGWIN32__)
#define COMPILER_CYGWIN

#elif defined (__MINGW32__)
#define COMPILER_MINGW

#else
#error "Unknown compiler - please report the compiler details to http://fieldtrip.fcdonders.nl"

#endif

#endif /* COMPILER_H */
