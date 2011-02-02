#ifndef LABVIEWDLL_H
#define LABVIEWDLL_H

#include <stdint.h>

#ifdef WIN32
#include <windows.h>
#else
/* Linux, OS X (?) */
typedef void* HANDLE;
#endif

typedef HANDLE (*OPEN_DRIVER_ASYNC_T)(void);
typedef int (*USB_WRITE_T)(HANDLE handle, char *data);
typedef int (*READ_MULTIPLE_SWEEPS_T)(HANDLE handle, char *data, intptr_t nBytesToRead);
typedef int (*READ_POINTER_T)(HANDLE handle, intptr_t *pointer);
typedef int (*CLOSE_DRIVER_ASYNC_T)(HANDLE handle);

#endif



