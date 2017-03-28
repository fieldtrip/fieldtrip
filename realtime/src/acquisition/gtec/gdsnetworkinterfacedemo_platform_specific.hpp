//COPYRIGHT © 2013 G.TEC MEDICAL ENGINEERING GMBH, AUSTRIA
#ifndef GDSNETWORKINTERFACE_PLATFORM_SPECIFIC_HPP_INCLUDED
#define GDSNETWORKINTERFACE_PLATFORM_SPECIFIC_HPP_INCLUDED

// these are from FieldTrip/realtime/buffer
#include "platform.h"
#include "compiler.h"
#include "platform_includes.h"

// Linux specific stuff here ...
//-------------------------------------------------------------------------------------
#if defined(PLATFORM_LINUX) || defined(PLATFORM_OSX)

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <arpa/inet.h>
#include <pthread.h>
#include <time.h>
#include <unistd.h>
#include <stdexcept>

typedef int socketd_t;
typedef pthread_t thread_t;
typedef pthread_mutex_t mutex_t;

typedef struct _event_t
{
	pthread_mutex_t mutex_;
	pthread_cond_t cond_;
	bool signaled_;
} event_t;

#define SOCKET_ERROR (-1)
#define TCP_RECV_BUF_SIZE (64*1024)

#define MKR_WRITE( descriptor, buffer, size ) write( descriptor, buffer, size )
#define MKR_READ( descriptor, buffer, size ) read( descriptor, buffer, size )
#define MKR_ACCEPT( descriptor, ptr_address, address_length ) accept( descriptor, ( struct sockaddr* ) ptr_address, address_length )
#define MKR_CLOSE_SOCKET( descriptor ) close( descriptor )
#define MKR_SLEEP( t ) sleep( t )

#define MKR_MUTEX_INIT( mutex ) pthread_mutex_init( &mutex, 0 )
#define MKR_MUTEX_LOCK( mutex ) pthread_mutex_lock( &mutex )
#define MKR_MUTEX_UNLOCK( mutex ) pthread_mutex_unlock( &mutex )
#define MKR_MUTEX_DESTROY( mutex ) pthread_mutex_destroy( &mutex )

#define MKR_THREAD_CREATE( handle, start_routine, data ) CreateThread( handle, start_routine, data)
#define MKR_THREAD_EXIT() pthread_exit(0)
#define MKR_THREAD_JOIN( handle ) pthread_join( handle, 0 )
#define MKR_THREAD_DESTROY( handle )

#define THREAD_CALLBACK_RETURN_TYPE void*
#define THREAD_CALLBACK_CALLING_CONVENTION

extern "C" {
	typedef THREAD_CALLBACK_RETURN_TYPE (THREAD_CALLBACK_CALLING_CONVENTION *StartRoutine)(void*);
}

inline bool CreateThread(thread_t* handle, StartRoutine start_routine, void* arg)
{
    if (handle == 0)
        throw std::invalid_argument("Parameter 'handle' was NULL.");

    return (pthread_create(handle, 0, start_routine, arg) == 0);
}

inline void InitializeEvent(event_t *ev)
{
	pthread_mutex_init(&ev->mutex_, 0);
	pthread_cond_init(&ev->cond_, 0);
	ev->signaled_ = false;
}

inline void SetEvent(event_t *ev)
{
	pthread_mutex_lock(&ev->mutex_);
	ev->signaled_ = true;
	pthread_mutex_unlock(&ev->mutex_);
	pthread_cond_signal(&ev->cond_);
}

inline void ResetEvent(event_t *ev)
{
	pthread_mutex_lock(&ev->mutex_);
	ev->signaled_ = false;
	pthread_mutex_unlock(&ev->mutex_);
}

inline bool WaitForEvent(event_t *ev, unsigned int timeoutMilliseconds)
{
	struct timespec absoluteTimeout;
	clock_gettime(CLOCK_REALTIME, &absoluteTimeout);
	absoluteTimeout.tv_sec += timeoutMilliseconds / 1000;
	absoluteTimeout.tv_nsec += (timeoutMilliseconds % 1000) * 1e6;

	bool timeout = false;
	pthread_mutex_lock(&ev->mutex_);
	while (!ev->signaled_ && !timeout)
		timeout = (pthread_cond_timedwait(&ev->cond_, &ev->mutex_, &absoluteTimeout) != 0);
	ev->signaled_ = false;	//auto-reset signal
	pthread_mutex_unlock(&ev->mutex_);

	return !timeout;
}

inline void DestroyEvent(event_t *ev)
{
	pthread_mutex_destroy(&ev->mutex_);
	pthread_cond_destroy(&ev->cond_);
}

// MS Windows specific stuff here ...
//-------------------------------------------------------------------------------------
#elif defined(PLATFORM_WINDOWS)

#include <winsock2.h>
#include <ws2tcpip.h>
#include <process.h>
#include <stdexcept>
#pragma comment(lib, "ws2_32.lib")

typedef SOCKET socketd_t;
typedef HANDLE thread_t;
typedef HANDLE mutex_t;
typedef HANDLE event_t;

#define TCP_RECV_BUF_SIZE (64*1024)

#define MKR_WRITE( descriptor, buffer, size ) send( descriptor, buffer, (int) size, 0 )
#define MKR_READ( descriptor, buffer, size ) recv( descriptor, buffer, (int) size, 0 )
#define MKR_ACCEPT( descriptor, ptr_address, address_length ) accept( descriptor, NULL, NULL );
#define MKR_CLOSE_SOCKET( descriptor ) closesocket( descriptor )
#define MKR_SLEEP( t ) Sleep(t * 1000)

#define MKR_MUTEX_INIT( mutex ) { mutex = CreateMutex( NULL, false, NULL ); }
#define MKR_MUTEX_LOCK( mutex ) WaitForSingleObject( mutex, INFINITE )
#define MKR_MUTEX_UNLOCK( mutex ) ReleaseMutex( mutex )
#define MKR_MUTEX_DESTROY( mutex ) CloseHandle( mutex )

#define MKR_THREAD_CREATE( handle, start_routine, data ) CreateThread( handle, start_routine, data )
#define MKR_THREAD_EXIT() _endthreadex(0)
#define MKR_THREAD_JOIN( handle ) WaitForSingleObject( handle, INFINITE )
#define MKR_THREAD_DESTROY( handle ) CloseHandle( handle )

#define THREAD_CALLBACK_RETURN_TYPE unsigned
#define THREAD_CALLBACK_CALLING_CONVENTION __stdcall

inline bool CreateThread( thread_t *handle, THREAD_CALLBACK_RETURN_TYPE (THREAD_CALLBACK_CALLING_CONVENTION *start_routine)(void*), void *arg )
{
    if (handle == NULL)
    throw std::invalid_argument("Parameter 'handle' was NULL.");

    *handle = (HANDLE) _beginthreadex(NULL, 0, start_routine, arg, 0, NULL);

    return (*handle != 0);
}

inline void InitializeEvent(event_t *ev)
{
	*ev = CreateEvent(NULL, FALSE, FALSE, NULL);
}

inline void SetEvent(event_t *ev)
{
	SetEvent(*ev);
}

inline void ResetEvent(event_t *ev)
{
	ResetEvent(*ev);
}

inline bool WaitForEvent(event_t *ev, unsigned int timeoutMilliseconds)
{
	return (WaitForSingleObject(*ev, timeoutMilliseconds) == WAIT_OBJECT_0);
}

inline void DestroyEvent(event_t *ev)
{
	CloseHandle(*ev);
}

#endif // platform specific stuff
#endif // GDSNETWORKINTERFACE_PLATFORM_SPECIFIC
