/* This is the POSIX thread API (POSIX 1003).
 * 
 * Pthreads-win32 - POSIX Threads Library for Win32
 * Copyright (C) 1998
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public
 * License along with this library; if not, write to the Free
 * Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
 * MA 02111-1307, USA
 */

#if !defined( PTHREAD_H )
#define PTHREAD_H

/*
 * -------------------------------------------------------------
 *
 *
 * Module: pthread.h
 *
 * Purpose:
 *      Provides an implementation of PThreads based upon the
 *      standard:
 *
 *              POSIX 1003.1c-1995      (POSIX.1c)
 *
 * Authors:
 *      Contributors are listed in the file "MAINTAINERS".
 *
 * The following functions are implemented:
 *      ---------------------------
 *      PThreads
 *      ---------------------------
 *      pthread_attr_init
 *      pthread_attr_destroy
 *      pthread_attr_getdetachstate
 *      pthread_attr_getstackaddr
 *      pthread_attr_getstacksize
 *      pthread_attr_setdetachstate     
 *      pthread_attr_setstackaddr
 *      pthread_attr_setstacksize
 *
 *      pthread_create
 *      pthread_detach
 *      pthread_equal
 *      pthread_exit
 *      pthread_join
 *      pthread_self
 *      sched_yield
 *
 *      pthread_cancel
 *      pthread_cleanup_pop
 *      pthread_cleanup_push
 *      pthread_setcancelstate
 *      pthread_setcanceltype
 *      pthread_testcancel
 *
 *      ---------------------------
 *      Thread Specific Data
 *      ---------------------------
 *      pthread_key_create
 *      pthread_key_delete
 *      pthread_setspecific
 *      pthread_getspecific
 *
 *      ---------------------------
 *      Mutexes
 *      ---------------------------
 *      pthread_mutexattr_init
 *      pthread_mutexattr_destroy
 *      pthread_mutexattr_getpshared
 *      pthread_mutexattr_setpshared
 *
 *      pthread_mutex_init
 *      pthread_mutex_destroy
 *      pthread_mutex_lock
 *      pthread_mutex_trylock
 *      pthread_mutex_unlock
 *
 *      ---------------------------
 *      Condition Variables
 *      ---------------------------
 *      pthread_condattr_init
 *      pthread_condattr_destroy
 *      pthread_condattr_getpshared
 *      pthread_condattr_setpshared
 *
 *      pthread_cond_init
 *      pthread_cond_destroy
 *      pthread_cond_wait
 *      pthread_cond_timedwait
 *      pthread_cond_signal
 *      pthread_cond_broadcast
 *
 *      ---------------------------
 *      Protected Methods
 *      ---------------------------
 *      pthreadCancelableWait
 *
 *      ---------------------------
 *      RealTime Scheduling:
 *      ---------------------------
 *      pthread_attr_getschedparam
 *      pthread_attr_setschedparam
 *      pthread_getschedparam
 *      pthread_setschedparam
 *      sched_get_priority_max
 *      sched_get_priority_min
 *
 *      ---------------------------
 *      Signals:
 *      ---------------------------
 *      pthread_sigmask
 *
 * Limitations
 * ===========
 *      The following functions are not implemented:
 *
 *      ---------------------------
 *      RealTime Scheduling:
 *      ---------------------------
 *      pthread_attr_getinheritsched
 *      pthread_attr_getschedpolicy
 *      pthread_attr_getscope
 *      pthread_attr_setinheritsched
 *      pthread_attr_setschedpolicy
 *      pthread_attr_setscope
 *      pthread_mutex_getprioceiling
 *      pthread_mutex_setprioceiling
 *      pthread_mutex_attr_getprioceiling
 *      pthread_mutex_attr_getprotocol
 *      pthread_mutex_attr_setprioceiling
 *      pthread_mutex_attr_setprotocol
 *
 *      ---------------------------
 *      Fork Handlers:
 *      ---------------------------
 *      pthread_atfork
 *
 *      ---------------------------
 *      Stdio:
 *      ---------------------------
 *      flockfile
 *      ftrylockfile
 *      funlockfile
 *      getc_unlocked
 *      getchar_unlocked
 *      putc_unlocked
 *      putchar_unlocked
 *
 *      ---------------------------
 *      Thread-Safe C Runtime Library:
 *      ---------------------------
 *      readdir_r
 *      getgrgid_r
 *      getgrnam_r
 *      getpwuid_r
 *      getpwnam_r
 *
 *      ---------------------------
 *      Signals:
 *      ---------------------------
 *      pthread_kill
 *      sigtimedwait
 *      sigwait
 *      sigwaitinfo
 *
 *
 * -------------------------------------------------------------
 */

#ifdef _MSC_VER
/*
 * Disable following warnings when including Windows headers
 *
 * warning C4115: named type definition in parentheses
 * warning C4116: unnamed type definition in parentheses
 * warning C4127: conditional expression is constant
 * warning C4201: nonstandard extension used : nameless struct/union
 * warning C4214: nonstandard extension used : bit field types other than int
 * warning C4514: unreferenced inline function has been removed
 */
#pragma warning( disable : 4115 4116 4127 4201 4214 4514)
#endif

/*
 * -----------------
 * autoconf switches
 * -----------------
 */

#if HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */

#include <windows.h>
#include <time.h>

#if HAVE_SIGNAL_H
#include <signal.h>
#endif /* HAVE_SIGNAL_H */

#include <malloc.h>

#ifndef HAVE_STRUCT_TIMESPEC
struct timespec {
	int tv_sec;
	int tv_nsec;
};
#endif /* HAVE_STRUCT_TIMESPEC */

#ifndef SIG_BLOCK
#define SIG_BLOCK 0
#endif /* SIG_BLOCK */

#ifndef SIG_UNBLOCK 
#define SIG_UNBLOCK 1
#endif /* SIG_UNBLOCK */

#ifndef SIG_SETMASK
#define SIG_SETMASK 2
#endif /* SIG_SETMASK */


#include <process.h>
#include <errno.h>

#ifdef _WIN32
#ifndef ETIMEDOUT
#define ETIMEDOUT 19981220     /* Let's hope this is unique */
#endif
#endif /* _WIN32 */

#ifdef _MSC_VER
/*
 * Re-enable all but 4127, 4514
 */
#pragma warning( default : 4115 4116 4201 4214)
#endif

#if !defined( TRUE )
#define TRUE	!FALSE
#define FALSE	0
#endif /* !TRUE */

#ifdef __MINGW32__
#define PT_STDCALL
#else
#define PT_STDCALL __stdcall
#endif

/* 
 * This should perhaps be in autoconf or 
 * possibly fixed in Mingw32 to
 * correspond to the Windows headers.
 */
#ifdef __MINGW32__
#define _timeb timeb
#endif

#ifdef __cplusplus
extern "C"
{
#endif				/* __cplusplus */

/*
 * -------------------------------------------------------------
 *
 * POSIX 1003.1c-1995 Options
 * ===========================
 *
 * _POSIX_THREADS (set)
 *                      If set, you can use threads
 *
 * _POSIX_THREAD_ATTR_STACKSIZE (set)
 *                      If set, you can control the size of a thread's
 *                      stack
 *                              pthread_attr_getstacksize
 *                              pthread_attr_setstacksize
 *
 * _POSIX_THREAD_ATTR_STACKADDR (not set)
 *                      If set, you can allocate and control a thread's
 *                      stack. If not supported, the following functions
 *                      will return ENOSYS, indicating they are not
 *                      supported:
 *                              pthread_attr_getstackaddr
 *                              pthread_attr_setstackaddr
 *
 * _POSIX_THREAD_PRIORITY_SCHEDULING (not set)
 *                      If set, you can use realtime scheduling.
 *                      Indicates the availability of:
 *                              pthread_attr_getinheritsched
 *                              pthread_attr_getschedparam
 *                              pthread_attr_getschedpolicy
 *                              pthread_attr_getscope
 *                              pthread_attr_setinheritsched
 *                              pthread_attr_setschedparam
 *                              pthread_attr_setschedpolicy
 *                              pthread_attr_setscope
 *                              pthread_getschedparam
 *                              pthread_setschedparam
 *                              pthread_get_priority_max
 *                              pthread_get_priority_min
 *
 * _POSIX_THREAD_PRIO_INHERIT (not set)
 *                      If set, you can create priority inheritance
 *                      mutexes.
 *                              pthread_mutexattr_getprotocol +
 *                              pthread_mutexattr_setprotocol +
 *
 * _POSIX_THREAD_PRIO_PROTECT (not set)
 *                      If set, you can create priority ceiling mutexes
 *                      Indicates the availability of:
 *                              pthread_mutex_getprioceiling
 *                              pthread_mutex_setprioceiling
 *                              pthread_mutexattr_getprioceiling
 *                              pthread_mutexattr_getprotocol     +
 *                              pthread_mutexattr_setprioceiling
 *                              pthread_mutexattr_setprotocol     +
 *
 * _POSIX_THREAD_PROCESS_SHARED (not set)
 *                      If set, you can create mutexes and condition
 *                      variables that can be shared with another
 *                      process.If set, indicates the availability
 *                      of:
 *                              pthread_mutexattr_getpshared
 *                              pthread_mutexattr_setpshared
 *                              pthread_condattr_getpshared
 *                              pthread_condattr_setpshared
 *
 * _POSIX_THREAD_SAFE_FUNCTIONS (set)
 *                      If set you can use the special *_r library
 *                      functions that provide thread-safe behaviour
 *
 *      + These functions provide both 'inherit' and/or
 *        'protect' protocol, based upon these macro
 *        settings.
 *
 * POSIX 1003.1c-1995 Limits
 * ===========================
 *
 * PTHREAD_DESTRUCTOR_ITERATIONS
 *                      Maximum number of attempts to destroy
 *                      a thread's thread-specific data on
 *                      termination (must be at least 4)
 *
 * PTHREAD_KEYS_MAX
 *                      Maximum number of thread-specific data keys
 *                      available per process (must be at least 128)
 *
 * PTHREAD_STACK_MIN
 *                      Minimum supported stack size for a thread
 *
 * PTHREAD_THREADS_MAX
 *                      Maximum number of threads supported per
 *                      process (must be at least 64).
 *
 * -------------------------------------------------------------
 */

/*
 * POSIX Options
 */
#define _POSIX_THREADS
#define _POSIX_THREAD_SAFE_FUNCTIONS

#define _POSIX_THREAD_ATTR_STACKSIZE

#if defined( KLUDGE )
/*
 * The following are not supported
 */
#define _POSIX_THREAD_ATTR_STACKADDR
#define _POSIX_THREAD_PRIORITY_SCHEDULING
#define _POSIX_THREAD_PRIO_INHERIT
#define _POSIX_THREAD_PRIO_PROTECT
#define _POSIX_THREAD_PROCESS_SHARED

#endif				/* KLUDGE */

/*
 * POSIX Limits
 *
 *      PTHREAD_DESTRUCTOR_ITERATIONS
 *              Standard states this must be at least
 *              4.
 *
 *      PTHREAD_KEYS_MAX
 *              WIN32 permits only 64 TLS keys per process.
 *              This limitation could be worked around by
 *              simply simulating keys.
 *
 *      PTHREADS_STACK_MIN
 *              artibrarily chose 1K. By default, WIN32
 *              selects 1Meg stacks.
 *
 *      PTHREAD_THREADS_MAX
 *              Not documented by WIN32. Wrote a test program
 *              that kept creating threads until it failed
 *              revealed this approximate number.
 *
 */
#define PTHREAD_DESTRUCTOR_ITERATIONS	                   4
#define PTHREAD_KEYS_MAX				  64
#define PTHREAD_STACK_MIN				1024
#define PTHREAD_THREADS_MAX				2019


typedef struct pthread_t_ *pthread_t;
typedef struct pthread_attr_t_ *pthread_attr_t;
typedef struct pthread_once_t_ pthread_once_t;
typedef struct pthread_key_t_ *pthread_key_t;
typedef struct pthread_mutex_t_ *pthread_mutex_t;
typedef struct pthread_mutexattr_t_ *pthread_mutexattr_t;
typedef struct pthread_cond_t_ *pthread_cond_t;
typedef struct pthread_condattr_t_ *pthread_condattr_t;


/*
 * ====================
 * ====================
 * POSIX Threads
 * ====================
 * ====================
 */

/*
 * pthread_attr_{get,set}detachstate
 */
#define PTHREAD_CREATE_JOINABLE		0
#define PTHREAD_CREATE_DETACHED		1

/*
 * pthread_attr{get,set}inheritsched
 */
#define PTHREAD_INHERIT_SCHED		0
#define PTHREAD_EXPLICIT_SCHED		1

/*
 * pthread_setcancelstate paramters
 */
#define PTHREAD_CANCEL_ENABLE		0
#define PTHREAD_CANCEL_DISABLE		1

/*
 * pthread_setcanceltype parameters
 */
#define PTHREAD_CANCEL_ASYNCHRONOUS	0
#define PTHREAD_CANCEL_DEFERRED		1

/*
 * pthread_mutexattr_{get,set}pshared
 * pthread_condattr_{get,set}pshared
 */
#define PTHREAD_PROCESS_PRIVATE		0
#define PTHREAD_PROCESS_SHARED		1

/*
 * pthread_mutexattr_setforcecs_np
 */
#define PTHREAD_MUTEX_AUTO_CS_NP        0
#define PTHREAD_MUTEX_FORCE_CS_NP       1

/*
 * ====================
 * ====================
 * Cancelation
 * ====================
 * ====================
 */
#define PTHREAD_CANCELED       ((void *) -1)


/*
 * ====================
 * ====================
 * Once Key
 * ====================
 * ====================
 */
#define PTHREAD_ONCE_INIT       { FALSE, -1 }

struct pthread_once_t_
{
  int done;                 /* indicates if user function executed  */
  long started;             /* First thread to increment this value */
                            /* to zero executes the user function   */
};


#define PTHREAD_MUTEX_INITIALIZER ((pthread_mutex_t) -1)

#define PTHREAD_COND_INITIALIZER ((pthread_cond_t) -1)


/*
 * ====================
 * ====================
 * Scheduling
 * ====================
 * ====================
 */

/* Thread scheduling policies */

#define SCHED_OTHER 0
#define SCHED_FIFO  1
#define SCHED_RR    2

#define SCHED_MIN   SCHED_OTHER
#define SCHED_MAX   SCHED_RR

struct sched_param {
  int sched_priority;
};


/* There are three implementations of cancel cleanup.
 * Note that pthread.h is included in both application
 * compilation units and also internally for the library.
 * The code here and within the library aims to work
 * for all reasonable combinations of environments.
 * For example, although the library itself can't be
 * built (yet) in C, an application written in C can
 * be linked and run against a library built using
 * either WIN32 SEH or C++ EH.
 *
 * The three implementations are:
 *
 *   WIN32 SEH
 *   C
 *   C++
 */

typedef struct _pthread_cleanup_t _pthread_cleanup_t;

struct _pthread_cleanup_t
{
  void (*routine) (void *);
  void *arg;
#if !defined(_MSC_VER) && !defined(__cplusplus)
  _pthread_cleanup_t *prev;
#endif
};

#ifdef _MSC_VER
	/*
	 * WIN32 SEH version of cancel cleanup.
	 */

#define pthread_cleanup_push( _rout, _arg ) \
	{ \
	    _pthread_cleanup_t	_cleanup; \
	    \
            _cleanup.routine	= (_rout); \
	    _cleanup.arg	= (_arg); \
	    __try \
	      { \

#define pthread_cleanup_pop( _execute ) \
	      } \
	    __finally \
		{ \
		    if( _execute || AbnormalTermination()) \
		      { \
			  (*(_cleanup.routine))( _cleanup.arg ); \
		      } \
		} \
	}

#else /* _MSC_VER */

#ifndef __cplusplus

	/*
	 * C implementation of PThreads cancel cleanup
	 */

#define pthread_cleanup_push( _rout, _arg ) \
	{ \
	    _pthread_cleanup_t	_cleanup; \
            \
	    pthread_push_cleanup( &_cleanup, (_rout), (_arg) ); \

#define pthread_cleanup_pop( _execute ) \
	    (void) pthread_pop_cleanup( _execute ); \
	}

#else /* !__cplusplus */

	/*
	 * C++ version of cancel cleanup.
	 * - John E. Bossom.
	 */

	class PThreadCleanup {
	  /*
	   * PThreadCleanup
	   *
	   * Purpose
	   *      This class is a C++ helper class that is
	   *      used to implement pthread_cleanup_push/
	   *      pthread_cleanup_pop.
	   *      The destructor of this class automatically
	   *      pops the pushed cleanup routine regardless
	   *      of how the code exits the scope
	   *      (i.e. such as by an exception)
	   */
	  void            (*cleanUpRout)( void * );
	  void    *       obj;
	  int             executeIt;

	public:
	  PThreadCleanup() :
	    cleanUpRout( NULL ),
	    obj( NULL ),
	    executeIt( 0 )
	    /*
	     * No cleanup performed
	     */
	    {
	    }

	  PThreadCleanup(
			 void            (*routine)( void * ),
			 void    *       arg ) :
	    cleanUpRout( routine ),
	    obj( arg ),
	    executeIt( 1 )
	    /*
             * Registers a cleanup routine for 'arg'
             */
	    {
	    }

	  ~PThreadCleanup()
	    {
	      if ( executeIt && cleanUpRout != NULL )
		{
		  (*cleanUpRout)( obj );
		}
	    }

	  void execute( int exec )
	    {
	      executeIt = exec;
	    }
	};

	/*
	 * C++ implementation of PThreads cancel cleanup;
	 * This implementation takes advantage of a helper
	 * class who's destructor automatically calls the
	 * cleanup routine if we exit our scope weirdly
	 */
#define pthread_cleanup_push( _rout, _arg ) \
        { \
	    PThreadCleanup  cleanup((void (PT_STDCALL *)(void *))(_rout), \
				    (void *) (_arg) );

#define pthread_cleanup_pop( _execute ) \
	    cleanup.execute( _execute ); \
	}

#endif /* !__cplusplus */

#endif /* _MSC_VER */

/*
 * ===============
 * ===============
 * Methods
 * ===============
 * ===============
 */

/*
 * PThread Attribute Functions
 */
int pthread_attr_init (pthread_attr_t * attr);

int pthread_attr_destroy (pthread_attr_t * attr);

int pthread_attr_getdetachstate (const pthread_attr_t * attr,
				 int *detachstate);

int pthread_attr_getstackaddr (const pthread_attr_t * attr,
			       void **stackaddr);

int pthread_attr_getstacksize (const pthread_attr_t * attr,
			       size_t * stacksize);

int pthread_attr_setdetachstate (pthread_attr_t * attr,
				 int detachstate);

int pthread_attr_setstackaddr (pthread_attr_t * attr,
			       void *stackaddr);

int pthread_attr_setstacksize (pthread_attr_t * attr,
			       size_t stacksize);

/*
 * PThread Functions
 */
int pthread_create (pthread_t * tid,
		    const pthread_attr_t * attr,
		    void *(*start) (void *),
		    void *arg);

int pthread_detach (pthread_t tid);

int pthread_equal (pthread_t t1,
		   pthread_t t2);

void pthread_exit (void *value_ptr);

int pthread_join (pthread_t thread,
		  void **value_ptr);

pthread_t pthread_self (void);

int pthread_cancel (pthread_t thread);

int pthread_setcancelstate (int state,
			    int *oldstate);

int pthread_setcanceltype (int type,
			   int *oldtype);

void pthread_testcancel (void);

int pthread_once (pthread_once_t * once_control,
		  void (*init_routine) (void));

_pthread_cleanup_t *pthread_pop_cleanup (int execute);

void pthread_push_cleanup (_pthread_cleanup_t * cleanup,
			   void (*routine) (void *),
			   void *arg);

/*
 * Thread Specific Data Functions
 */
int pthread_key_create (pthread_key_t * key,
			void (*destructor) (void *));

int pthread_key_delete (pthread_key_t key);

int pthread_setspecific (pthread_key_t key,
			 const void *value);

void *pthread_getspecific (pthread_key_t key);


/*
 * Mutex Attribute Functions
 */
int pthread_mutexattr_init (pthread_mutexattr_t * attr);

int pthread_mutexattr_destroy (pthread_mutexattr_t * attr);

int pthread_mutexattr_getpshared (const pthread_mutexattr_t
				  * attr,
				  int *pshared);

int pthread_mutexattr_setpshared (pthread_mutexattr_t * attr,
				  int pshared);

int pthread_mutexattr_setforcecs_np(pthread_mutexattr_t *attr,
				    int forcecs);

/*
 * Mutex Functions
 */
int pthread_mutex_init (pthread_mutex_t * mutex,
			const pthread_mutexattr_t * attr);

int pthread_mutex_destroy (pthread_mutex_t * mutex);

int pthread_mutex_lock (pthread_mutex_t * mutex);

int pthread_mutex_trylock (pthread_mutex_t * mutex);

int pthread_mutex_unlock (pthread_mutex_t * mutex);

/*
 * Condition Variable Attribute Functions
 */
int pthread_condattr_init (pthread_condattr_t * attr);

int pthread_condattr_destroy (pthread_condattr_t * attr);

int pthread_condattr_getpshared (const pthread_condattr_t * attr,
				 int *pshared);

int pthread_condattr_setpshared (pthread_condattr_t * attr,
				 int pshared);

/*
 * Condition Variable Functions
 */
int pthread_cond_init (pthread_cond_t * cond,
		       const pthread_condattr_t * attr);

int pthread_cond_destroy (pthread_cond_t * cond);

int pthread_cond_wait (pthread_cond_t * cond,
		       pthread_mutex_t * mutex);

int pthread_cond_timedwait (pthread_cond_t * cond,
			    pthread_mutex_t * mutex,
			    const struct timespec *abstime);

int pthread_cond_signal (pthread_cond_t * cond);

int pthread_cond_broadcast (pthread_cond_t * cond);

/*
 * Scheduling
 */

int pthread_setschedparam (pthread_t thread,
			   int policy,
			   const struct sched_param *param);

int pthread_getschedparam (pthread_t thread,
			   int *policy,
			   struct sched_param *param);

int pthread_attr_getschedparam (const pthread_attr_t *attr,
				struct sched_param *param);

int pthread_attr_setschedparam (pthread_attr_t *attr,
				const struct sched_param *param);

/*
 * Protected Methods
 *
 * This function blocks until the given WIN32 handle
 * is signaled or pthread_cancel had been called.
 * This function allows the caller to hook into the
 * PThreads cancel mechanism. It is implemented using
 *
 *              WaitForMultipleObjects
 *
 * on 'waitHandle' and a manually reset WIN32 Event
 * used to implement pthread_cancel. The 'timeout'
 * argument to TimedWait is simply passed to
 * WaitForMultipleObjects.
 */
int pthreadCancelableWait (HANDLE waitHandle);
int pthreadCancelableTimedWait (HANDLE waitHandle, DWORD timeout);

/*
 * Thread-Safe C Runtime Library Mappings.
 */
#if ! defined( _REENTRANT ) && ! defined( _MT )
int * _errno( void );
#endif

/*
 * WIN32 C runtime library had been made thread-safe
 * without affecting the user interface. Provide
 * mappings from the UNIX thread-safe versions to
 * the standard C runtime library calls.
 * Only provide function mappings for functions that
 * actually exist on WIN32.
 */

#if !defined(__MINGW32__)
#define strtok_r( _s, _sep, _lasts ) \
	( *(_lasts) = strtok( (_s), (_sep) ) )
#endif /* !__MINGW32__ */

#define asctime_r( _tm, _buf ) \
	( strcpy( (_buf), asctime( (_tm) ) ), \
	  (_buf) )

#define ctime_r( _clock, _buf ) \
	( strcpy( (_buf), ctime( (_tm) ) ), \
	  (_buf) )

#define gmtime_r( _clock, _result ) \
	( *(_result) = *gmtime( (_clock) ), \
	  (_result) )

#define localtime_r( _clock, _result ) \
	( *(_result) = *localtime( (_clock) ), \
	  (_result) )

#define rand_r( _seed ) \
	rand()


#ifdef __cplusplus
}				/* End of extern "C" */
#endif				/* __cplusplus */

#endif /* PTHREAD_H */
