/*
 * Module: semaphore.h
 *
 * Purpose:
 *      Semaphores aren't actually part of the PThreads standard.
 *      They are defined by the POSIX Standard:
 *
 *              POSIX 1003.1b-1993      (POSIX.1b)
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
#if !defined( SEMAPHORE_H )
#define SEMAPHORE_H

#include <process.h>
#include <errno.h>

#define _POSIX_SEMAPHORES

#ifdef __cplusplus
extern "C"
{
#endif                          /* __cplusplus */

#if defined(_MSC_VER)
typedef unsigned int mode_t;
#endif

typedef HANDLE sem_t;

int sem_init (sem_t * sem,
	      int pshared,
	      unsigned int value
	      );

int sem_destroy (sem_t * sem
		 );

int sem_trywait (sem_t * sem
		 );

int sem_wait (sem_t * sem
	      );

int sem_post (sem_t * sem
	      );

int sem_open (const char * name,
	      int oflag,
	      mode_t mode,
            unsigned int value
	      );

int sem_close (sem_t * sem
	       );

int sem_unlink (const char * name
		);

int sem_getvalue (sem_t * sem,
		  int * sval
		  );

#ifdef __cplusplus
}                               /* End of extern "C" */
#endif                          /* __cplusplus */

#endif                          /* !SEMAPHORE_H */
