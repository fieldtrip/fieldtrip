/*
 * Module: sched.h
 *
 * Purpose:
 *      Provides an implementation of POSIX realtime extensions
 *      as defined in 
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
#ifndef _SCHED_H
#define _SCHED_H

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

#include <windows.h>
#include <process.h>
#include <errno.h>

#ifdef _MSC_VER
/*
 * Re-enable all but 4127, 4514
 */
#pragma warning( default : 4115 4116 4201 4214)
#endif

#ifdef __cplusplus
extern "C"
{
#endif                          /* __cplusplus */

int sched_yield (void);

int sched_get_priority_min (int policy);

int sched_get_priority_max (int policy);


#ifdef __cplusplus
}                               /* End of extern "C" */
#endif                          /* __cplusplus */


#endif                          /* !_SCHED_H */
