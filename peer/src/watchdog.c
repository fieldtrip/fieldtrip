/*
 * WATCHDOG is a MATLAB mex file which is used to exit Matlab
 * when the running peer process should be aborted.
 *
 * Copyright (C) 2010, Robert Oostenveld
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */

#include "mex.h"
#include "matrix.h"
#include <pthread.h>

#include "peer.h"
#include "extern.h"
#include "platform_includes.h"

int peerInitialized = 0;

/* the thread IDs are needed for thread cancelation at cleanup */
pthread_t discoverThread;
pthread_t expireThread;

/* this function will be called upon unloading of the mex file */
void exitFun(void) {

		if (!peerInitialized)
				return;

		/* disable the watchdog */
		pthread_mutex_lock(&mutexwatchdog);
		watchdog.enabled  = 0;
		watchdog.evidence = 0;
		watchdog.controllerid = 0;
		watchdog.memory   = 0;
		watchdog.time     = 0;
		mexPrintf("watchdog: disabled\n");
		pthread_mutex_unlock(&mutexwatchdog);

		/* stop the maintenance thread */
		pthread_mutex_lock(&mutexstatus);
		if (discoverStatus) {
				pthread_mutex_unlock(&mutexstatus);
				mexPrintf("watchdog: requesting cancelation of discover thread\n");
				pthread_cancel(discoverThread);
				pthread_join(discoverThread, NULL);
		}
		else {
				pthread_mutex_unlock(&mutexstatus);
		}

		/* stop the maintenance thread */
		pthread_mutex_lock(&mutexstatus);
		if (expireStatus) {
				pthread_mutex_unlock(&mutexstatus);
				mexPrintf("watchdog: requesting cancelation of expire thread\n");
				pthread_cancel(expireThread);
				pthread_join(expireThread, NULL);
		}
		else {
				pthread_mutex_unlock(&mutexstatus);
		}

		/* free the shared dynamical memory */
		peerexit(NULL);
		peerInitialized = 0;

		return;

} /* exitFun */

void mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]) {
		int rc, enabled;
		UINT32_T controllerid = 0;
		time_t timallow = 0;
		UINT64_T memallow = 0;
		UINT64_T rss, vs;

		/* this function will be called upon unloading of the mex file */
		mexAtExit(exitFun);

		if (nrhs<3)
				mexErrMsgTxt ("invalid number of input arguments");

		if (mxIsScalar(prhs[0]))
				controllerid = mxGetScalar(prhs[0]);
		else if (mxIsEmpty(prhs[0]))
				controllerid = 0;
		else
				mexErrMsgTxt ("invalid input argument #1");

		if (mxIsScalar(prhs[1]))
				timallow = mxGetScalar(prhs[1]);
		else if (mxIsEmpty(prhs[1]))
				timallow = 0;
		else
				mexErrMsgTxt ("invalid input argument #2");

		if (mxIsScalar(prhs[2]))
				memallow = mxGetScalar(prhs[2]);
		else if (mxIsEmpty(prhs[2]))
				memallow = 0;
		else
				mexErrMsgTxt ("invalid input argument #3");

		if (controllerid!=0 || timallow!=0 || memallow!=0) {
				enabled = 1;
				/* in this case the mex file is not allowed to be cleared from memory */
				if (!mexIsLocked())
						mexLock(); 
		}

		if (controllerid==0 && timallow==0 && memallow==0) {
				enabled = 0;
				/* in this case the mex file is allowed to be cleared from memory */
#ifdef SOLUTION_FOR_UNEXPLAINED_CRASH
				if (mexIsLocked())
						mexUnlock(); 
#endif
		}

		if (!peerInitialized) {
				mexPrintf("watchdog: init\n");
				peerinit(NULL);
				peerInitialized = 1;
		}

		/* start the discover thread */
		pthread_mutex_lock(&mutexstatus);
		if (!discoverStatus) {
				pthread_mutex_unlock(&mutexstatus);
				mexPrintf("watchdog: spawning discover thread\n");
				rc = pthread_create(&discoverThread, NULL, discover, (void *)NULL);
				if (rc)
						mexErrMsgTxt("problem with return code from pthread_create()");
				else {
						/* wait until the thread has properly started */
						pthread_mutex_lock(&mutexstatus);
						if (!discoverStatus)
								pthread_cond_wait(&condstatus, &mutexstatus);
						pthread_mutex_unlock(&mutexstatus);
				}
		}
		else {
				pthread_mutex_unlock(&mutexstatus);
		}

		/* start the expire thread */
		pthread_mutex_lock(&mutexstatus);
		if (!expireStatus) {
				pthread_mutex_unlock(&mutexstatus);
				mexPrintf("watchdog: spawning expire thread\n");
				rc = pthread_create(&expireThread, NULL, expire, (void *)NULL);
				if (rc)
						mexErrMsgTxt("problem with return code from pthread_create()");
				else {
						/* wait until the thread has properly started */
						pthread_mutex_lock(&mutexstatus);
						if (!expireStatus)
								pthread_cond_wait(&condstatus, &mutexstatus);
						pthread_mutex_unlock(&mutexstatus);
				}
		}
		else {
				pthread_mutex_unlock(&mutexstatus);
		}

		if (timallow>0) {
				/* timallow should be relative to now */
				timallow += time(NULL);
		}

		if (memallow>0) {
				/* memallow should be in absolute numbers, add the current memory footprint */
				getmem(&rss, &vs);
				memallow += rss;
		}

		/* enable the watchdog: the expire thread will exit if the controller is not seen any more */
		pthread_mutex_lock(&mutexwatchdog);
		watchdog.enabled  = enabled;
		watchdog.evidence = 0;
		watchdog.controllerid = controllerid;
		watchdog.memory   = memallow;
		watchdog.time     = timallow;
		pthread_mutex_unlock(&mutexwatchdog);

		if (enabled)
				mexPrintf("watchdog: enabled for controllerid = %lu, time = %d, memory = %lu\n", controllerid, timallow, memallow);
		else
				mexPrintf("watchdog: disabled for controllerid = %lu, time = %d, memory = %lu\n", controllerid, timallow, memallow);

		return;
} /* main */

