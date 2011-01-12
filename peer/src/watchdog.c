/*
 * DELAYEDEXIT is a MATLAB mex file which can be used to exit Matlab
 * after a specified time. A thread is started in the background,
 * which starts a timer and exits when the timer elapses.
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
		watchdog.masterid = 0;
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

		pthread_cond_destroy(&condstatus);
		pthread_mutex_destroy(&mutexstatus);
		return;

} /* exitFun */

void mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]) {
		int rc;
		unsigned int time = 0, masterid = 0;

		/* this function will be called upon unloading of the mex file */
		mexAtExit(exitFun);

		if (nrhs<2)
				mexErrMsgTxt ("invalid number of input arguments");

		if (mxIsScalar(prhs[0]))
				masterid = mxGetScalar(prhs[0]);
		else if (mxIsEmpty(prhs[0]))
				masterid = 0;
		else
				mexErrMsgTxt ("invalid input argument #1");

		if (mxIsScalar(prhs[1]))
				time = mxGetScalar(prhs[1]);
		else if (mxIsEmpty(prhs[1]))
				time = 0;
		else
				mexErrMsgTxt ("invalid input argument #2");

		if (masterid!=0 || time!=0) {
				/* in this case the mex file is not allowed to be cleared from memory */
				if (!mexIsLocked())
						mexLock(); 
		}

		if (masterid==0 && time==0) {
				/* in this case the mex file is again allowed to be cleared from memory */
				if (mexIsLocked())
						mexUnlock(); 
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

		/* enable the watchdog: the expire thread will exit if the master is not seen any more */
		pthread_mutex_lock(&mutexwatchdog);
		watchdog.enabled  = 1;
		watchdog.masterid = masterid;
		watchdog.time     = time;
		mexPrintf("watchdog: enabled for masterid = %lu and time = %d\n", masterid, time);
		pthread_mutex_unlock(&mutexwatchdog);

		return;
} /* main */

