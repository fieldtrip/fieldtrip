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

		/* disable the kill switch */
		pthread_mutex_lock(&mutexkillswitch);
		killswitch.enabled  = 0;
		killswitch.masterid = 0;
		mexPrintf("killswitched: disabled\n");
		pthread_mutex_unlock(&mutexkillswitch);

		/* stop the maintenance thread */
		pthread_mutex_lock(&mutexstatus);
		if (discoverStatus) {
				pthread_mutex_unlock(&mutexstatus);
				mexPrintf("killswitch: requesting cancelation of discover thread\n");
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
				mexPrintf("killswitch: requesting cancelation of expire thread\n");
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
		int rc;

		/* this function will be called upon unloading of the mex file */
		mexAtExit(exitFun);

		if (nrhs<1)
				mexErrMsgTxt ("invalid number of input arguments");

		if (!mxIsScalar(prhs[0]))
				mexErrMsgTxt ("invalid input argument #1");

		if (!peerInitialized) {
				mexPrintf("killswitch: init\n");
				peerinit(NULL);
				peerInitialized = 1;
		}

		/* start the maintenance thread */
		pthread_mutex_lock(&mutexstatus);
		if (!discoverStatus) {
				pthread_mutex_unlock(&mutexstatus);
				mexPrintf("killswitch: spawning discover thread\n");
				rc = pthread_create(&discoverThread, NULL, discover, (void *)NULL);
				if (rc)
						mexErrMsgTxt("problem with return code from pthread_create()");
		}
		else {
				pthread_mutex_unlock(&mutexstatus);
		}

		/* start the maintenance thread */
		pthread_mutex_lock(&mutexstatus);
		if (!expireStatus) {
				pthread_mutex_unlock(&mutexstatus);
				mexPrintf("killswitch: spawning expire thread\n");
				rc = pthread_create(&expireThread, NULL, expire, (void *)NULL);
				if (rc)
						mexErrMsgTxt("problem with return code from pthread_create()");
		}
		else {
				pthread_mutex_unlock(&mutexstatus);
		}

		/* wait some time for the discovery thread to pick up the other peers */
		usleep(2*ANNOUNCESLEEP);

		/* enable the kill switch: the expire thread will exit if the master is not seen any more */
		pthread_mutex_lock(&mutexkillswitch);
		killswitch.enabled  = 1;
		killswitch.masterid = mxGetScalar(prhs[0]);
		mexPrintf("killswitch: enabled\n");
		pthread_mutex_unlock(&mutexkillswitch);

		return;
} /* main */

