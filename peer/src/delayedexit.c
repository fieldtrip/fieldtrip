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

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "mex.h"
#include "matrix.h"
#include "platform.h"
#include "platform_includes.h"

#define PAUSE 1 /* in seconds */
#define FREE(x) {if (x) {free(x); x=NULL;}}

time_t timer = 0;
pthread_mutex_t mutextimer  = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutexstatus = PTHREAD_MUTEX_INITIALIZER;
pthread_t timerThread;
int timerStatus = 0;

/* this is the thread that checks the timer */
void *checktimer(void *arg) {

		pthread_mutex_lock(&mutexstatus);	
		if (timerStatus) {
				/* only a single instance should be running */
				pthread_mutex_unlock(&mutexstatus);	
				return;
		}
		else {
				timerStatus = 1;
				mexPrintf("starting delayedexit\n");
				pthread_mutex_unlock(&mutexstatus);	
		}

		while (1) {
				/* test whether the timer has elapsed */
				pthread_mutex_lock(&mutextimer);
				if (difftime(time(NULL), timer)>0) {
						/* do a brute force exit */
						exit(0);
				}
				pthread_mutex_unlock(&mutextimer);

				pthread_testcancel();
				sleep(PAUSE);
		}

		return NULL;
}

/* this function will be called upon unloading of the mex file */
void exitFun(void) {
		pthread_mutex_lock(&mutexstatus);	
		if (timerStatus) {
				timerStatus = 0;
				mexPrintf("stopping delayedexit\n");
				pthread_mutex_unlock(&mutexstatus);	
				pthread_cancel(timerThread);
				pthread_join(timerThread, NULL);
		}
		else {
				pthread_mutex_unlock(&mutexstatus);	
		}
}

void mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]) {
		int rc;
		double delay;

		/* this function will be called upon unloading of the mex file */
		mexAtExit(exitFun);

		if (nrhs<1) {
				/* show the status of the delayed exit */

				pthread_mutex_lock(&mutextimer);
				pthread_mutex_lock(&mutexstatus);
				if (timerStatus) {
						delay = difftime(timer-time(NULL));
						if (nlhs<1)
								mexPrintf("delayed exit scheduled at %d seconds\n", (int)delay);
						else
								plhs[0] = mxCreateDoubleScalar(delay);
				}
				else {
						if (nlhs<1)
								mexPrintf("no delayed exit scheduled\n");
						else
								plhs[0] = mxCreateDoubleScalar(mxGetInf());
				}
				pthread_mutex_unlock(&mutexstatus);
				pthread_mutex_unlock(&mutextimer);

		}

		else {
				/* the first argument is the delay in seconds */
				if (!mxIsScalar(prhs[0]))
						mexErrMsgTxt ("invalid input argument #1");

				delay = mxGetScalar(prhs[0]);

				/* set or update the timer */
				pthread_mutex_lock(&mutextimer);
				timer = time(NULL) + (unsigned int)delay;
				pthread_mutex_unlock(&mutextimer);

				/* start the thread */
				rc = pthread_create(&timerThread, NULL, checktimer, (void *)NULL);
				if (rc)
						mexErrMsgTxt("problem with return code from pthread_create()");
		}

} /* main */

