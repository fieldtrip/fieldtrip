/*
 * Copyright (C) 2012, Robert Oostenveld
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "mex.h"
#include "engine.h"
#include "platform.h"

/*
#define DEBUG_PRINT(string, ...) {fprintf(stderr, string, ## __VA_ARGS__);}
#define DEBUG_PRINT(string, ...) {mexPrintf(string, ## __VA_ARGS__);}
*/
#define DEBUG_PRINT(string, ...) {}
#define FREE(x)                  {if (x) {free(x); x=NULL;}}
#define TRUE            1
#define FALSE           0
#define MAX_THREADS     256
#define STRLEN          256

/************************************************************************/
#ifdef USE_PTHREADS
#include <pthread.h>

pthread_t                 threadid    [MAX_THREADS];
pthread_mutex_t           mutex_start [MAX_THREADS];
pthread_mutex_t           mutex_finish[MAX_THREADS];
pthread_cond_t            cond_start	[MAX_THREADS];
pthread_cond_t            cond_finish	[MAX_THREADS];
#define MUTEX_INIT(x,y)   pthread_mutex_init(x, y)
#define MUTEX_LOCK(x)     pthread_mutex_lock(x)
#define MUTEX_UNLOCK(x)   pthread_mutex_unlock(x)
#define MUTEX_DESTROY(x)  pthread_mutex_destroy(x)
#define COND_INIT(x,y)    pthread_cond_init(x, y)
#define COND_WAIT(x,y)    pthread_cond_wait(x, y)
#define COND_SIGNAL(x)    pthread_cond_signal(x)
#define COND_DESTROY(x)   pthread_cond_destroy(x)
#define THREAD_JOIN(x)    pthread_join(x, NULL)
#define THREAD_EXIT(x)    pthread_exit(x)

#else
#include <windows.h>

HANDLE                    threadid    [MAX_THREADS];
CRITICAL_SECTION          mutex_start [MAX_THREADS];
CRITICAL_SECTION          mutex_finish[MAX_THREADS];
CONDITION_VARIABLE        cond_start	[MAX_THREADS];
CONDITION_VARIABLE        cond_finish	[MAX_THREADS];
#define MUTEX_INIT(x,y)   InitializeCriticalSection(x)
#define MUTEX_LOCK(x)     EnterCriticalSection(x)
#define MUTEX_UNLOCK(x)   LeaveCriticalSection(x)
#define MUTEX_DESTROY(x)  CloseHandle(x)
#define COND_INIT(x,y)    InitializeConditionVariable(x)
#define COND_WAIT(x,y)    SleepConditionVariableCS(x, y, INFINITE)
#define COND_SIGNAL(x)    WakeConditionVariable(x)
#define COND_DESTROY(x)   CloseHandle(x)
#define THREAD_JOIN(x)    WaitForSingleObject(x, INFINITE)
#define THREAD_EXIT(x)

#define  sleep(x)         (Sleep((x)*1000)) /* in seconds      */
#define usleep(x)         (Sleep((x)/1000)) /* in microseconds */

#endif
/************************************************************************/

/* these are used in start and finish */
#define ENGINE_IDLE     1
#define ENGINE_BUSY     2
#define ENGINE_INVALID  3
/* these are used in start */
#define ENGINE_INFO     10
#define ENGINE_PUT      11
#define ENGINE_EVAL     12
#define ENGINE_GET      13
#define ENGINE_CLOSE    14

unsigned int    start [MAX_THREADS]; /* this is protexted by mutex_start */
unsigned int    finish[MAX_THREADS]; /* this is protexted by mutex_start */
unsigned int    poolsize      = 0;
unsigned int    initialized   = FALSE;
void           *input1        = NULL;
void           *input2        = NULL;
void           *output        = NULL;
char            matlabcmd[STRLEN];

/* this is called the first time that the mex-file is loaded */
void initFun(void)
{
		unsigned int engine;
		if (initialized==FALSE) {
				DEBUG_PRINT("Init\n");
				for (engine=0; engine<MAX_THREADS; engine++) {
						start [engine] = ENGINE_INVALID;
						finish[engine] = ENGINE_INVALID;
						MUTEX_INIT(&mutex_start[engine], NULL);
						MUTEX_INIT(&mutex_finish[engine], NULL);
						COND_INIT(&cond_start[engine], NULL);
						COND_INIT(&cond_finish[engine], NULL);
				}
				poolsize    = 0;
				initialized = TRUE;
		}
		return;
}

/* this function is called upon unloading of the mex-file */
void exitFun(void)
{
		unsigned int engine;
		if (initialized == TRUE) {
				DEBUG_PRINT("Exit\n");
				for (engine=0; engine<poolsize; engine++) {
						MUTEX_LOCK(&mutex_start[engine]);
						start[engine] = ENGINE_CLOSE;
						MUTEX_UNLOCK(&mutex_start[engine]);
						COND_SIGNAL(&cond_start[engine]);
						THREAD_JOIN(threadid[engine]);
				}
				for (engine=0; engine<MAX_THREADS; engine++) {
						start [engine] = ENGINE_INVALID;
						finish[engine] = ENGINE_INVALID;
						MUTEX_DESTROY(&mutex_start[engine]);
						MUTEX_DESTROY(&mutex_finish[engine]);
						COND_DESTROY(&cond_start[engine]);
						COND_DESTROY(&cond_finish[engine]);
				}
				initialized = FALSE;
				poolsize    = 0;
		}
		if (mexIsLocked())
				mexUnlock();
		return;
}

/* this function will be started as a separate thread */
void engineThread(void *argin)
{
		Engine         *ep = NULL;
		char            str[STRLEN];
		unsigned int    engine;
    int             operation, retval;

		if (!mexIsLocked()) {
				DEBUG_PRINT("Locking mex file\n");
				mexLock();
		}

		/* store the identifier for this thread, it is used for indexing the global arrays */
		engine = (int)argin;

		DEBUG_PRINT("Starting engine thread %d\n", engine);

		/* start the MATLAB engine, this thread will remain responsible for it */
#ifdef PLATFORM_WINDOWS
		ep = engOpenSingleUse(NULL, NULL, &retval);
    if (ep)
      engSetVisible(ep, 0);
#else
		ep = engOpen(matlabcmd);
#endif

		COND_SIGNAL(&cond_finish[engine]);

    if (!ep) {
      /* failed to start the engine */
      MUTEX_LOCK(&mutex_finish[engine]);
      finish[engine] = ENGINE_INVALID;
      MUTEX_UNLOCK(&mutex_finish[engine]);
      THREAD_EXIT(NULL);
      return;
    }
    
		while (1) {
				DEBUG_PRINT("Entering while loop in engine thread %d\n", engine);

				MUTEX_LOCK(&mutex_finish[engine]);
				finish[engine] = ENGINE_IDLE;
				MUTEX_UNLOCK(&mutex_finish[engine]);

				/* wait until the start condition is signalled */
				MUTEX_LOCK(&mutex_start[engine]);
				while (start[engine]==ENGINE_IDLE)
						COND_WAIT(&cond_start[engine], &mutex_start[engine]);
				/* make a local copy of the command */
				operation  = start[engine];
				start[engine] = ENGINE_IDLE;
				MUTEX_UNLOCK(&mutex_start[engine]);

				MUTEX_LOCK(&mutex_finish[engine]);
				finish[engine] = ENGINE_BUSY;
				MUTEX_UNLOCK(&mutex_finish[engine]);

				DEBUG_PRINT("Engine thread executing %d in %d\n", operation, engine);

				if (operation==ENGINE_CLOSE) {
						/* close the engine belonging to this thread and exit the thread */
						retval = engClose(ep);
						COND_SIGNAL(&cond_finish[engine]);
						THREAD_EXIT(NULL);
						return;
				}
				else if (operation==ENGINE_PUT) {
						/* copy a variable from the local matlab workspace to the remote engine */
						retval = engPutVariable(ep, mxArrayToString(input1), input2);
						DEBUG_PRINT("put variable %s, retval = %d\n", mxArrayToString(input1), retval);
						input1 = NULL;
						input2 = NULL;
						COND_SIGNAL(&cond_finish[engine]);
				}
				else if (operation==ENGINE_GET) {
						/* copy a variable from the remote engine workspace to the local matlab */
						output = engGetVariable(ep, mxArrayToString(input1));
						DEBUG_PRINT("get variable %s, retval = %d\n", mxArrayToString(input1), output);
						input1 = NULL;
						COND_SIGNAL(&cond_finish[engine]);
				}
				else if (operation==ENGINE_EVAL) {
						/* evaluate the string on the remote engine */
						memset(str, 0, STRLEN);
						strncpy(str, mxArrayToString(input1), STRLEN);
						input1 = NULL;
						COND_SIGNAL(&cond_finish[engine]);
						DEBUG_PRINT("Starting evaluation of '%s')\n", str);
						retval = engEvalString(ep, str);
						DEBUG_PRINT("Finished evaluation, retval = %d\n", retval);
				}
				else if (operation==ENGINE_INFO) {
						/* give some information about the engine */
						DEBUG_PRINT("Information for engine %d, retval = %d\n", engine, retval);
						COND_SIGNAL(&cond_finish[engine]);
				}
				else {
						COND_SIGNAL(&cond_finish[engine]);
				} /* if operation */

		} /* while */
}

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
		char            str[STRLEN];
		char           *ptr;
		int             retval;
    unsigned int    engine;

		initFun();
		mexAtExit(exitFun);

		/* the first argument is always the command string */
		if (nrhs < 1)
				mexErrMsgTxt("invalid number of input arguments");

		if (!mxIsChar(prhs[0]))
				mexErrMsgTxt("invalid input argument #1");
		if (mxGetString(prhs[0], str, STRLEN - 1))
				mexErrMsgTxt("invalid input argument #1");

		/* convert to lower case */
		ptr = str;
		while (*ptr) {
				*ptr = tolower(*ptr);
				ptr++;
		}

		DEBUG_PRINT("Starting MEX main loop for %s\n", str);

		/****************************************************************************/
		if (strcmp(str, "open") == 0) {
				/* engine open num matlabcmd */

				if (poolsize != 0)
						mexErrMsgTxt("There are already engines running");

				if (nrhs < 2)
						mexErrMsgTxt("Invalid number of input arguments");

				if (nrhs > 1) {
						if (!mxIsNumeric(prhs[1]))
								mexErrMsgTxt("Invalid input argument #2, should be numeric");
						poolsize = mxGetScalar(prhs[1]);
				}

				if (nrhs > 2) {
						if (!mxIsChar(prhs[2]))
								mexErrMsgTxt("Invalid input argument #3, should be a string");
						if (mxGetNumberOfElements(prhs[2])>STRLEN-1)
								mexErrMsgTxt("Invalid input argument #3, matlab command is too long");
						mxGetString(prhs[2], matlabcmd, STRLEN - 1);
				}
				else {
						sprintf(matlabcmd, "matlab -nosplash -nodisplay");
				}

				if (poolsize < 1)
						mexErrMsgTxt("The poolsize should be positive");

				if (poolsize > MAX_THREADS)
						mexErrMsgTxt("The poolsize is too large\n");

				for (engine=0; engine<poolsize; engine++) {
						MUTEX_LOCK(&mutex_finish[engine]);
#if USE_PTHREADS
						retval = pthread_create(&threadid[engine], NULL, (void *) engineThread, (void *) engine);
#else
						threadid[engine] = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE) engineThread, (void *) (engine), 0, NULL);
						retval = (threadid[engine]==NULL);
#endif
						while (finish[engine]!=ENGINE_IDLE)
								COND_WAIT(&cond_finish[engine], &mutex_finish[engine]);
						DEBUG_PRINT("The engine thread %d returned control\n", engine);
						MUTEX_UNLOCK(&mutex_finish[engine]);

						DEBUG_PRINT("Started engine %d, retval = %d\n", engine+1, retval);
						if (retval)
								break; /* a problem was detected in starting the threads, don't continue */
				};

				if (retval) {
						exitFun();	/* this cleans up all engines */
						mexErrMsgTxt("Failed to start engine threads");
				}
				else {
						DEBUG_PRINT("Succeeded in starting engine threads\n");
				}
		}
		/****************************************************************************/
		else if (strcmp(str, "put") == 0) {
				/* engine put number name value */

				if (nrhs != 4)
						mexErrMsgTxt("Exactly four input arguments needed");

				if (!mxIsNumeric(prhs[1]))
						mexErrMsgTxt("Argument #2 should be numeric");
				engine = mxGetScalar(prhs[1]);

				if (engine < 1 || engine > poolsize)
						mexErrMsgTxt("Invalid engine number");
				engine--; /* switch from MATLAB to C indexing */

				if (!mxIsChar(prhs[2]))
						mexErrMsgTxt("Argument #3 should be a string");

				MUTEX_LOCK(&mutex_finish[engine]);
				if (finish[engine]!=ENGINE_IDLE) {
						MUTEX_UNLOCK(&mutex_finish[engine]);
						/* wait for 1 second and try once more */
						usleep(1000000);
						MUTEX_LOCK(&mutex_finish[engine]);
				}
				if (finish[engine]!=ENGINE_IDLE) {
						MUTEX_UNLOCK(&mutex_finish[engine]);
						mexErrMsgTxt("The specified engine is not available");
				}

				input1 = (void *)prhs[2]; /* name */
				input2 = (void *)prhs[3]; /* value */

				MUTEX_LOCK(&mutex_start[engine]);
				start[engine] = ENGINE_PUT;
				MUTEX_UNLOCK(&mutex_start[engine]);

				DEBUG_PRINT("Signalling ENGINE_PUT in engine thread %d\n", engine);
				COND_SIGNAL(&cond_start[engine]);
				while (finish[engine]!=ENGINE_IDLE)
						COND_WAIT(&cond_finish[engine], &mutex_finish[engine]);
				while (input1!=NULL)
						COND_WAIT(&cond_finish[engine], &mutex_finish[engine]);
				while (input2!=NULL)
						COND_WAIT(&cond_finish[engine], &mutex_finish[engine]);
				DEBUG_PRINT("The engine thread %d returned control\n", engine);
				MUTEX_UNLOCK(&mutex_finish[engine]);

		}
		/****************************************************************************/
		else if (strcmp(str, "get") == 0) {
				/* engine get number name */

				if (nrhs != 3)
						mexErrMsgTxt("Exactly three input arguments needed");

				if (!mxIsNumeric(prhs[1]))
						mexErrMsgTxt("Argument #2 should be numeric");

				if (!mxIsChar(prhs[2]))
						mexErrMsgTxt("Argument #3 should be a string");
				engine = mxGetScalar(prhs[1]);

				if (engine < 1 || engine > poolsize)
						mexErrMsgTxt("Invalid engine number");
				engine--; /* switch from MATLAB to C indexing */

				MUTEX_LOCK(&mutex_finish[engine]);
				if (finish[engine]!=ENGINE_IDLE) {
						MUTEX_UNLOCK(&mutex_finish[engine]);
						/* wait for 1 second and try once more */
						usleep(1000000);
						MUTEX_LOCK(&mutex_finish[engine]);
				}
				if (finish[engine]!=ENGINE_IDLE) {
						MUTEX_UNLOCK(&mutex_finish[engine]);
						mexErrMsgTxt("The specified engine is not available");
				}

				input1 = (void *)prhs[2]; /* name */

				MUTEX_LOCK(&mutex_start[engine]);
				start[engine] = ENGINE_GET;
				MUTEX_UNLOCK(&mutex_start[engine]);

				DEBUG_PRINT("Signalling ENGINE_GET in engine thread %d\n", engine);
				COND_SIGNAL(&cond_start[engine]);
				while (finish[engine]!=ENGINE_IDLE)
						COND_WAIT(&cond_finish[engine], &mutex_finish[engine]);
				while (input1!=NULL)
						COND_WAIT(&cond_finish[engine], &mutex_finish[engine]);
				DEBUG_PRINT("The engine thread %d returned control\n", engine);
				plhs[0] = (mxArray *)output;
				output  = NULL;
				MUTEX_UNLOCK(&mutex_finish[engine]);
		}
		/****************************************************************************/
		else if (strcmp(str, "eval") == 0) {
				/* engine eval num str */

				if (nrhs < 3)
						mexErrMsgTxt("At least three input arguments needed");

				if (!mxIsNumeric(prhs[1]))
						mexErrMsgTxt("Argument #2 should be numeric");
				engine = mxGetScalar(prhs[1]);

				if (!mxIsChar(prhs[2]))
						mexErrMsgTxt("Invalid input argument #3, should be a string");
				if (mxGetNumberOfElements(prhs[2])>STRLEN-1)
						mexErrMsgTxt("Invalid input argument #3, command is too long");

				if (engine < 1 || engine > poolsize)
						mexErrMsgTxt("Invalid engine number");
				engine--; /* switch from MATLAB to C indexing */

				MUTEX_LOCK(&mutex_finish[engine]);
				if (finish[engine]!=ENGINE_IDLE) {
						MUTEX_UNLOCK(&mutex_finish[engine]);
						/* wait for 1 second and try once more */
						usleep(1000000);
						MUTEX_LOCK(&mutex_finish[engine]);
				}
				if (finish[engine]!=ENGINE_IDLE) {
						MUTEX_UNLOCK(&mutex_finish[engine]);
						mexErrMsgTxt("The specified engine is not available");
				}

				input1 = (void *)prhs[2]; /* command */

				MUTEX_LOCK(&mutex_start[engine]);
				start[engine] = ENGINE_EVAL;
				MUTEX_UNLOCK(&mutex_start[engine]);

				DEBUG_PRINT("Signalling ENGINE_EVAL in engine thread %d\n", engine);
				COND_SIGNAL(&cond_start[engine]);
				while (finish[engine]!=ENGINE_IDLE)
						COND_WAIT(&cond_finish[engine], &mutex_finish[engine]);
				while (input1!=NULL)
						COND_WAIT(&cond_finish[engine], &mutex_finish[engine]);
				DEBUG_PRINT("The engine thread %d returned control\n", engine);
				MUTEX_UNLOCK(&mutex_finish[engine]);
		}
		/****************************************************************************/
		else if (strcmp(str, "info") == 0) {
				/* engine info */

				for (engine=0; engine<poolsize; engine++) {

						MUTEX_LOCK(&mutex_finish[engine]);
						if (finish[engine]!=ENGINE_IDLE) {
								MUTEX_UNLOCK(&mutex_finish[engine]);
								/* wait for 1 second and try once more */
								usleep(1000000);
								MUTEX_LOCK(&mutex_finish[engine]);
						}
						if (finish[engine]!=ENGINE_IDLE) {
								MUTEX_UNLOCK(&mutex_finish[engine]);
								mexErrMsgTxt("The specified engine is not available");
						}

						MUTEX_LOCK(&mutex_start[engine]);
						start[engine] = ENGINE_INFO;
						MUTEX_UNLOCK(&mutex_start[engine]);

						DEBUG_PRINT("Signalling start for ENGINE_INFO\n");
						COND_SIGNAL(&cond_start[engine]);
						while (finish[engine]!=ENGINE_IDLE)
								COND_WAIT(&cond_finish[engine], &mutex_finish[engine]);
						DEBUG_PRINT("The engine thread %d returned control\n", engine);
						MUTEX_UNLOCK(&mutex_finish[engine]);
				}
		}
		/****************************************************************************/
		else if (strcmp(str, "isbusy") == 0) {
				/* engine isbusy num */

				if (nrhs != 2)
						mexErrMsgTxt("Exactly two input arguments needed");

				if (!mxIsNumeric(prhs[1]))
						mexErrMsgTxt("Argument #2 should be numeric");
				engine = mxGetScalar(prhs[1]);

				if (engine < 1 || engine > poolsize)
						mexErrMsgTxt("Invalid engine number");
				engine--; /* switch from MATLAB to C indexing */

				MUTEX_LOCK(&mutex_finish[engine]);
        retval = (finish[engine]==ENGINE_IDLE);
				MUTEX_UNLOCK(&mutex_finish[engine]);

				plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
				mxGetPr(plhs[0])[0] = (!retval);
		}
		/****************************************************************************/
		else if (strcmp(str, "poolsize") == 0) {
				/* engine poolsize */

				plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
				mxGetPr(plhs[0])[0] = poolsize;
		}
		/****************************************************************************/
		else if (strcmp(str, "close") == 0) {
				/* engine close */

				exitFun();	/* this cleans up all engines */
		}
		/****************************************************************************/
		else {
				mexErrMsgTxt("Unknown command");
				return;
		}

		DEBUG_PRINT("Finished MEX main loop\n");
		return;
}
