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

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "mex.h"
#include "engine.h"
#include "platform.h"

typedef struct {
		Engine         *ep;
		pthread_t       tid;
		int             busy;
		char           *cmd;
		int             retval;
} engine_t;

#define FREE(x)        {if (x) {free(x); x=NULL;}}
#define TRUE            1
#define FALSE           0
#define MAX_THREADS     256
#define STRLEN          256

/* these are used in status */
#define ENGINE_UNKNOWN  0
#define ENGINE_IDLE     1
#define ENGINE_BUSY     2

/* these are used in command */
#define ENGINE_INFO     10
#define ENGINE_PUT      11
#define ENGINE_EVAL     12
#define ENGINE_GET      13
#define ENGINE_CLOSE    14

pthread_mutex_t mutex_start [MAX_THREADS];
pthread_mutex_t mutex_finish[MAX_THREADS];
pthread_cond_t  cond_start	[MAX_THREADS];
pthread_cond_t  cond_finish	[MAX_THREADS];
pthread_t       threadid    [MAX_THREADS];
unsigned int    command     [MAX_THREADS];
unsigned int    status      [MAX_THREADS];
unsigned int    poolsize      = 0;
unsigned int    initialized   = FALSE;
unsigned int    retval        = 0;
void           *ptr1          = NULL;
void           *ptr2          = NULL;
char            matlabcmd[STRLEN];

/* this is called the first time that the mex-file is loaded */
void initFun(void)
{
		int engine;
		if (initialized==FALSE) {
				fprintf(stderr, "engine init()\n");
				for (engine=0; engine<MAX_THREADS; engine++) {
						status[engine] = ENGINE_UNKNOWN;
						pthread_mutex_init(&mutex_start[engine], NULL);
						pthread_mutex_init(&mutex_finish[engine], NULL);
						pthread_cond_init(&cond_start[engine], NULL);
						pthread_cond_init(&cond_finish[engine], NULL);
				}
				poolsize    = 0;
				initialized = TRUE;
		}
		return;
}

/* this function is called upon unloading of the mex-file */
void exitFun(void)
{
		int engine;
		fprintf(stderr, "engine exit()\n");
		for (engine=0; engine<poolsize; engine++) {
				pthread_mutex_lock(&mutex_start[engine]);
				command[engine] = ENGINE_CLOSE;
				pthread_mutex_unlock(&mutex_start[engine]);
				pthread_cond_signal(&cond_start[engine]);
				pthread_join(threadid[engine], NULL);
		}
		for (engine=0; engine<MAX_THREADS; engine++) {
				status[engine] = ENGINE_UNKNOWN;
				pthread_mutex_destroy(&mutex_start[engine]);
				pthread_mutex_destroy(&mutex_finish[engine]);
				pthread_cond_destroy(&cond_start[engine]);
				pthread_cond_destroy(&cond_finish[engine]);
		}
		initialized = FALSE;
		poolsize = 0;
		if (mexIsLocked())
				mexUnlock();
		return;
}

/* this function will be started as a seperate thread */
void engineThread(void *argin)
{
		Engine         *ep = NULL;
		char            cmd[STRLEN];
		int             thisId, retval;

		if (!mexIsLocked()) {
				fprintf(stderr, "Locking mex file\n");
				mexLock();
		}

		fprintf(stderr, "Starting thread\n");

		/* store the identifier for this thread, it is used for indexing the cmd array */
		thisId = (int)argin;

		/* start the MATLAB engine, this thread will remain responsible for it */
#ifdef PLATFORM_WINDOWS
		ep = engOpenSingleUse(NULL, NULL, &retval);
#else
		ep = engOpen(matlabcmd);
		retval = (ep!=NULL);
#endif

		status[thisId] = ENGINE_IDLE;
		pthread_cond_signal(&cond_finish[thisId]);

		while (1) {
				fprintf(stderr, "Entering while loop\n");
				pthread_mutex_lock(&mutex_start[thisId]);
				/* wait until the start condition is signalled */
				pthread_cond_wait(&cond_start[thisId], &mutex_start[thisId]);

				fprintf(stderr, "Engine thread executing %d in %d\n", command[thisId], thisId);
				status[thisId] = ENGINE_BUSY;
				retval = 0;

				if (command[thisId]==ENGINE_CLOSE) {
						/* close the engine belonging to this thread and exit the thread */
						retval = engClose(ep);
						status[thisId] = ENGINE_UNKNOWN;
						pthread_mutex_unlock(&mutex_start[thisId]); /* return to the main mex loop */
						pthread_cond_signal(&cond_finish[thisId]);
						pthread_exit(NULL);
						return;
				}
				else if (command[thisId]==ENGINE_PUT) {
						/* copy a variable from the local matlab workspace to the remote engine */
						retval = engPutVariable(ep, mxArrayToString(ptr1), ptr2);
						pthread_mutex_unlock(&mutex_start[thisId]); /* return to the main loop */
						pthread_cond_signal(&cond_finish[thisId]);
				}
				else if (command[thisId]==ENGINE_GET) {
						/* copy a variable from the remote engine workspace to the local matlab */
						ptr2 = engGetVariable(ep, mxArrayToString(ptr1));
						retval = (ptr2==NULL);
						pthread_mutex_unlock(&mutex_start[thisId]); /* return to the main loop */
						pthread_cond_signal(&cond_finish[thisId]);
				}
				else if (command[thisId]==ENGINE_EVAL) {
						/* evaluate the string on the remote engine */
						bzero(cmd, STRLEN);
						strncpy(cmd, mxArrayToString(ptr1), STRLEN);
						pthread_mutex_unlock(&mutex_start[thisId]); /* return to the main loop */
						pthread_cond_signal(&cond_finish[thisId]);
						fprintf(stderr, "Starting evaluation of '%s'\n", cmd);
						engEvalString(ep, cmd);
						fprintf(stderr, "Finished evaluation\n");
				}
				else if (command[thisId]==ENGINE_INFO) {
						/* give some information about the engine */
						fprintf(stderr, "Information for engine %d\n", thisId+1);
						pthread_mutex_unlock(&mutex_start[thisId]); /* return to the main loop */
						pthread_cond_signal(&cond_finish[thisId]);
				} /* if operation */

				status[thisId] = ENGINE_IDLE;

		} /* while */
}

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
		char            cmd[STRLEN];
		char           *ptr;
		int             retval, block, engine;

		initFun();
		mexAtExit(exitFun);

		/* the first argument is always the cmd string */
		if (nrhs < 1)
				mexErrMsgTxt("invalid number of input arguments");

		if (!mxIsChar(prhs[0]))
				mexErrMsgTxt("invalid input argument #1");
		if (mxGetString(prhs[0], cmd, STRLEN - 1))
				mexErrMsgTxt("invalid input argument #1");

		/* convert to lower case */
		ptr = cmd;
		while (*ptr) {
				*ptr = tolower(*ptr);
				ptr++;
		}

		/****************************************************************************/
		if (strcmp(cmd, "open") == 0) {
				/* engine open num cmd */

				if (poolsize > 0)
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
						mxGetString(prhs[2], matlabcmd, STRLEN - 1);
				}
				else {
						sprintf(matlabcmd, "matlab");
				}

				if (poolsize < 1)
						mexErrMsgTxt("The poolsize should be positive");

				if (poolsize > MAX_THREADS)
						mexErrMsgTxt("The poolsize is too large\n");

				for (engine=0; engine<poolsize; engine++) {
						pthread_mutex_lock(&mutex_finish[engine]);
						retval = pthread_create(&threadid[engine], NULL, (void *) engineThread, (void *) engine);
						pthread_cond_wait(&cond_finish[engine], &mutex_finish[engine]);
						pthread_mutex_unlock(&mutex_finish[engine]);

						fprintf(stderr, "Started engine %d, retval = %d\n", engine+1, retval);
						if (retval)
								break;
				};

				if (retval!=0) {
						exitFun();	/* this cleans up all engines */
						mexErrMsgTxt("Failed to start engine threads");
				}
				else {
						mexPrintf("Succeeded in starting engine threads\n");
				}
		}
		/****************************************************************************/
		else if (strcmp(cmd, "close") == 0) {
				/* engine close */

				exitFun();	/* this cleans up all engines */
				plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
				mxGetPr(plhs[0])[0] = 0;
		}
		/****************************************************************************/
		else if (strcmp(cmd, "put") == 0) {
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

				pthread_mutex_lock(&mutex_start[engine]);

				if (status[engine]!=ENGINE_IDLE) {
						pthread_mutex_unlock(&mutex_start[engine]);
						mexErrMsgTxt("The specified engine is not initialized or busy");
				}

				ptr1 = (void *)prhs[2]; /* name */
				ptr2 = (void *)prhs[3]; /* value */
				command[engine] = ENGINE_PUT;
				pthread_mutex_unlock(&mutex_start[engine]);

				pthread_mutex_lock(&mutex_finish[engine]);
				fprintf(stderr, "signalling thread for ENGINE_PUT\n");
				pthread_cond_signal(&cond_start[engine]);
				pthread_cond_wait(&cond_finish[engine], &mutex_finish[engine]);
				retval = 0;
				plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
				mxGetPr(plhs[0])[0] = retval;
				pthread_mutex_unlock(&mutex_finish[engine]);

		}
		/****************************************************************************/
		else if (strcmp(cmd, "get") == 0) {
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

				pthread_mutex_lock(&mutex_start[engine]);

				if (status[engine]!=ENGINE_IDLE) {
						pthread_mutex_unlock(&mutex_start[engine]);
						mexErrMsgTxt("The specified engine is not initialized or busy");
				}

				ptr1 = (void *)prhs[2]; /* name */
				command[engine] = ENGINE_GET;
				pthread_mutex_unlock(&mutex_start[engine]);

				pthread_mutex_lock(&mutex_finish[engine]);
				fprintf(stderr, "signalling thread for ENGINE_GET\n");
				pthread_cond_signal(&cond_start[engine]);
				pthread_cond_wait(&cond_finish[engine], &mutex_finish[engine]);
				plhs[0] = (mxArray *)ptr2;
				pthread_mutex_unlock(&mutex_finish[engine]);

		}
		/****************************************************************************/
		else if (strcmp(cmd, "eval") == 0) {
				/* engine eval num str */

				if (nrhs < 3)
						mexErrMsgTxt("At least three input arguments needed");

				if (!mxIsNumeric(prhs[1]))
						mexErrMsgTxt("Argument #2 should be numeric");
				engine = mxGetScalar(prhs[1]);

				if (!mxIsChar(prhs[2]))
						mexErrMsgTxt("Invalid input argument #3, should be a string");

				if (engine < 1 || engine > poolsize)
						mexErrMsgTxt("Invalid engine number");
				engine--; /* switch from MATLAB to C indexing */

				pthread_mutex_lock(&mutex_start[engine]);

				if (status[engine]!=ENGINE_IDLE) {
						pthread_mutex_unlock(&mutex_start[engine]);
						mexErrMsgTxt("The specified engine is not initialized or busy");
				}

				ptr1 = (void *)prhs[2]; /* command */
				command[engine] = ENGINE_EVAL;
				pthread_mutex_unlock(&mutex_start[engine]);

				pthread_mutex_lock(&mutex_finish[engine]);
				fprintf(stderr, "signalling thread for ENGINE_EVAL\n");
				pthread_cond_signal(&cond_start[engine]);
				pthread_cond_wait(&cond_finish[engine], &mutex_finish[engine]);
				retval = 0;
				plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
				mxGetPr(plhs[0])[0] = retval;
				pthread_mutex_unlock(&mutex_finish[engine]);
		}
		/****************************************************************************/
		else if (strcmp(cmd, "isbusy") == 0) {
				/* engine isbusy num */

				if (nrhs != 2)
						mexErrMsgTxt("Exactly two input arguments needed");

				if (!mxIsNumeric(prhs[1]))
						mexErrMsgTxt("Argument #2 should be numeric");
				engine = mxGetScalar(prhs[1]);

				if (engine < 1 || engine > poolsize)
						mexErrMsgTxt("Invalid engine number");
				engine--; /* switch from MATLAB to C indexing */

				pthread_mutex_lock(&mutex_start[engine]);
				plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
				mxGetPr(plhs[0])[0] = (status[engine] == ENGINE_BUSY);
				pthread_mutex_unlock(&mutex_start[engine]);
		}
		/****************************************************************************/
		else if (strcmp(cmd, "poolsize") == 0) {
				/* engine poolsize */

				plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
				mxGetPr(plhs[0])[0] = poolsize;
		}
		/****************************************************************************/
		else if (strcmp(cmd, "info") == 0) {
				/* engine info */

				for (engine=0; engine<poolsize; engine++) {
						pthread_mutex_lock(&mutex_start[engine]);
						command[engine] = ENGINE_INFO;
						pthread_mutex_unlock(&mutex_start[engine]);
						pthread_cond_signal(&cond_start[engine]);

				}
				retval = 0;
				plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
				mxGetPr(plhs[0])[0] = retval;
		}
		/****************************************************************************/
		else {
				mexErrMsgTxt("Unknown command");
				return;
		}

		return;
}
