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
		Engine *ep;
		pthread_t tid;
		int busy; 
		char *cmd;
		int retval; 
} engine_t;

#define STRLEN 256
#define FREE(x) {if (x) {free(x); x=NULL;}}

pthread_mutex_t enginemutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t busycondition = PTHREAD_COND_INITIALIZER;

engine_t *enginepool = NULL;
unsigned int initialized = 0;
unsigned int poolsize = 0;

/* this is called the first time that the mex-file is loaded */
void initFun(void) {
		if (!initialized) {
				initialized = 1;
				mexPrintf("engine init()\n");
		}
		return;
}

/* this function is called upon unloading of the mex-file */
void exitFun(void) {
		mexPrintf("engine cleanup()\n");
		/* FIXME list over the engines, check whether any one is busy, give warning */
		while (poolsize) {
				if (enginepool[poolsize-1].busy)
						mexWarnMsgTxt("Closing busy engine");
				if (enginepool[poolsize-1].ep!=NULL)
						engClose(enginepool[poolsize].ep);
        poolsize--;
		}
		FREE(enginepool);
		initialized = 0;
		return;
}

/* this function will be started as a seperate thread */
void *evalString(void *argin) {
		engine_t *engine;
		Engine *ep;
		char cmd[STRLEN];
		int retval;

		engine = (engine_t *)argin;
		pthread_mutex_lock(&enginemutex);
		if (engine[0].busy) {
				pthread_mutex_unlock(&enginemutex);
				pthread_exit(NULL);
				mexErrMsgTxt("The specified engine is already busy");
		}

		engine[0].busy = 1;
		ep = engine[0].ep;
		strncpy(cmd, engine[0].cmd, STRLEN-1);

		pthread_cond_signal(&busycondition);
		pthread_mutex_unlock(&enginemutex);

		retval = engEvalString(ep, cmd);

		pthread_mutex_lock(&enginemutex);
		engine[0].retval = retval;
		engine[0].busy   = 0;
#ifndef PLATFORM_WINDOWS
		engine[0].tid    = NULL;
#endif
		FREE(engine[0].cmd); 
		pthread_mutex_unlock(&enginemutex);

		pthread_exit(NULL);
		return;
}

void mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]) {
		char command[STRLEN];
		char matlabcmd[STRLEN];
		char *ptr;
		int retval, block, num, i, status;

		initFun();
		mexAtExit(exitFun);

		/* the first argument is always the command string */
		if (nrhs<1)
				mexErrMsgTxt("invalid number of input arguments");

		if (!mxIsChar(prhs[0]))
				mexErrMsgTxt("invalid input argument #1");
		if (mxGetString(prhs[0], command, STRLEN-1))
				mexErrMsgTxt("invalid input argument #1");

		/* convert to lower case */
		ptr = command;
		while (*ptr) {
				*ptr = tolower(*ptr);
				ptr++;
		}

		/****************************************************************************/
		if (strcmp(command, "open")==0) {
				/* engine open num cmd */

				if (enginepool!=NULL)
						mexErrMsgTxt("There are already engines running");

				if (nrhs<2)
						mexErrMsgTxt("Invalid number of input arguments");

				if (nrhs>1) {
						if (!mxIsNumeric(prhs[1]))
								mexErrMsgTxt("Invalid input argument #2, should be numeric");
						poolsize = mxGetScalar(prhs[1]);
				}

				if (nrhs>2) {
						if (!mxIsChar(prhs[2]))
								mexErrMsgTxt("Invalid input argument #3, should be a string");
						mxGetString(prhs[2], matlabcmd, STRLEN-1);
				}
				else {
						sprintf(matlabcmd, "matlab");
				}

				if (poolsize<1)
						mexErrMsgTxt("The number of engines in the pool should be positive");

				pthread_mutex_lock(&enginemutex);
				enginepool = (engine_t *)malloc(poolsize*sizeof(engine_t));
				status = 1;
				for (i=0; i<poolsize; i++) {
#ifdef PLATFORM_WINDOWS
                        enginepool[i].ep     = engOpen(NULL); /* returns NULL on failure */
#else
						enginepool[i].ep     = engOpen(matlabcmd); /* returns NULL on failure */
						enginepool[i].tid    = NULL;
#endif
						enginepool[i].busy   = 0;
						enginepool[i].retval = 0;
						enginepool[i].cmd    = NULL;
						status = status & (enginepool[i].ep!=NULL);
				}
				pthread_mutex_unlock(&enginemutex);
				if (!status) {
						exitFun(); /* this cleans up all engines */
						mexErrMsgTxt("failed to open MATLAB engine");
				}
		}

		/****************************************************************************/
		else if (strcmp(command, "close")==0) {
				/* engine close */

				if (enginepool==NULL)
						mexErrMsgTxt("There are no engines running");

				pthread_mutex_lock(&enginemutex);
				/* FIXME list over the engines, check whether any one is busy, give error */
				while (poolsize--) {
						if (enginepool[poolsize].busy)
								mexWarnMsgTxt("Closing busy engine");
						engClose(enginepool[poolsize].ep);
				}
				FREE(enginepool);

				retval = 0;
				plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
				mxGetPr(plhs[0])[0] = retval;

				pthread_mutex_unlock(&enginemutex);
		}

		/****************************************************************************/
		else if (strcmp(command, "put")==0) {
				/* engine put num name val */

				if (nrhs!=4)
						mexErrMsgTxt("Exactly four input arguments needed");

				if (!mxIsNumeric(prhs[1]))
						mexErrMsgTxt("Argument #2 should be numeric");
				num = mxGetScalar(prhs[1]);

				if (num<1 || num>poolsize)
						mexErrMsgTxt("Invalid engine number");

				if (!mxIsChar(prhs[2]))
						mexErrMsgTxt("Argument #3 should be a string");

				pthread_mutex_lock(&enginemutex);

				if (enginepool[num-1].busy) {
						pthread_mutex_unlock(&enginemutex);
						mexErrMsgTxt("The specified engine is busy");
				}

				retval = engPutVariable(enginepool[num-1].ep, mxArrayToString(prhs[2]), prhs[3]);

				plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
				mxGetPr(plhs[0])[0] = retval;

				pthread_mutex_unlock(&enginemutex);
		}

		/****************************************************************************/
		else if (strcmp(command, "get")==0) {
				/* engine get num name */

				if (nrhs!=3)
						mexErrMsgTxt("Exactly three input arguments needed");

				if (!mxIsNumeric(prhs[1]))
						mexErrMsgTxt("Argument #2 should be numeric");

				if (!mxIsChar(prhs[2]))
						mexErrMsgTxt("Argument #3 should be a string");

				num = mxGetScalar(prhs[1]);

				if (num<1 || num>poolsize)
						mexErrMsgTxt("Invalid engine number");

				pthread_mutex_lock(&enginemutex);

				if (enginepool[num-1].busy) {
						pthread_mutex_unlock(&enginemutex);
						mexErrMsgTxt("The specified engine is busy");
				}

				plhs[0] = engGetVariable(enginepool[num-1].ep,mxArrayToString(prhs[2]));

				pthread_mutex_unlock(&enginemutex);

		}

		/****************************************************************************/
		else if (strcmp(command, "isbusy")==0) {
				/* engine isbusy num */

				if (nrhs!=2)
						mexErrMsgTxt("Exactly two input arguments needed");

				if (!mxIsNumeric(prhs[1]))
						mexErrMsgTxt("Argument #2 should be numeric");
				num = mxGetScalar(prhs[1]);

				if (num<1 || num>poolsize)
						mexErrMsgTxt("Invalid engine number");

				pthread_mutex_lock(&enginemutex);
				retval = enginepool[num-1].busy;

				plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
				mxGetPr(plhs[0])[0] = retval;

				pthread_mutex_unlock(&enginemutex);
		}

		/****************************************************************************/
		else if (strcmp(command, "eval")==0) {
				/* engine eval num str block */

				if (nrhs<3)
						mexErrMsgTxt("At least three input arguments needed");

				if (!mxIsNumeric(prhs[1]))
						mexErrMsgTxt("Argument #2 should be numeric");
				num = mxGetScalar(prhs[1]);

				if (!mxIsChar(prhs[2]))
						mexErrMsgTxt("Invalid input argument #3, should be a string");
				mxGetString(prhs[2], matlabcmd, STRLEN-1);

				if (nrhs>3) {
						if (!mxIsNumeric(prhs[3]))
								mexErrMsgTxt("Invalid input argument #4, should be numeric");
						block = mxGetScalar(prhs[3]);
				}
				else {
						block = 0;
				}

				if (num<1 || num>poolsize)
						mexErrMsgTxt("Invalid engine number");

				pthread_mutex_lock(&enginemutex);

				if (enginepool[num-1].busy) {
						pthread_mutex_unlock(&enginemutex);
						mexErrMsgTxt("The specified engine is already busy");
				}

				if (!block) {
						enginepool[num-1].cmd = malloc(STRLEN);
						strncpy(enginepool[num-1].cmd, matlabcmd, STRLEN-1);
						retval = pthread_create(&(enginepool[num-1].tid), NULL, evalString, (void *)(enginepool+num-1));

						while (enginepool[num-1].busy==0)
								pthread_cond_wait(&busycondition, &enginemutex);
				}
				else {
						retval = engEvalString(enginepool[num-1].ep, matlabcmd);
				}

				plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
				mxGetPr(plhs[0])[0] = retval; 

				pthread_mutex_unlock(&enginemutex);
		}

		/****************************************************************************/
		else if (strcmp(command, "poolsize")==0) {
				/* engine poolsize */

				pthread_mutex_lock(&enginemutex);
				plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
				mxGetPr(plhs[0])[0] = poolsize; 
				pthread_mutex_unlock(&enginemutex);
		}

		/****************************************************************************/
		else if (strcmp(command, "info")==0) {
				/* engine info */

				pthread_mutex_lock(&enginemutex);
				for (i=0; i<poolsize; i++) {
						mexPrintf("engine[%d].ep     = %d\n",   i+1, enginepool[i].ep);
						mexPrintf("engine[%d].tid    = %d\n",   i+1, enginepool[i].tid);
						mexPrintf("engine[%d].busy   = %d\n",   i+1, enginepool[i].busy);
						mexPrintf("engine[%d].cmd    = '%s'\n", i+1, enginepool[i].cmd);
						mexPrintf("engine[%d].retval = %d\n",   i+1, enginepool[i].retval);
				}
				pthread_mutex_unlock(&enginemutex);
		}

		/****************************************************************************/
		else {
				mexErrMsgTxt("unknown command");
				return;
		}

		return;
}

