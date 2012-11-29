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

#include "mex.h"
#include "engine.h"
#include "platform.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#ifdef USE_PTHREADS
#include <pthread.h>
#else
#include <windows.h>
#endif

#define STRLEN 256
#define FREE(x) {if (x) {free(x); x=NULL;}}

/**************************************************************************************************/
#ifdef USE_PTHREADS
pthread_mutex_t enginemutex   = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t  busycond      = PTHREAD_COND_INITIALIZER;
#define ENGINEMUTEX_LOCK        pthread_mutex_lock(&enginemutex)
#define ENGINEMUTEX_UNLOCK      pthread_mutex_unlock(&enginemutex)
#define BUSYCOND_WAIT           pthread_cond_wait(&busycond, &enginemutex)
#define BUSYCOND_RELEASE        pthread_cond_signal(&busycond)
#define THREAD_EXIT             pthread_exit(NULL);

#else
CRITICAL_SECTION enginemutex;
#define ENGINEMUTEX_LOCK   EnterCriticalSection(&enginemutex)
#define ENGINEMUTEX_UNLOCK LeaveCriticalSection(&enginemutex)
#define BUSYCOND_WAIT
#define BUSYCOND_RELEASE
#define THREAD_EXIT

#endif
/**************************************************************************************************/

Engine         *engine = NULL;
unsigned int    initialized = 0;

/* this is called the first time that the mex-file is loaded */
void initFun(void) {
  if (!initialized) {
    ENGINEMUTEX_LOCK;
    initialized = 1;
    ENGINEMUTEX_UNLOCK;
    mexPrintf("engine init()\n");
  }
  return;
}

/* this function is called upon unloading of the mex-file */
void exitFun(void) {
  ENGINEMUTEX_LOCK;
  initialized = 0;
  ENGINEMUTEX_UNLOCK;
  mexPrintf("engine exit()\n");
  return;
}

/* this function will be started as a seperate thread */
void evalString(void *argin) {
  int   retval;
  char *cmd;
  
  mexPrintf("Starting thread\n");
  ENGINEMUTEX_LOCK;
  
  cmd = (char *) argin;
  
  if (engine!=0) {
    mexPrintf("Executing command \"%s\" on %d\n", cmd, engine);
    retval = engEvalString(engine, cmd);
    mexPrintf("Finished executing command, retval = %d\n", retval);
  }
  
  ENGINEMUTEX_UNLOCK;
  
  mexPrintf("Finished thread\n");
  THREAD_EXIT;
  return;
}

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]) {
  char            command[STRLEN];
  char            matlabcmd[STRLEN];
  char            cmd[STRLEN];
  char           *ptr;
  int             poolsize, retval, block, num, i, status;
  #ifndef USE_PTHREADS
          HANDLE TName;
  InitializeCriticalSection(&enginemutex);
  #endif
          
          initFun();
  mexAtExit(exitFun);
  
  /* the first argument is always the command string */
  if (nrhs < 1)
    mexErrMsgTxt("invalid number of input arguments");
  
  if (!mxIsChar(prhs[0]))
    mexErrMsgTxt("invalid input argument #1");
  if (mxGetString(prhs[0], command, STRLEN - 1))
    mexErrMsgTxt("invalid input argument #1");
  
  /* convert to lower case */
  ptr = command;
  while (*ptr) {
    *ptr = tolower(*ptr);
    ptr++;
  }
  
  /****************************************************************************/
  if (strcmp(command, "open") == 0) {
    /* engine open num cmd */
    
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
    } else {
      sprintf(matlabcmd, "matlab");
    }
    
    if (poolsize < 1)
      mexErrMsgTxt("The number of engines in the pool should be positive");
    
    ENGINEMUTEX_LOCK;
    engine = engOpenSingleUse(NULL, NULL, &retval);	/* returns NULL on failure */
    ENGINEMUTEX_UNLOCK;
    
    if (!engine) {
      exitFun();	/* this cleans up all engines */
      mexErrMsgTxt("failed to open MATLAB engine");
    }
  }
  /****************************************************************************/
  else if (strcmp(command, "close") == 0) {
    /* engine close */
    
    exitFun();
    
    retval = 0;
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxGetPr(plhs[0])[0] = retval;
    
  }
  /****************************************************************************/
  else if (strcmp(command, "put") == 0) {
    /* engine put num name val */
    
    if (nrhs != 4)
      mexErrMsgTxt("Exactly four input arguments needed");
    
    if (!mxIsNumeric(prhs[1]))
      mexErrMsgTxt("Argument #2 should be numeric");
    num = mxGetScalar(prhs[1]);
    
    if (num < 1 || num > poolsize)
      mexErrMsgTxt("Invalid engine number");
    
    if (!mxIsChar(prhs[2]))
      mexErrMsgTxt("Argument #3 should be a string");
    
    ENGINEMUTEX_LOCK;
    
    retval = engPutVariable(engine, mxArrayToString(prhs[2]), prhs[3]);
    
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxGetPr(plhs[0])[0] = retval;
    
    ENGINEMUTEX_UNLOCK;
  }
  /****************************************************************************/
  else if (strcmp(command, "get") == 0) {
    /* engine get num name */
    
    if (nrhs != 3)
      mexErrMsgTxt("Exactly three input arguments needed");
    
    if (!mxIsNumeric(prhs[1]))
      mexErrMsgTxt("Argument #2 should be numeric");
    
    if (!mxIsChar(prhs[2]))
      mexErrMsgTxt("Argument #3 should be a string");
    
    num = mxGetScalar(prhs[1]);
    
    ENGINEMUTEX_LOCK;
    
    plhs[0] = engGetVariable(engine, mxArrayToString(prhs[2]));
    
    ENGINEMUTEX_UNLOCK;
    
  }
  /****************************************************************************/
  else if (strcmp(command, "isbusy") == 0) {
    /* engine isbusy num */
    
    if (nrhs != 2)
      mexErrMsgTxt("Exactly two input arguments needed");
    
    if (!mxIsNumeric(prhs[1]))
      mexErrMsgTxt("Argument #2 should be numeric");
    num = mxGetScalar(prhs[1]);
    
    ENGINEMUTEX_LOCK;
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxGetPr(plhs[0])[0] = 1;
    ENGINEMUTEX_UNLOCK;
  }
  /****************************************************************************/
  else if (strcmp(command, "eval") == 0) {
    /* engine eval num str block */
    
    if (nrhs < 3)
      mexErrMsgTxt("At least three input arguments needed");
    
    if (!mxIsNumeric(prhs[1]))
      mexErrMsgTxt("Argument #2 should be numeric");
    num = mxGetScalar(prhs[1]);
    
    if (!mxIsChar(prhs[2]))
      mexErrMsgTxt("Invalid input argument #3, should be a string");
    mxGetString(prhs[2], cmd, STRLEN - 1);
    
    if (nrhs > 3) {
      if (!mxIsNumeric(prhs[3]))
        mexErrMsgTxt("Invalid input argument #4, should be numeric");
      block = mxGetScalar(prhs[3]);
    } else {
      block = 0;
    }
    
    ENGINEMUTEX_LOCK;
    
    if (!block) {
      
#ifdef USE_PTHREADS
      retval = pthread_create(&(enginepool[num - 1].tid), NULL, (void *) &evalString, (void *) (enginepool + num - 1));
#else
      TName = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE) evalString, (void *) (cmd), 0, NULL);
#endif

    } else {
      retval = engEvalString(engine, cmd);
    }
    mexPrintf("Finished executing command, retval = %d\n", retval);
    ENGINEMUTEX_UNLOCK;
    
    WaitForSingleObject(TName, INFINITE);
    
    /* wait for the thread to become busy */
    mexPrintf("Waiting for engine to become busy\n");
    BUSYCOND_WAIT;
    mexPrintf("The engine has become busy\n");
    
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxGetPr(plhs[0])[0] = retval;
  }
  /****************************************************************************/
  else if (strcmp(command, "poolsize") == 0) {
    /* engine poolsize */
    
    ENGINEMUTEX_LOCK;
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxGetPr(plhs[0])[0] = (engine!=0);
    ENGINEMUTEX_UNLOCK;
  }
  /****************************************************************************/
  else if (strcmp(command, "info") == 0) {
    /* engine info */
    
    ENGINEMUTEX_LOCK;
    mexPrintf("engine = %d\n",   engine);
    ENGINEMUTEX_UNLOCK;
  }
  /****************************************************************************/
  else {
    mexErrMsgTxt("unknown command");
    return;
  }
  
  return;
}
