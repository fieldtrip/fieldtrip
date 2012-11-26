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

#define DEBUG_PRINT(string, ...) {mexPrintf(string, ## __VA_ARGS__);}
#define FREE(x)                  {if (x) {free(x); x=NULL;}}
#define TRUE            1
#define FALSE           0
#define MAX_THREADS     256
#define STRLEN          256

/************************************************************************/
#ifdef USE_PTHREADS
#include <pthread.h>

pthread_mutex_t mutex_start [MAX_THREADS];
pthread_mutex_t mutex_finish[MAX_THREADS];
pthread_cond_t  cond_start	[MAX_THREADS];
pthread_cond_t  cond_finish	[MAX_THREADS];
pthread_t       threadid    [MAX_THREADS];
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

CRITICAL_SECTION    mutex_start [MAX_THREADS];
CRITICAL_SECTION    mutex_finish[MAX_THREADS];
CONDITION_VARIABLE  cond_start	[MAX_THREADS];
CONDITION_VARIABLE  cond_finish	[MAX_THREADS];
HANDLE              threadid    [MAX_THREADS];
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

#endif
/************************************************************************/

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

unsigned int    status[MAX_THREADS];
unsigned int    command[MAX_THREADS];
unsigned int    poolsize      = 0;
unsigned int    initialized   = FALSE;
unsigned int    retval        = 0;
void           *ptr1          = NULL;
void           *ptr2          = NULL;
char            matlabcmd[STRLEN];

/* this is called the first time that the mex-file is loaded */
void initFun(void)
{
    unsigned int engine;
    if (initialized==FALSE) {
        DEBUG_PRINT("Init\n");
        for (engine=0; engine<MAX_THREADS; engine++) {
            status[engine] = ENGINE_UNKNOWN;
            command[engine] = ENGINE_UNKNOWN;
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
            command[engine] = ENGINE_CLOSE;
            MUTEX_UNLOCK(&mutex_start[engine]);
            COND_SIGNAL(&cond_start[engine]);
            THREAD_JOIN(threadid[engine]);
        }
        for (engine=0; engine<MAX_THREADS; engine++) {
            status[engine] = ENGINE_UNKNOWN;
            command[engine] = ENGINE_UNKNOWN;
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

/* this function will be started as a seperate thread */
void engineThread(void *argin)
{
    Engine         *ep = NULL;
    char            str[STRLEN];
    int             thisId, thisCommand, retval;
    
    if (!mexIsLocked()) {
        DEBUG_PRINT("Locking mex file\n");
        mexLock();
    }
    
    /* store the identifier for this thread, it is used for indexing the global arrays */
    thisId = (int)argin;
    
    DEBUG_PRINT("Starting engine thread %d\n", thisId);
    
    /* start the MATLAB engine, this thread will remain responsible for it */
#ifdef PLATFORM_WINDOWS
ep = engOpenSingleUse(NULL, NULL, &retval);
#else
ep = engOpen(matlabcmd);
retval = (ep!=NULL);
#endif

COND_SIGNAL(&cond_finish[thisId]);

while (1) {
    DEBUG_PRINT("Entering while loop in engine thread %d\n", thisId);
    MUTEX_LOCK(&mutex_start[thisId]);
    status[thisId]  = ENGINE_IDLE;
    
    /* wait until the start condition is signalled */
    COND_WAIT(&cond_start[thisId], &mutex_start[thisId]);
    
    thisCommand     = command[thisId];
    status[thisId]  = ENGINE_BUSY;
    MUTEX_UNLOCK(&mutex_start[thisId]);
    
    DEBUG_PRINT("Engine thread executing %d in %d\n", thisCommand, thisId);
    retval          = 0;
    
    if (thisCommand==ENGINE_CLOSE) {
        /* close the engine belonging to this thread and exit the thread */
        retval = engClose(ep);
        status[thisId] = ENGINE_UNKNOWN;
        COND_SIGNAL(&cond_finish[thisId]);
        THREAD_EXIT(NULL);
        return;
    }
    else if (thisCommand==ENGINE_PUT) {
        /* copy a variable from the local matlab workspace to the remote engine */
        retval = engPutVariable(ep, mxArrayToString(ptr1), ptr2);
        DEBUG_PRINT("retval = %d\n", retval);
        COND_SIGNAL(&cond_finish[thisId]);
    }
    else if (thisCommand==ENGINE_GET) {
        /* copy a variable from the remote engine workspace to the local matlab */
        ptr2 = engGetVariable(ep, mxArrayToString(ptr1));
        retval = (ptr2==NULL);
        DEBUG_PRINT("retval = %d\n", retval);
        COND_SIGNAL(&cond_finish[thisId]);
    }
    else if (thisCommand==ENGINE_EVAL) {
        /* evaluate the string on the remote engine */
        memset(str, 0, STRLEN);
        strncpy(str, mxArrayToString(ptr1), STRLEN);
        COND_SIGNAL(&cond_finish[thisId]);
        DEBUG_PRINT("Starting evaluation of '%s')\n", str);
        retval = engEvalString(ep, str);
        DEBUG_PRINT("Finished evaluation, retval = %d\n", retval);
    }
    else if (thisCommand==ENGINE_INFO) {
        /* give some information about the engine */
        DEBUG_PRINT("Information for engine %d, retval = %d\n", thisId, retval);
        COND_SIGNAL(&cond_finish[thisId]);
    }
    else {
        COND_SIGNAL(&cond_finish[thisId]);
    } /* if operation */
    
} /* while */
}

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
    char            str[STRLEN];
    char           *ptr;
    int             retval, block, engine;
    
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
			COND_WAIT(&cond_finish[engine], &mutex_finish[engine]);
			DEBUG_PRINT("The engine thread finished\n");
			MUTEX_UNLOCK(&mutex_finish[engine]);

			DEBUG_PRINT("Started engine %d, retval = %d\n", engine+1, retval);
			if (retval)
			    break;
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
        if (status[engine]!=ENGINE_IDLE) {
            MUTEX_UNLOCK(&mutex_finish[engine]);
            mexErrMsgTxt("The specified engine is not available");
        }
        
        ptr1 = (void *)prhs[2]; /* name */
        ptr2 = (void *)prhs[3]; /* value */
        command[engine] = ENGINE_PUT;
        
        DEBUG_PRINT("Signalling thread for ENGINE_PUT\n");
        COND_SIGNAL(&cond_start[engine]);
        COND_WAIT(&cond_finish[engine], &mutex_finish[engine]);
        DEBUG_PRINT("The engine thread finished\n");
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
        if (status[engine]!=ENGINE_IDLE) {
            MUTEX_UNLOCK(&mutex_finish[engine]);
            mexErrMsgTxt("The specified engine is not available");
        }
        
        ptr1 = (void *)prhs[2]; /* name */
        command[engine] = ENGINE_GET;
        
        DEBUG_PRINT("Signalling start for ENGINE_GET\n");
        COND_SIGNAL(&cond_start[engine]);
        COND_WAIT(&cond_finish[engine], &mutex_finish[engine]);
        DEBUG_PRINT("The engine thread finished\n");
        plhs[0] = (mxArray *)ptr2;
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
        if (status[engine]!=ENGINE_IDLE) {
            MUTEX_UNLOCK(&mutex_finish[engine]);
            mexErrMsgTxt("The specified engine is not available");
        }
        
        ptr1 = (void *)prhs[2]; /* command */
        command[engine] = ENGINE_EVAL;
        
        DEBUG_PRINT("Signalling start for ENGINE_EVAL\n");
        COND_SIGNAL(&cond_start[engine]);
        COND_WAIT(&cond_finish[engine], &mutex_finish[engine]);
        DEBUG_PRINT("The engine thread returned control\n");
        MUTEX_UNLOCK(&mutex_finish[engine]);
    }
    /****************************************************************************/
    else if (strcmp(str, "info") == 0) {
        /* engine info */
        
        for (engine=0; engine<poolsize; engine++) {
            
            MUTEX_LOCK(&mutex_finish[engine]);
            if (status[engine]!=ENGINE_IDLE) {
                MUTEX_UNLOCK(&mutex_finish[engine]);
                mexErrMsgTxt("The specified engine is not available");
            }
            
            command[engine] = ENGINE_INFO;
            
            DEBUG_PRINT("Signalling start for ENGINE_INFO\n");
            COND_SIGNAL(&cond_start[engine]);
            COND_WAIT(&cond_finish[engine], &mutex_finish[engine]);
            DEBUG_PRINT("The engine thread finished\n");
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
        
        plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
        mxGetPr(plhs[0])[0] = (status[engine] == ENGINE_BUSY);
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
    
    return;
}
