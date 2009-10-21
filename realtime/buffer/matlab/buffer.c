/*
 * Copyright (C) 2008, Robert Oostenveld
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 * $Log: buffer.c,v $
 * Revision 1.9  2008/10/22 10:45:43  roboos
 * added buffer_putdat
 *
 * Revision 1.8  2008/07/09 10:04:55  roboos
 * added event_thread, removed other acquisition devices (because sinewave and event are enough as examples)
 *
 * Revision 1.7  2008/07/08 20:31:34  roboos
 * extended the list of acquisition threads in the mex file, also added all of them to thread stopping
 *
 * Revision 1.6  2008/06/19 20:48:32  roboos
 * added support for flushing header, data and events
 *
 * Revision 1.5  2008/05/22 13:35:55  roboos
 * some small changes to make it compile on Matlab75 with borland, commented out putdat
 *
 * Revision 1.4  2008/05/22 09:25:53  roboos
 * fixed some minor issues with Borland compiler, added xxxThreadRunning variables instead of looking at thread value
 *
 * Revision 1.3  2008/03/23 17:08:25  roboos
 * improved thread accounting, implemented thread cancelation
 * added sleep_thread for testing the cancelation
 *
 * Revision 1.2  2008/03/23 12:50:10  roboos
 * multiple changes related to spawning and closing of threads
 *
 * Revision 1.1  2008/03/09 22:36:37  roboos
 * implemented new buffer mex-function, which acts as dispatcher for all functionality
 * changed buffer_xxx code from mex-function into c-function to be included in dispatcher
 *
 */

#include <pthread.h>

#include "mex.h"
#include "matrix.h"

#include "message.h"
#include "buffer.h"
#include "unix_includes.h"

#define DEFAULT_HOST "localhost"
#define DEFAULT_PORT 1972

/* these variables keep track of the threads */
pthread_t tcpserverThread;
pthread_t sleepThread;
int tcpserverThreadRunning  = 0;
int sleepThreadRunning      = 0;

/* these threads can be used to perform the acquisition by this mex file */
pthread_t sinewaveThread;
pthread_t eventThread;
int sinewaveThreadRunning   = 0;	/* usefull as demo implementation */
int eventThreadRunning      = 0;	/* usefull as demo implementation */
int acquisitionClient       = 0;	/* this counts the total number of acquisition clients, and should be 0 or 1 */

/* these are for testing the threading using "do_thread" */
#define MAX_NUM_THREADS 100
pthread_t threads[MAX_NUM_THREADS];

/* these subfunctions are used in the switch-yard below */
void buffer_gethdr(char *, int, mxArray **, const mxArray **);
void buffer_getdat(char *, int, mxArray **, const mxArray **);
void buffer_getevt(char *, int, mxArray **, const mxArray **);
void buffer_getprp(char *, int, mxArray **, const mxArray **);
void buffer_puthdr(char *, int, mxArray **, const mxArray **);
void buffer_putdat(char *, int, mxArray **, const mxArray **);
void buffer_putevt(char *, int, mxArray **, const mxArray **);
void buffer_putprp(char *, int, mxArray **, const mxArray **);

/* this is a test function for the multithreading */
void sleep_cleanup(void *arg) {
  mexPrintf("stopping sleep thread\n");
  return;
}

/* this is a test function for the multithreading */
void *sleep_thread(void *arg) {
  int oldcancelstate, oldcanceltype;
  pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, &oldcancelstate);
  pthread_setcanceltype(PTHREAD_CANCEL_DEFERRED, &oldcanceltype);
  pthread_cleanup_push(sleep_cleanup, NULL);
  mexPrintf("starting sleep thread\n");
  while (1) {
    pthread_testcancel();
    mexPrintf("sleeping for 4 seconds...\n");
    sleep(4);
  }
  pthread_cleanup_pop(1);
  return NULL;
}

/* this is a test function for the multithreading */
void *do_thread(void *threadid) {
  int tid;
  tid = (int)threadid;
  mexPrintf("In thread: just started thread #%d\n", tid);
  sleep(4); /* pretend that the code is busy for a few seconds */
  mexPrintf("In thread: about to exit from thread #%d\n", tid);
  pthread_exit(NULL);
  return NULL;
}

/* this function is called upon unloading of the mex-file */
void cleanupMex(void)
{
  if (tcpserverThreadRunning) {
    mexWarnMsgTxt("stopping tcpserver thread");
    tcpserverThreadRunning = 0;
  }
  if (sleepThreadRunning) {
    mexWarnMsgTxt("stopping sleep thread");
    sleepThreadRunning = 0;
  }
  if (sinewaveThreadRunning) {
    mexWarnMsgTxt("stopping sinewave thread");
    acquisitionClient--;
    sinewaveThreadRunning = 0;
  }
  if (eventThreadRunning) {
    mexWarnMsgTxt("stopping event thread");
    acquisitionClient--;
    eventThreadRunning = 0;
  }
  return;
}

void mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
  int port;
  char *command = NULL, *argument = NULL, *hostname = NULL;
  int num;
  int rc, t;
  host_t *host; /* WARNING this is passed on to the thread */
  
  if (nrhs <2)
    mexErrMsgTxt ("Invalid number of input arguments");
  
  if (nrhs>0) {
    /* first argument is command, e.g. GET_HDR or TCPSERVER */
    command = mxArrayToString(prhs[0]);
    if (command==NULL)
      mexErrMsgTxt ("invalid input argument #1");
  }
  
  if (nrhs>1) {
    /* second argument depends on the command */
    /* this will be processed below           */
  }
  
  if (nrhs>2) {
    /* third argument is the hostname (optional) */
    hostname = mxArrayToString(prhs[2]);
    if (hostname==NULL)
      mexErrMsgTxt ("invalid input argument #3");
  }
  else {
    /* use the default value */
    hostname = malloc(strlen(DEFAULT_HOST)+1);
    strncpy(hostname, DEFAULT_HOST, strlen(DEFAULT_HOST)+1);
  }
  
  if (nrhs>3) {
    /* fourth argument is the port number (optional) */
    if (mxIsEmpty(prhs[3]) || !mxIsNumeric(prhs[3]))
      mexErrMsgTxt ("invalid input argument #4");
    else
      port = (int)mxGetScalar(prhs[3]);
  }
  else {
    /* use the default value */
    port = DEFAULT_PORT;
  }
  
  /* FIXME clean up the various host/port representations */
  /* FIXME this leaks some memory, perhaps call free in the thread? */
  host = malloc(sizeof(host_t));
  sprintf(host->name, "%s", hostname);
  host->port = port;
  
  if (strcasecmp(command, "tcpserver")==0) {
    argument = mxArrayToString(prhs[1]);
    if (argument==NULL)
      mexErrMsgTxt ("invalid input argument #2");
    if (strcasecmp(argument, "init")==0) {
      if (tcpserverThreadRunning)
        mexErrMsgTxt("thread is already running");
      mexPrintf("In main: spawning tcpserver thread\n");
      rc = pthread_create(&tcpserverThread, NULL, tcpserver, (void *)host);
      if (rc)
        mexErrMsgTxt("problem with return code from pthread_create()");
      else
        tcpserverThreadRunning = 1;
      
    }
    else if (strcasecmp(argument, "exit")==0) {
      if (!tcpserverThreadRunning)
        mexErrMsgTxt("thread is not running");
      mexPrintf("In main: requesting cancelation of tcpserver thread\n");
      rc = pthread_cancel(tcpserverThread);
      if (rc)
        mexErrMsgTxt("problem with return code from pthread_cancel()");
      else
        tcpserverThreadRunning = 0;
    }
  }
  
  else if (strcasecmp(command, "sleep")==0) {
    argument = mxArrayToString(prhs[1]);
    if (argument==NULL)
      mexErrMsgTxt ("invalid input argument #2");
    if (strcasecmp(argument, "init")==0) {
      if (sleepThreadRunning)
        mexErrMsgTxt("thread is already running");
      mexPrintf("In main: spawning sleep thread\n");
      rc = pthread_create(&sleepThread, NULL, sleep_thread, NULL);
      if (rc)
        mexErrMsgTxt("problem with return code from pthread_create()");
      else
        sleepThreadRunning = 1;
    }
    else if (strcasecmp(argument, "exit")==0) {
      if (!sleepThreadRunning)
        mexErrMsgTxt("thread is already running");
      mexPrintf("In main: requesting cancelation of sleep thread\n");
      rc = pthread_cancel(sleepThread);
      if (rc)
        mexErrMsgTxt("problem with return code from pthread_cancel()");
      else
        sleepThreadRunning = 0;
    }
  }
  
  else if (strcasecmp(command, "sinewave")==0) {
    argument = mxArrayToString(prhs[1]);
    if (argument==NULL)
      mexErrMsgTxt ("invalid input argument #2");
    if (strcasecmp(argument, "init")==0) {
      if (acquisitionClient>0)
        mexErrMsgTxt("cannot create more than one acquisition client");
      if (sinewaveThreadRunning)
        mexErrMsgTxt("thread is already running");
      mexPrintf("In main: spawning sinewave thread\n");
      rc = pthread_create(&sinewaveThread, NULL, sinewave_thread, (void *)host);
      if (rc)
        mexErrMsgTxt("problem with return code from pthread_create()");
      else {
        sinewaveThreadRunning = 1;
        acquisitionClient++;
      }
    }
    else if (strcasecmp(argument, "exit")==0) {
      if (!sinewaveThreadRunning)
        mexErrMsgTxt("thread is not running");
      mexPrintf("In main: requesting cancelation of sinewave thread\n");
      rc = pthread_cancel(sinewaveThread);
      if (rc)
        mexErrMsgTxt("problem with return code from pthread_cancel()");
      else {
        sinewaveThreadRunning = 0;
        acquisitionClient--;
      }
    }
  }
  
  else if (strcasecmp(command, "event")==0) {
    argument = mxArrayToString(prhs[1]);
    if (argument==NULL)
      mexErrMsgTxt ("invalid input argument #2");
    if (strcasecmp(argument, "init")==0) {
      if (acquisitionClient>0)
        mexErrMsgTxt("cannot create more than one acquisition client");
      if (eventThreadRunning)
        mexErrMsgTxt("thread is already running");
      mexPrintf("In main: spawning event thread\n");
      rc = pthread_create(&eventThread, NULL, event_thread, (void *)host);
      if (rc)
        mexErrMsgTxt("problem with return code from pthread_create()");
      else {
        eventThreadRunning = 1;
        acquisitionClient++;
      }
    }
    else if (strcasecmp(argument, "exit")==0) {
      if (!eventThreadRunning)
        mexErrMsgTxt("thread is not running");
      mexPrintf("In main: requesting cancelation of event thread\n");
      rc = pthread_cancel(eventThread);
      if (rc)
        mexErrMsgTxt("problem with return code from pthread_cancel()");
      else {
        eventThreadRunning = 0;
        acquisitionClient--;
      }
    }
  }
  
  else if (strcasecmp(command, "get_hdr")==0) {
    buffer_gethdr(hostname, port, &(plhs[0]), &(prhs[1]));
  }
  
  else if (strcasecmp(command, "get_dat")==0) {
    buffer_getdat(hostname, port, &(plhs[0]), &(prhs[1]));
  }
  
  else if (strcasecmp(command, "get_evt")==0) {
    buffer_getevt(hostname, port, &(plhs[0]), &(prhs[1]));
  }
  
  else if (strcasecmp(command, "get_prp")==0) {
    buffer_getprp(hostname, port, &(plhs[0]), &(prhs[1]));
  }
  
  else if (strcasecmp(command, "put_hdr")==0) {
    buffer_puthdr(hostname, port, &(plhs[0]), &(prhs[1]));
  }
  
  else if (strcasecmp(command, "put_dat")==0) {
    buffer_putdat(hostname, port, &(plhs[0]), &(prhs[1]));
  }
  
  else if (strcasecmp(command, "put_evt")==0) {
    buffer_putevt(hostname, port, &(plhs[0]), &(prhs[1]));
  }
  
  else if (strcasecmp(command, "put_prp")==0) {
    buffer_putprp(hostname, port, &(plhs[0]), &(prhs[1]));
  }
  
  else if (strcasecmp(command, "flush_hdr")==0) {
    buffer_flushhdr(hostname, port, &(plhs[0]), &(prhs[1]));
  }
  
  else if (strcasecmp(command, "flush_dat")==0) {
    buffer_flushdat(hostname, port, &(plhs[0]), &(prhs[1]));
  }
  
  else if (strcasecmp(command, "flush_evt")==0) {
    buffer_flushevt(hostname, port, &(plhs[0]), &(prhs[1]));
  }
  
  else if (strcasecmp(command, "do_thread")==0) {
    /* get the second input argument, should be single number */
    if (mxIsEmpty(prhs[1]) || !mxIsNumeric(prhs[1]))
      mexErrMsgTxt ("invalid input argument #2");
    else
      num = mxGetScalar(prhs[1]);
    
    if (num>MAX_NUM_THREADS)
      mexErrMsgTxt ("number of requested threads too large");
    
    for(t=0; t<num; t++){
      mexPrintf("In main: spawning thread %d\n", t);
      rc = pthread_create(&threads[t], NULL, do_thread, (void *)t);
      if (rc){
        mexErrMsgTxt("problem with return code from pthread_create()");
      }
      /* sleep for one seconds */
      sleep(1);
    }
  }
  
  else {
    mexErrMsgTxt ("unknown command");
  }
  
  /* FIXME freeing the ones returned by mx... gives an error
  FREE(command);
  FREE(argument);
  FREE(hostname);
   */
  
  return;
}

