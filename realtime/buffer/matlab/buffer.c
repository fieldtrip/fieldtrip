/*
 * Copyright (C) 2008, Robert Oostenveld
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 */

#include <pthread.h>

#include "mex.h"
#include "matrix.h"
#include "buffer.h"
#include <string.h>

#define DEFAULT_HOST "localhost"
#define DEFAULT_PORT 1972


#ifdef WIN32
#define sleep(x)   Sleep((x)*1000)
#endif


/* This struct is used for keeping a linked list of hostnames / ports / sockets.
   The list is initially empty and added to whenever a new connection is opened.
   sock>0 indicates an open TCP connection.
   sock=0 corresponds to a local fieldtrip buffer spawned within this MEX file.
   The list is kept until the MEX-file is unloaded and automatically cleaned up 
   inside the MEX-file exit routine. TCP communication errors on a specific socket
   will trigger closing the connection and removing the corresponding list item.
*/
typedef struct host_port_sock_list_item {
	char *hostname;
	int port;
	int sock; 	 							
	struct host_port_sock_list_item *next;	/* NULL if last elemenent  */
} host_port_sock_list_item_t;

/* This is the head of the list */
host_port_sock_list_item_t *firstHostPortSock = NULL;

/* these variables keep track of the threads */
pthread_t tcpserverThread;
pthread_t sleepThread;
int tcpserverThreadRunning  = 0;
int sleepThreadRunning      = 0;
/* a pointer to this is passed on to the thread, we only need it once
   and therefore use a static variable */
host_t host_for_server; 

/* these are for testing the threading using "do_thread" */
#define MAX_NUM_THREADS 100
pthread_t threads[MAX_NUM_THREADS];

/* these subfunctions are used in the switch-yard below 
** SK: they now return integer error codes like the following
**  <0 : some low-level communication error occured (especially -3 from the tcprequest)
**   0 : no error
**  >0 : an error on the buffer level occured (as returned from the buffer
** To keep track of connection failures (for closing the sockets), 
** mexErrMsgTxt should only be used to report pure Matlab errors (e.g., invalid 
** arguments).
*/
int buffer_gethdr(int, mxArray **, const mxArray **);
int buffer_getdat(int, mxArray **, const mxArray **);
int buffer_getevt(int, mxArray **, const mxArray **);
int buffer_getprp(int, mxArray **, const mxArray **);
int buffer_puthdr(int, mxArray **, const mxArray **);
int buffer_putdat(int, mxArray **, const mxArray **);
int buffer_putevt(int, mxArray **, const mxArray **);
int buffer_putprp(int, mxArray **, const mxArray **);
int buffer_flushhdr(int, mxArray **, const mxArray **);
int buffer_flushdat(int, mxArray **, const mxArray **);
int buffer_flushevt(int, mxArray **, const mxArray **);



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
void cleanupMex(void) {
	int verbose = 1;
  
	if (verbose) {
		printf("Entering cleanupMex() routine\n");
	}
  
	if (tcpserverThreadRunning) {
		mexWarnMsgTxt("stopping tcpserver thread");
		tcpserverThreadRunning = 0;
	}
	if (sleepThreadRunning) {
		mexWarnMsgTxt("stopping sleep thread");
		sleepThreadRunning = 0;
	}
	/* clean up host/address/socket list and close open sockets */
	while (firstHostPortSock != NULL) {
		host_port_sock_list_item_t *hpsli = firstHostPortSock;

		if (hpsli->sock > 0) {
			if (verbose) {
				printf("Closing socket and ");
			}
			close_connection(hpsli->sock);
		}
		if (verbose) {
			printf("cleaning up list entry %s:%i\n",hpsli->hostname, hpsli->port);
		}
		FREE(hpsli->hostname);

		firstHostPortSock = hpsli->next;
		free(hpsli);
	}
	return;
}

/** This function searches through the linked list to retrieve the socket number
	of the first (and only) item with matching hostname and port. 
	Returns -1 if no item matched.
*/
int lookup_hps_item(const char *hostname, int port) {
	host_port_sock_list_item_t *hpsli = firstHostPortSock;
	
	while (hpsli != NULL) {
		if ((strcmp(hpsli->hostname, hostname)==0) && (hpsli->port == port)) {
			return hpsli->sock;
		}
		hpsli = hpsli->next;
	}
	return -1;
}

/* This function adds a new item to the host/port/socket list without checking if
	the same element is already there. Also registers the cleanup routine.
*/
int add_hps_item(const char *hostname, int port, int sock) {
	host_port_sock_list_item_t *hpsli;
	int n;
	
	if (hostname == NULL) return 0;
	n = strlen(hostname);

		
	hpsli = (host_port_sock_list_item_t *) malloc(sizeof(host_port_sock_list_item_t));
	if (hpsli == NULL) return 0;	/* out of memory - probably never */
	
	hpsli->hostname = (char *) malloc(n+1);
	if (hpsli->hostname == NULL) {
		/* out of memory - probably never */
		free(hpsli);
		return 0;    
	}
	memcpy(hpsli->hostname, hostname, n+1);
	
	hpsli->port = port;
	hpsli->sock = sock;
	hpsli->next = firstHostPortSock;
	firstHostPortSock = hpsli;
	
	mexAtExit(cleanupMex);   /* register cleanup routine so the list get's properly destroyed */
	return 1;
}

/** This is a MEX-specific wrapper for the open_connection call.
	If the requested host+port combination is already in the linked list,
	it means that the connection is still open (at least from this side), and
	we can just return the corresponding socket.
	Otherwise, we try to connect using the normal open_connection, and
	in case of success, add the new host/port/sock triple to the list.
*/
int open_connection_with_list(const char *hostname, int port) {
	int verbose = 0;
	int sock;
			
	/* first perform lookup of already known hostnames and connections */
	if (verbose > 0) 
		printf("Looking for list item...\n");
	sock = lookup_hps_item(hostname, port);
	/* hostname/port combination found? just return the socket */
	if (sock >= 0) {
		if (verbose > 0) 
			printf("Found connection on %i\n",sock);
		return sock;
	}
	
	if (verbose>0) 
		printf("Trying to open new connection...\n");

	sock = open_connection(hostname, port);
	if (sock>0) {
		if (verbose>0)
			printf("open_connection: connected to %s:%d on socket %d\n", hostname, port, sock);
		if (!add_hps_item(hostname, port, sock)) {
			/* out of memory for creating a list entry - IF this ever happens, we better close the socket right away */
			close_connection(sock);
			/* return with an error further down */
		} else {
			return sock;
		}
	}
	mexErrMsgTxt("ERROR: failed to create socket (1)");
	return -1;
}

/* This function searches the item with matching socket number and removes it from the list 
	Note: The socket itself is not closed here.
*/
void remove_hps_item(int sock) {
	host_port_sock_list_item_t *hpsli = firstHostPortSock;
	/* the "next" pointer of the previous list item, or initially a pointer to the head of the list */
	host_port_sock_list_item_t **prev_next = &firstHostPortSock;   
	
	while (hpsli != NULL) {
		if (hpsli->sock == sock) {
			*prev_next = hpsli->next;
			free(hpsli->hostname);			
			free(hpsli);
			return;
		}
		prev_next = &(hpsli->next);
		hpsli = hpsli->next;
	}
}

/* for debugging */
void dump_hps_list() {
	host_port_sock_list_item_t *hpsli = firstHostPortSock;
	
	printf("----- Open connections -----\n");
	while (hpsli != NULL) {
		printf("%s:%i -> %i\n",hpsli->hostname, hpsli->port, hpsli->sock);
		hpsli = hpsli->next;
	}
	printf("----------------------------\n");
}


void mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
  int port;
  char *command = NULL, *argument = NULL, *hostname = NULL;
  int num;
  int rc, t;
  int server, errorCode = 0;
  
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
    hostname = (char *) mxMalloc(strlen(DEFAULT_HOST)+1);
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
      
  if (strcasecmp(command, "tcpserver")==0) {
    argument = mxArrayToString(prhs[1]);
    if (argument==NULL)
      mexErrMsgTxt ("invalid input argument #2");
    if (strcasecmp(argument, "init")==0) {
		
      if (tcpserverThreadRunning)
        mexErrMsgTxt("thread is already running");
      mexPrintf("In main: spawning tcpserver thread\n");
	  
      strncpy(host_for_server.name, hostname, 256);  /* the 256 is a limit inside the host_t structure */
      host_for_server.port = port;
	  
      rc = pthread_create(&tcpserverThread, NULL, tcpserver, (void *)&host_for_server);
      if (rc)
        mexErrMsgTxt("problem with return code from pthread_create()");
		
      tcpserverThreadRunning = 1;

	  /* put this as a "connection" with sock=0 into the list */
	  add_hps_item(hostname, port, 0);
    }
    else if (strcasecmp(argument, "exit")==0) {
      if (!tcpserverThreadRunning)
        mexErrMsgTxt("thread is not running");
      mexPrintf("In main: requesting cancelation of tcpserver thread\n");
      rc = pthread_cancel(tcpserverThread);
      if (rc)
        mexErrMsgTxt("problem with return code from pthread_cancel()");
      tcpserverThreadRunning = 0;
	  /* remove sock=0 connection from host/port/socket list */
	  remove_hps_item(0);
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
  
  else if (strcasecmp(command, "get_hdr")==0) {
	server = open_connection_with_list(hostname, port);
    errorCode = buffer_gethdr(server, &(plhs[0]), &(prhs[1]));
  }
  
  else if (strcasecmp(command, "get_dat")==0) {
    server = open_connection_with_list(hostname, port);
    errorCode = buffer_getdat(server, &(plhs[0]), &(prhs[1]));
  }
  
  else if (strcasecmp(command, "get_evt")==0) {
    server = open_connection_with_list(hostname, port);
    errorCode = buffer_getevt(server, &(plhs[0]), &(prhs[1]));
  }
  
  else if (strcasecmp(command, "get_prp")==0) {
    server = open_connection_with_list(hostname, port);
    errorCode = buffer_getprp(server, &(plhs[0]), &(prhs[1]));
  }
  
  else if (strcasecmp(command, "put_hdr")==0) {
    server = open_connection_with_list(hostname, port);
    errorCode = buffer_puthdr(server, &(plhs[0]), &(prhs[1]));
  }
  
  else if (strcasecmp(command, "put_dat")==0) {
    server = open_connection_with_list(hostname, port);
    errorCode = buffer_putdat(server, &(plhs[0]), &(prhs[1]));
  }
  
  else if (strcasecmp(command, "put_evt")==0) {
    server = open_connection_with_list(hostname, port);
    errorCode = buffer_putevt(server, &(plhs[0]), &(prhs[1]));
  }
  
  else if (strcasecmp(command, "put_prp")==0) {
    server = open_connection_with_list(hostname, port);
    errorCode = buffer_putprp(server, &(plhs[0]), &(prhs[1]));
  }
  
  else if (strcasecmp(command, "flush_hdr")==0) {
    server = open_connection_with_list(hostname, port);
    errorCode = buffer_flushhdr(server, &(plhs[0]), &(prhs[1]));
  }
  
  else if (strcasecmp(command, "flush_dat")==0) {
    server = open_connection_with_list(hostname, port);
    errorCode = buffer_flushdat(server, &(plhs[0]), &(prhs[1]));
  }
  
  else if (strcasecmp(command, "flush_evt")==0) {
    server = open_connection_with_list(hostname, port);
    errorCode = buffer_flushevt(server, &(plhs[0]), &(prhs[1]));
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
  
  else if (strcasecmp(command, "list_connections")==0) {
	dump_hps_list();
  }
  
  else if (strcasecmp(command, "lookup_connection")==0) {
	plhs[0] = mxCreateDoubleScalar(lookup_hps_item(hostname, port));
  }
  
  else if (strcasecmp(command, "close_connection")==0) {
	int sock = lookup_hps_item(hostname, port);
	/* Note that we ignore request to close the dma "connection" */
	if (sock > 0) {
		close_connection(sock);
		remove_hps_item(sock);
	}
	dump_hps_list();
  }
  
  else {
    mexErrMsgTxt ("unknown command");
  }
  
  /* These are either allocated using mxMalloc'ed or mxArrayToString.
	 Deallocating them manually is not strictly necessary, but cleaner (says Mathworks)
  */
  if (argument) mxFree(argument); 
  if (command)  mxFree(command);
  if (hostname) mxFree(hostname); 
  
  if (errorCode > 0) {
    char msg[512];
    /* valid connection, but invalid request */
    sprintf(msg, "ERROR: the buffer returned an error (%d)\n", errorCode);
    mexErrMsgTxt(msg);
  } 
  else if (errorCode == -3) {
    /* valid connection, but invalid request, we'll close the socket and remove the list item */
    printf("Trying to close socket %i\n", server);
	close_connection(server);
	remove_hps_item(server);
    mexErrMsgTxt("ERROR: tcp connection to buffer failed.");
  }
  
  return;
}

