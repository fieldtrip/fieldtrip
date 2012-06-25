/*
 * Copyright (C) 2008, Robert Oostenveld
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 */

#include <string.h>

#include "mex.h"
#include "matrix.h"
#include "buffer.h"
#include <pthread.h>
#include "extern.h"

#define DEFAULT_HOST "localhost"
#define DEFAULT_PORT 1972


#ifdef WIN32 
#  ifndef COMPILER_LCC
#    define sleep(x)   Sleep((x)*1000)
#  endif
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

/* This keeps track of the running thread */
pthread_t tcpserverThread;

/* a pointer to this is passed on to the thread, we only need it once
   and therefore use a static variable */
host_t host_for_server; 

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
int buffer_waitdat(int, mxArray **, const mxArray **);

/* this function is called upon unloading of the mex-file */
void exitFun(void) {
	int verbose = 1;
	int rc;
  
	if (verbose) {
		printf("Entering exitFun() routine\n");
	}

    /* tell the tcpserver thread to stop */
    pthread_mutex_lock(&mutexstatus);
    if (tcpserverStatus) {
            pthread_mutex_unlock(&mutexstatus);
            mexPrintf("requesting cancelation of tcpserver thread\n");
            pthread_cancel(tcpserverThread);
            pthread_join(tcpserverThread, NULL);
    }
    else {
            pthread_mutex_unlock(&mutexstatus);
    }

	/* free the memory that is used for the header, data and events */
	free_event();
	free_data();
	free_header();

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
	
	mexAtExit(exitFun);   /* register cleanup routine so the list get's properly destroyed */
	return 1;
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
	if (sock == 0) {
		/* dma connection to internal buffer */
		if (verbose > 0) 
			printf("Found (dma) connection to internal buffer\n");
		return sock;
	}
	if (sock > 0) {
		/* First check whether that socket can be read from, indicating a 
			closed connection from the remote side.
		*/
		fd_set readSet;
		struct timeval tv;
		int sel;
		
		tv.tv_sec = 0;
		tv.tv_usec = 0;
		
		FD_ZERO(&readSet);
		FD_SET(sock, &readSet);
		
		sel = select(sock + 1, &readSet, NULL, NULL, &tv);
		
		if (sel == 1) {
			mexWarnMsgTxt("Existing connection seems to have been closed from remote side - trying to re-open...\n");
			remove_hps_item(sock);
			/* go on further down with opening new connection */
		} else {
			/* reading not possible means that connection is still open */
			if (verbose > 0) 
				printf("Found connection on %i\n",sock);
			return sock;
		}
	}
	
	if (verbose>0) 
		printf("Trying to open new connection...\n");

	if (port == 0) {
		sock = open_unix_connection(hostname);
	} else {
		sock = open_connection(hostname, port);
	}
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
		
      if (tcpserverStatus)
        mexErrMsgTxt("thread is already running");
      mexPrintf("In main: spawning tcpserver thread\n");
	  
      strncpy(host_for_server.name, hostname, 256);  /* the 256 is a limit inside the host_t structure */
      host_for_server.port = port;
	  
      rc = pthread_create(&tcpserverThread, NULL, tcpserver, (void *)&host_for_server);
      if (rc)
        mexErrMsgTxt("problem with return code from pthread_create()");
		
	  /* put this as a "connection" with sock=0 into the list */
	  add_hps_item(hostname, port, 0);
    }
    else if (strcasecmp(argument, "exit")==0) {
      if (!tcpserverStatus)
        mexErrMsgTxt("thread is not running");
      mexPrintf("In main: requesting cancelation of tcpserver thread\n");
      rc = pthread_cancel(tcpserverThread);
      if (rc)
        mexErrMsgTxt("problem with return code from pthread_cancel()");
	  /* remove sock=0 connection from host/port/socket list */
	  remove_hps_item(0);
    }
	else if (strcasecmp(argument, "status")==0) {
		plhs[0] = mxCreateDoubleScalar(tcpserverStatus);
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
  
  else if (strcasecmp(command, "wait_dat")==0) {
    server = open_connection_with_list(hostname, port);
    errorCode = buffer_waitdat(server, &(plhs[0]), &(prhs[1]));
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

