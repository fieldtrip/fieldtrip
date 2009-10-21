/*
 * Copyright (C) 2008, Robert Oostenveld & Christian Hesse
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 * $Log: tcpserver.c,v $
 * Revision 1.15  2009/01/23 19:47:50  roboos
 * changed verbosity
 *
 * Revision 1.14  2009/01/23 08:26:44  roboos
 * fixed a serious bug that caused a lot of memory to leak (in fact all packets that were sent over the socket would eventually leack away), both on the client and server side
 *
 * Revision 1.13  2009/01/21 20:54:35  roboos
 * added some debugging code to keep track of number of sockets and threads (both seem ok)
 * cleaned up the debugging: more consistent use of verbose flag and function name in fprintf
 *
 * Revision 1.12  2009/01/21 08:47:22  roboos
 * fixed bug in thread counting
 *
 * Revision 1.11  2009/01/21 08:33:01  roboos
 * added some code for debugging the leaking of sockets and threads
 *
 * Revision 1.10  2008/05/22 13:40:06  roboos
 * moved declarations to top
 *
 * Revision 1.9  2008/05/22 10:00:31  roboos
 * some small changes related to compatibility with Borland, thanks to Jurgen
 *
 * Revision 1.8  2008/05/21 19:45:20  roboos
 * only some changes in perror/fprintf output
 *
 * Revision 1.7  2008/03/26 14:36:12  thohar
 * made new code to unblock sockets portable
 *
 * Revision 1.6  2008/03/23 16:59:21  roboos
 * added a cleanup label -> goto there in case something is wrong: the thread cancelation is handled there
 *
 * Revision 1.5  2008/03/23 16:43:56  roboos
 * implemented thread cancelation: first put socket in non-blocking
 * mode, when accept returns an error do pthread_testcancel and pause
 * 1 milisecond, then continue with next accept. The socket that is
 * opened by the client and passed to tcpsocket is put back into
 * blocking mode, otherwise tcpclient/bufread do not detect correctly
 * that it gets closed at a certain point.
 *
 * Revision 1.4  2008/03/13 13:38:31  roboos
 * added cancel options for thread
 *
 * Revision 1.3  2008/03/13 12:34:22  thohar
 * added headers for win32 build
 * added functions to start sockets on win32
 *
 * Revision 1.2  2008/03/10 10:03:47  roboos
 * use *host_t as input argument
 *
 * Revision 1.1  2008/03/08 09:41:02  roboos
 * renamed buffer_thread function to tcpserver, idem for the c-file
 * renamed socket_thread function to tcpsocket, idem for the c-file
 *
 * Revision 1.7  2008/02/27 10:13:15  roboos
 * disabled fprintf statement
 *
 * Revision 1.6  2008/02/20 13:36:52  roboos
 * added counter for incoming connections
 *
 * Revision 1.5  2008/02/20 07:10:25  roboos
 * tried out some low-level tcp specific details
 *
 * Revision 1.4  2008/02/19 10:24:42  roboos
 * removed some old documentation from the comments
 *
 * Revision 1.3  2008/02/19 10:22:55  roboos
 * added consistent copyrights and log message to each of the files
 *
 * Revision 1.2  2008/02/18 12:13:46  roboos
 * moved executable from buffer to demo
 * fixed bugs in sinewave and socket for events
 * stripped down the eventdef_t fields
 * many small changes
 *
 * Revision 1.1  2008/02/18 10:05:25  roboos
 * restructured the directory layout, copied all code into src, added directory for external code
 *
 * Revision 1.8  2008/02/13 13:07:56  roboos
 * fixed numerous bugs
 * implemented append function, which after all does not seem to be neccessary (since the problem is in the TCP buffer size)
 *
 * Revision 1.7  2008/02/11 21:12:22  roboos
 * added check for each datatype wordsize
 * deleted global numevents numsamples
 *
 * Revision 1.6  2008/02/10 20:26:25  roboos
 * only some small changes
 *
 * Revision 1.5  2008/02/10 10:48:43  roboos
 * removed init function, changed data from double into void, reindented the complete file
 *
 * Revision 1.4  2008/02/09 16:45:47  roboos
 * moved declarations that are shared between threads to seperate header file
 * moved code that deals with incoming socket connection to seperate file (socket.c)
 *
 * Revision 1.3  2008/02/08 10:32:25  roboos
 * replaced the time and date code by some real buffer handling
 * this version compiles, but may not work (untested)
 *
 * Revision 1.2  2008/02/05 20:32:04  roboos
 * First woring implementation, based on daytimed example code from
 * http://www.freebsd.org/doc/en_US.ISO8859-1/books/developers-handbook/sockets-essential-functions.html
 * and on threads example code from https://computing.llnl.gov/tutorials/pthreads/fork_vs_thread.txt
 *
 */
 
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* these are for changing the socket to non-blocking mode */
#include <fcntl.h>
#include <errno.h>

#include "message.h"
#include "buffer.h"
#include "socket_includes.h"
#include "unix_includes.h"

#define ACCEPTSLEEP 1000

//extern int errno;

//pthread_attr_t attr; // this one would be passed to the thread

extern pthread_mutex_t mutexthreadcount;
extern int threadcount;

/***********************************************************************
 * this thread listens to incoming TCP connections
 * if a connection is made by a client, it starts the tcpsocket function
 ***********************************************************************/
void *tcpserver(void *arg) {
    int verbose = 0;

	/* these variables are for the socket */
	struct sockaddr_in sa;
	int s, c;
    unsigned int b;
	int optval;
	struct timeval timeout;
	int oldcancelstate, oldcanceltype;

#ifdef WIN32
    unsigned long enable = 0;
	WSADATA wsa;
#endif

	/* these variables are for the threading */
	int rc;
	pthread_t tid;

	/* this determines the port on which the server will listen */
	host_t *host = (host_t *)arg;
	if (!arg)
		exit(1);

	if (verbose>0) fprintf(stderr, "tcpserver: host.name =  %s\n", host->name);
	if (verbose>0) fprintf(stderr, "tcpserver: host.port =  %d\n", host->port);

	/* the thread is allowed to be canceled at almost any moment */
	pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, &oldcancelstate);
	pthread_setcanceltype(PTHREAD_CANCEL_DEFERRED, &oldcanceltype);
	pthread_cleanup_push(cleanup_socket, (void *)&s);

#ifdef WIN32
 	if(WSAStartup(MAKEWORD(1, 1), &wsa))
	{
		if (verbose>0) fprintf(stderr, "tcpserver: cannot start sockets\n");
	}
#endif

	/* setup socket */
	if ((s = socket(PF_INET, SOCK_STREAM, 0)) < 0) {
		perror("tcpserver socket");
		goto cleanup;
	}

	/* place the socket in non-blocking mode, required to do thread cancelation */
#ifdef WIN32
	enable = 0;
	ioctlsocket(s, FIONBIO, &enable);
#else
	optval = fcntl(s, F_GETFL, NULL);
	optval = optval | O_NONBLOCK;
	if (fcntl(s, F_SETFL, optval)<0) {
		perror("tcpserver fcntl");
		goto cleanup;
	}
#endif

	/* change the receive timeout */
	/*
	   timeout.tv_sec  = 1;
	   timeout.tv_usec = 1;
	   if (setsockopt(s, SOL_SOCKET, SO_RCVTIMEO, &timeout, sizeof(optval)) < 0) {
	   perror("tcpserver setsockopt");
	   if (verbose>0) fprintf(stderr, "tcpserver: errno = %d\n", errno);
	   goto cleanup;
	   }
	 */

	/* prevend "bind: address already in use" */
	optval = 1;
	if (setsockopt(s, SOL_SOCKET, SO_REUSEADDR, (const char*)&optval, sizeof(optval)) < 0) {
		perror("tcpserver setsockopt");
		goto cleanup;
	}

	bzero(&sa, sizeof sa);

	sa.sin_family = AF_INET;
	sa.sin_port   = htons(host->port);

	if (INADDR_ANY)
		sa.sin_addr.s_addr = htonl(INADDR_ANY);

	if (bind(s, (struct sockaddr *)&sa, sizeof sa) < 0) {
		perror("tcpserver bind");
		goto cleanup;
	}

	if (listen(s, BACKLOG)<0) {
		perror("tcpserver listen");
		goto cleanup;
	}

	for (;;) {
		/*
		 * If no pending connections are present on the queue, and the socket
		 * is not marked as non-blocking, accept() blocks the caller until a
		 * connection is present.  If the socket is marked non-blocking and
		 * no pending connections are present on the queue, accept() returns
		 * an error as described below.
		 */

		b = sizeof sa;
		c = accept(s, (struct sockaddr *)&sa, &b);

	if (c<0) {
#ifdef WIN32
			if(errno == 0) {
				pthread_testcancel();
				usleep(ACCEPTSLEEP);
			}
			else {
				perror("tcpserver accept");
				goto cleanup;
			}
#else
			if (errno==EWOULDBLOCK) {
				pthread_testcancel();
				usleep(ACCEPTSLEEP);
			}
			else {
				perror("tcpserver accept");
				goto cleanup;
			}
#endif
		}

		else {
			if (verbose>0) fprintf(stderr, "tcpserver: opened connection to client on socket %d\n", c);

			/* set larger buffer */
			optval = SO_RCVBUF_SIZE;
			if (setsockopt(c, SOL_SOCKET, SO_RCVBUF, (const char*)&optval, sizeof(optval)) < 0) {
				perror("tcpserver setsockopt");
				goto cleanup;
			}

			/* set larger buffer */
			optval = SO_SNDBUF_SIZE;
			if (setsockopt(c, SOL_SOCKET, SO_SNDBUF, (const char*)&optval, sizeof(optval)) < 0) {
				perror("tcpserver setsockopt");
				goto cleanup;
			}

			/* place the socket back in blocking mode, this is needed for tcpsocket  */
#ifdef WIN32
			enable = 0;
			ioctlsocket(s, FIONBIO, &enable);
#else
			optval = fcntl(c, F_GETFL, NULL);
			optval = optval & (!O_NONBLOCK);
			if (fcntl(c, F_SETFL, optval)<0) {
				perror("tcpserver fcntl");
				goto cleanup;
			}
#endif

			/* disable the Nagle buffering algorithm */
			/*
			   optval = 1;
			   if (setsockopt(c, IPPROTO_TCP, TCP_NODELAY, &optval, sizeof(optval)) < 0) {
			   perror("tcpserver setsockopt");
			   goto cleanup;
			   }
			 */

			/* deal with the incoming connection on the TCP socket in a seperate thread */
			/* rc = pthread_create(&tid, &attr, tcpsocket, (void *)c); */
			rc = pthread_create(&tid, NULL, tcpsocket, (void *)c);

			if (rc) {
				if (verbose>0) fprintf(stderr, "tcpserver: return code from pthread_create() is %d\n", rc);
				goto cleanup;
			}
			else {
                /* this is for debugging */
                pthread_mutex_lock(&mutexthreadcount);
                threadcount++;
                pthread_mutex_unlock(&mutexthreadcount);
			    if (verbose>0) fprintf(stderr, "tcpserver: c = %d, threadcount = %d\n", c, threadcount);
				pthread_detach(tid);
			}
		}
	}

cleanup:
	/* from now on it is safe to cancel the thread */
	pthread_setcancelstate(oldcancelstate, NULL);
	pthread_setcanceltype(oldcanceltype, NULL);

	pthread_cleanup_pop(1);  // socket
	pthread_exit(NULL);
	return NULL;
}

