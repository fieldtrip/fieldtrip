/*
 * Copyright (C) 2008-2010, Robert Oostenveld
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

/* these are for changing the socket to non-blocking mode */
#include <fcntl.h>
#include <errno.h>

#include "peer.h"
#include "extern.h"
#include "platform_includes.h"

typedef struct {
		int *fd;
} threadlocal_t;

void cleanup_tcpserver(void *arg) {
		threadlocal_t *threadlocal;
		threadlocal = (threadlocal_t *)arg;

        DEBUG(LOG_DEBUG, "cleanup_tcpserver()");

		if (tcpserverStatus==0)
				return;

		if (threadlocal && (*threadlocal->fd)>0) {
				closesocket(*threadlocal->fd);
				*threadlocal->fd = 0;
		}

		pthread_mutex_lock(&mutexhost);
		host->port = 0;
		pthread_mutex_unlock(&mutexhost);

		pthread_mutex_lock(&mutexstatus);
		tcpserverStatus = 0;
		pthread_mutex_unlock(&mutexstatus);

		pthread_mutex_lock(&mutexthreadcount);
		threadcount--;
		pthread_mutex_unlock(&mutexthreadcount);
}

/***********************************************************************
 * this thread listens to incoming TCP connections
 * if a connection is made by a client, it starts the tcpsocket function
 ***********************************************************************/
void *tcpserver(void *arg) {
		int c, retry;
		int fd = 0;

		/* these variables are for the socket */
		struct sockaddr_in sa;
		unsigned int b;
		int optval;

		/* these variables are for the threading */
		int rc;
		pthread_t tid;

#ifdef WIN32
		unsigned long enable = 0;
		WSADATA wsa;
#endif

		threadlocal_t threadlocal;
		threadlocal.fd = &fd;

		/* this is for debugging */
		pthread_mutex_lock(&mutexthreadcount);
		threadcount++;
		pthread_mutex_unlock(&mutexthreadcount);

		pthread_cleanup_push(cleanup_tcpserver, &threadlocal);

		/* the status contains the thread id when running, or zero when not running */
		pthread_mutex_lock(&mutexstatus);
		if (tcpserverStatus==0) {
				tcpserverStatus = 1;
				/* signal that this thread has started */
				pthread_cond_signal(&condstatus);
				pthread_mutex_unlock(&mutexstatus);
		}
		else {
				pthread_mutex_unlock(&mutexstatus);
				goto cleanup;
		}

#ifdef WIN32
		if(WSAStartup(MAKEWORD(1, 1), &wsa))
		{
				DEBUG(LOG_ERR, "tcpserver: cannot start sockets");
				/* FIXME should this be handled more explicitely? */
		}
#endif

		/* setup socket */
		if ((fd = socket(PF_INET, SOCK_STREAM, 0)) < 0) {
				perror("tcpserver socket");
				DEBUG(LOG_ERR, "error: tcpserver socket");
				goto cleanup;
		}

		/* place the socket in non-blocking mode, required to do thread cancelation */
#ifdef WIN32
		enable = 0;
		ioctlsocket(fd, FIONBIO, &enable);
#else
		optval = fcntl(fd, F_GETFL, NULL);
		optval = optval | O_NONBLOCK;
		if (fcntl(fd, F_SETFL, optval)<0) {
				perror("tcpserver fcntl");
				DEBUG(LOG_ERR, "error: tcpserver fcntl");
				goto cleanup;
		}
#endif

		/* change the receive timeout */
		/*
		   timeout.tv_sec  = 1;
		   timeout.tv_usec = 1;
		   if (setsockopt(fd, SOL_SOCKET, SO_RCVTIMEO, &timeout, sizeof(optval)) < 0) {
		   perror("tcpserver setsockopt");
		   DEBUG(LOG_ERR, "error: tcpserver setsockopt");
		   goto cleanup;
		   }
		 */

		/* prevend "bind: Address already in use" */
		/*
		   optval = 1;
		   if (setsockopt(fd, SOL_SOCKET, SO_REUSEADDR, (const char*)&optval, sizeof(optval)) < 0) {
		   perror("tcpserver setsockopt");
		   DEBUG(LOG_ERR, "error: tcpserver setsockopt");
		   goto cleanup;
		   }
		 */

		/* set larger buffer */
		optval = SO_RCVBUF_SIZE;
		if (setsockopt(fd, SOL_SOCKET, SO_RCVBUF, (const char*)&optval, sizeof(optval)) < 0) {
				perror("tcpserver setsockopt");
				DEBUG(LOG_ERR, "error: tcpserver setsockopt");
				goto cleanup;
		}

		/* set larger buffer */
		optval = SO_SNDBUF_SIZE;
		if (setsockopt(fd, SOL_SOCKET, SO_SNDBUF, (const char*)&optval, sizeof(optval)) < 0) {
				perror("tcpserver setsockopt");
				DEBUG(LOG_ERR, "error: tcpserver setsockopt");
				goto cleanup;
		}

		/* disable the Nagle buffering algorithm */
		/*
		   optval = 1;
		   if (setsockopt(fd, IPPROTO_TCP, TCP_NODELAY, &optval, sizeof(optval)) < 0) {
		   perror("tcpserver setsockopt");
		   DEBUG(LOG_ERR, "error: tcpserver setsockopt");
		   goto cleanup;
		   }
		 */

		bzero(&sa, sizeof sa);
		b = sizeof sa;

		pthread_mutex_lock(&mutexhost);
		if (host->port==0)
		  host->port = DEFAULT_PORT;
		sa.sin_family = AF_INET;
		sa.sin_port   = htons(host->port);

		if (INADDR_ANY)
				sa.sin_addr.s_addr = htonl(INADDR_ANY);

		retry = 1000;
		while (retry>0) {
				if (bind(fd, (struct sockaddr *)&sa, sizeof sa)<0) {
						/* increment the port number and try again */
						host->port++;
						sa.sin_port = htons(host->port);
						retry--;
				}
				else {
						/* this signals that it was successful */
						retry = -1;
				}
		}
		DEBUG(LOG_DEBUG, "tcpserver: started on port %d, id = %\n", host->port, host->id);
		pthread_mutex_unlock(&mutexhost);

		if (retry==0) {
				/* it failed on mutliple attempts, give up */
				perror("tcpserver bind");
				DEBUG(LOG_ERR, "error: tcpserver bind");
				goto cleanup;
		}

		if (listen(fd, BACKLOG)<0) {
				perror("tcpserver listen");
				DEBUG(LOG_ERR, "error: tcpserver listen");
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

				c = accept(fd, (struct sockaddr *)&sa, &b);

				if (c<0) {
#ifdef WIN32
						if(errno==0) {
								pthread_testcancel();
								threadsleep(ACCEPTSLEEP);
						}
						else {
								perror("tcpserver accept");
								DEBUG(LOG_ERR, "error: tcpserver accept");
								goto cleanup;
						}
#else
						if (errno==EWOULDBLOCK) {
								pthread_testcancel();
								threadsleep(ACCEPTSLEEP);
						}
						else {
								perror("tcpserver accept");
								perror("tcpserver accept");
								goto cleanup;
						}
#endif
				}

				else {
						DEBUG(LOG_DEBUG, "tcpserver: opened connection to client on socket %d", c);

						/* place the socket back in blocking mode, this is needed for tcpsocket  */
#ifdef WIN32
						enable = 0;
						ioctlsocket(c, FIONBIO, &enable);
#else
						optval = fcntl(c, F_GETFL, NULL);
						optval = optval & (!O_NONBLOCK);
						if (fcntl(c, F_SETFL, optval)<0) {
								perror("tcpserver fcntl");
								perror("tcpserver fcntl");
								goto cleanup;
						}
#endif

						/* deal with the incoming connection on the TCP socket in a separate thread */
						/* rc = pthread_create(&tid, &attr, tcpsocket, (void *)c); */
						rc = pthread_create(&tid, NULL, tcpsocket, (void *)c);

						if (rc) {
								/* the code should never arrive here */
								DEBUG(LOG_ERR, "tcpserver: return code from pthread_create() is %d", rc);
								goto cleanup;
						}

						DEBUG(LOG_DEBUG, "tcpserver: c = %d, threadcount = %d", c, threadcount);
						pthread_detach(tid);
				}
		}

cleanup:
		printf(""); /* otherwise the pthread_cleanup_pop won't compile */
		pthread_cleanup_pop(1);
		return NULL;
}

