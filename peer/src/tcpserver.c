/* these are for changing the socket to non-blocking mode */
#include <fcntl.h>
#include <errno.h>

#include "peer.h"
#include "extern.h"
#include "platform_includes.h"

typedef struct {
		int fd;
} threadlocal_t;

void cleanup_tcpserver(void *arg) {
		threadlocal_t *threadlocal;
		threadlocal = (threadlocal_t *)arg;
		if (threadlocal && threadlocal->fd>0) {
				close(threadlocal->fd);
				threadlocal->fd = -1;
		}

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
		int verbose = 0;
		int c, fd, retry;

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
		threadlocal.fd = -1;

		/* this is for debugging */
		pthread_mutex_lock(&mutexthreadcount);
		threadcount++;
		pthread_mutex_unlock(&mutexthreadcount);

		pthread_cleanup_push(cleanup_tcpserver, &threadlocal);

		/* the status contains the thread id when running, or zero when not running */
		pthread_mutex_lock(&mutexstatus);
		if (tcpserverStatus==0) {
				tcpserverStatus = 1;
				pthread_mutex_unlock(&mutexstatus);
		}
		else {
				pthread_mutex_unlock(&mutexstatus);
				goto cleanup;
		}

#ifdef WIN32
		if(WSAStartup(MAKEWORD(1, 1), &wsa))
		{
				if (verbose>0) fprintf(stderr, "tcpserver: cannot start sockets\n");
		}
#endif

		if (verbose>1)
				fprintf(stderr, "tcpserver: port = %d\n", host->port);

		/* setup socket */
		if ((fd = socket(PF_INET, SOCK_STREAM, 0)) < 0) {
				perror("tcpserver socket");
				goto cleanup;
		}

		/* this will be closed at cleanup */
		threadlocal.fd = fd;

		/* place the socket in non-blocking mode, required to do thread cancelation */
#ifdef WIN32
		enable = 0;
		ioctlsocket(fd, FIONBIO, &enable);
#else
		optval = fcntl(fd, F_GETFL, NULL);
		optval = optval | O_NONBLOCK;
		if (fcntl(fd, F_SETFL, optval)<0) {
				perror("tcpserver fcntl");
				goto cleanup;
		}
#endif

		/* change the receive timeout */
		/*
		   timeout.tv_sec  = 1;
		   timeout.tv_usec = 1;
		   if (setsockopt(fd, SOL_SOCKET, SO_RCVTIMEO, &timeout, sizeof(optval)) < 0) {
		   perror("tcpserver setsockopt");
		   if (verbose>0) fprintf(stderr, "tcpserver: errno = %d\n", errno);
		   goto cleanup;
		   }
		 */

		/* prevend "bind: Address already in use" */
		optval = 1;
		if (setsockopt(fd, SOL_SOCKET, SO_REUSEADDR, (const char*)&optval, sizeof(optval)) < 0) {
				perror("tcpserver setsockopt");
				goto cleanup;
		}

		/* set larger buffer */
		optval = SO_RCVBUF_SIZE;
		if (setsockopt(fd, SOL_SOCKET, SO_RCVBUF, (const char*)&optval, sizeof(optval)) < 0) {
				perror("tcpserver setsockopt");
				goto cleanup;
		}

		/* set larger buffer */
		optval = SO_SNDBUF_SIZE;
		if (setsockopt(fd, SOL_SOCKET, SO_SNDBUF, (const char*)&optval, sizeof(optval)) < 0) {
				perror("tcpserver setsockopt");
				goto cleanup;
		}

		/* disable the Nagle buffering algorithm */
		/*
		   optval = 1;
		   if (setsockopt(fd, IPPROTO_TCP, TCP_NODELAY, &optval, sizeof(optval)) < 0) {
		   perror("tcpserver setsockopt");
		   goto cleanup;
		   }
		 */

		bzero(&sa, sizeof sa);
		b = sizeof sa;

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
		if (retry==0) {
				/* it failed on mutliple attempts, give up */
				perror("tcpserver bind");
				goto cleanup;
		}

		if (listen(fd, BACKLOG)<0) {
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

				c = accept(fd, (struct sockaddr *)&sa, &b);

				if (verbose>1)
						fprintf(stderr, "tcpserver: c = %d\n", c);

				if (c<0) {
#ifdef WIN32
						if(errno==0) {
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

						/* place the socket back in blocking mode, this is needed for tcpsocket  */
#ifdef WIN32
						enable = 0;
						ioctlsocket(c, FIONBIO, &enable);
#else
						optval = fcntl(c, F_GETFL, NULL);
						optval = optval & (!O_NONBLOCK);
						if (fcntl(c, F_SETFL, optval)<0) {
								perror("tcpserver fcntl");
								goto cleanup;
						}
#endif

						/* deal with the incoming connection on the TCP socket in a seperate thread */
						/* rc = pthread_create(&tid, &attr, tcpsocket, (void *)c); */
						rc = pthread_create(&tid, NULL, tcpsocket, (void *)c);

						if (rc) {
								/* the code should never arrive here */
								if (verbose>0) fprintf(stderr, "tcpserver: return code from pthread_create() is %d\n", rc);
								goto cleanup;
						}

						if (verbose>0) fprintf(stderr, "tcpserver: c = %d, threadcount = %d\n", c, threadcount);
						pthread_detach(tid);
				}
		}

cleanup:
		printf(""); /* otherwise the pthread_cleanup_pop won't compile */
		pthread_cleanup_pop(1);
		return NULL;
}

