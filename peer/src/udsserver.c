/*
 * Copyright (C) 2010, Robert Oostenveld
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

#ifdef WIN32
/* this is not yet implemented on windows */
#else
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/stat.h> /* for chmod */
#include <sys/un.h>
#include <pthread.h>
#endif

#define SOCK_PATH "/tmp/peer"

#include "peer.h"
#include "extern.h"
#include "platform_includes.h"

typedef struct {
		int fd;
} threadlocal_t;

void cleanup_udsserver(void *arg) {
		threadlocal_t *threadlocal;
		threadlocal = (threadlocal_t *)arg;
		if (threadlocal && threadlocal->fd>0) {
				closesocket(threadlocal->fd);
				threadlocal->fd = -1;
		}

		pthread_mutex_lock(&mutexhost);
		/* prevent the invalid unix domain socket from being announced */
		bzero(host->socket, STRLEN);
		pthread_mutex_unlock(&mutexhost);

		pthread_mutex_lock(&mutexstatus);
		udsserverStatus = 0;
		pthread_mutex_unlock(&mutexstatus);

		pthread_mutex_lock(&mutexthreadcount);
		threadcount--;
		pthread_mutex_unlock(&mutexthreadcount);
}

/***********************************************************************
 * this thread listens to incoming UDS connections
 * if a connection is made by a client, it starts the udssocket function
 ***********************************************************************/
void *udsserver(void *arg) {
#ifdef WIN32
		/* this is not yet implemented on windows */
#else
		int verbose = 0;
		int c, fd, retry;

		/* these variables are for the socket */
		struct sockaddr_un local, remote;
		int len;
		int optval;

		/* these variables are for the threading */
		int rc;
		pthread_t tid;

		threadlocal_t threadlocal;
		threadlocal.fd = -1;

		/* this is for debugging */
		pthread_mutex_lock(&mutexthreadcount);
		threadcount++;
		pthread_mutex_unlock(&mutexthreadcount);

		pthread_cleanup_push(cleanup_udsserver, &threadlocal);

		/* the status contains the thread id when running, or zero when not running */
		pthread_mutex_lock(&mutexstatus);
		if (udsserverStatus==0) {
				udsserverStatus = 1;
				pthread_mutex_unlock(&mutexstatus);
		}
		else {
				pthread_mutex_unlock(&mutexstatus);
				goto cleanup;
		}

		if (verbose>1)
				fprintf(stderr, "udsserver: port = %d\n", host->port);

		/* setup socket */
		if ((fd = socket(AF_UNIX, SOCK_STREAM, 0)) < 0) {
				perror("udsserver socket");
				goto cleanup;
		}

		/* this will be closed at cleanup */
		threadlocal.fd = fd;

		bzero(&local, sizeof local);
		bzero(&remote, sizeof remote);

		pthread_mutex_lock(&mutexhost);
		local.sun_family = AF_UNIX;
		sprintf(local.sun_path, "%s-%s.%d", SOCK_PATH, host->user, getpid());  /* for example "/tmp/peer-roboos.1436" */
		strncpy(host->socket, local.sun_path, STRLEN);
		pthread_mutex_unlock(&mutexhost);

		/* remove the file if it already exists */
		unlink(local.sun_path);

		len = strlen(local.sun_path) + sizeof(local.sun_family);
#ifdef USE_ABSTRACT_UDS_NAMES
		/* abstract unix domain socket namea do not show up in the file system */
		local.sun_path[0] = 0;
#endif
		if (bind(fd, (struct sockaddr *)&local, len) == -1) {
				perror("udsserver bind");
				goto cleanup;
		}

#ifndef USE_ABSTRACT_UDS_NAMES
		/* this is required to allow other users to read and write to the socket through the file system*/
		if (chmod(local.sun_path, 0777)!=0)
				perror('chmod');
#endif

		if (listen(fd, BACKLOG)<0) {
				perror("udsserver listen");
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

				len = sizeof remote;
				c = accept(fd, (struct sockaddr *)&remote, &len);

				if (verbose>1)
						fprintf(stderr, "udsserver: c = %d\n", c);

				if (c<0) {
						if (errno==EWOULDBLOCK) {
								pthread_testcancel();
								usleep(ACCEPTSLEEP);
						}
						else {
								perror("udsserver accept");
								goto cleanup;
						}
				}

				else {
						if (verbose>0) fprintf(stderr, "udsserver: opened connection to client on socket %d\n", c);

						/* deal with the incoming connection on the UDS socket in a seperate thread */
						/* rc = pthread_create(&tid, &attr, udssocket, (void *)c); */
						rc = pthread_create(&tid, NULL, tcpsocket, (void *)c);

						if (rc) {
								/* the code should never arrive here */
								if (verbose>0) fprintf(stderr, "udsserver: return code from pthread_create() is %d\n", rc);
								goto cleanup;
						}

						if (verbose>0) fprintf(stderr, "udsserver: c = %d, threadcount = %d\n", c, threadcount);
						pthread_detach(tid);
				}
		}

cleanup:
		printf(""); /* otherwise the pthread_cleanup_pop won't compile */
		pthread_cleanup_pop(1);
		return NULL;
#endif
}

