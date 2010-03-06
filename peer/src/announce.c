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

#include <sys/types.h>
#include <time.h>
#include <string.h>
#include <stdio.h>
#include <pthread.h>

#include "peer.h"
#include "extern.h"
#include "platform_includes.h"

typedef struct {
		void *message;
		int fd;
} threadlocal_t;

void cleanup_announce(void *arg) {
		threadlocal_t *threadlocal;
		threadlocal = (threadlocal_t *)arg;
		if (threadlocal && threadlocal->message) {
				FREE(threadlocal->message);
		}
		if (threadlocal && threadlocal->fd>0) {
				close(threadlocal->fd);
				threadlocal->fd = -1;
		}

		pthread_mutex_lock(&mutexstatus);
		announceStatus = 0;
		pthread_mutex_unlock(&mutexstatus);

		pthread_mutex_lock(&mutexthreadcount);
		threadcount--;
		pthread_mutex_unlock(&mutexthreadcount);
}

void *announce(void *arg) {
		int fd = 0;
		int verbose = 0;
		struct sockaddr_in multicast;
		hostdef_t *message = NULL;
		unsigned char ttl = 3;
		unsigned char one = 1;

		threadlocal_t threadlocal;
		threadlocal.message = NULL;
		threadlocal.fd = -1;

		/* this is for debugging */
		pthread_mutex_lock(&mutexthreadcount);
		threadcount++;
		pthread_mutex_unlock(&mutexthreadcount);

		pthread_cleanup_push(cleanup_announce, &threadlocal);

		/* the status contains the thread id when running, or zero when not running */
		pthread_mutex_lock(&mutexstatus);
		if (announceStatus==0) {
				announceStatus = 1;
				pthread_mutex_unlock(&mutexstatus);
		}
		else {
				pthread_mutex_unlock(&mutexstatus);
				goto cleanup;
		}

		/* allocate memory for the message */
		if ((message = (hostdef_t *)malloc(sizeof(hostdef_t))) == NULL) {
				perror("announce malloc");
				goto cleanup;
		}

		/* this will be deallocated at cleanup */
		threadlocal.message = message;

		if (verbose>1)
				fprintf(stderr, "announce: threadcount = %d\n", threadcount);

		if (verbose>0) {
				pthread_mutex_lock(&mutexhost);
				fprintf(stderr, "announce: host.name = %s\n", host->name);
				fprintf(stderr, "announce: host.port = %d\n", host->port);
				fprintf(stderr, "announce: host.id   = %d\n", host->id);
				pthread_mutex_unlock(&mutexhost);
		}

		/* create what looks like an ordinary UDP socket */
		if ((fd=socket(AF_INET,SOCK_DGRAM,IPPROTO_UDP)) < 0) {
				perror("announce socket");
				goto cleanup;
		}

		/* this will be closed at cleanup */
		threadlocal.fd = fd;

		/* set up destination address for multicasting */
		memset(&multicast,0,sizeof(multicast));
		multicast.sin_family      = AF_INET;
		multicast.sin_addr.s_addr = inet_addr(ANNOUNCE_GROUP);
		multicast.sin_port        = htons(ANNOUNCE_PORT);

		/* now just sendto() our destination */
		while (1) {

				pthread_mutex_lock(&mutexhost);

				/* the host details can change over time */
				memcpy(message, host, sizeof(hostdef_t));

				if (strncasecmp(host->name, "localhost", STRLEN)==0)
						ttl = 0;
				else
						ttl = 1;

				/*  the TTL (time to live/hop count) for the send can be
				 *  ----------------------------------------------------------------------
				 *  0      Restricted to the same host. Won't be output by any interface.
				 *  1      Restricted to the same subnet. Won't be forwarded by a router.
				 *  <32    Restricted to the same site, organization or department.
				 *  <64    Restricted to the same region.
				 *  <128   Restricted to the same continent.
				 *  <255   Unrestricted in scope. Global.
				 */

				if ((setsockopt(fd, IPPROTO_IP, IP_MULTICAST_TTL, (void*)&ttl, sizeof(ttl))) < 0)
						perror("setsockopt() failed");

				pthread_mutex_unlock(&mutexhost);

				if (verbose>1)
						fprintf(stderr, "announce\n");

				/* note that this is a thread cancelation point */
				if (sendto(fd,message,sizeof(hostdef_t),0,(struct sockaddr *) &multicast,sizeof(multicast)) < 0) {
						perror("announce sendto");
						goto cleanup;
				}

				/* note that this is a thread cancelation point */
				pthread_testcancel();
				usleep(ANNOUNCESLEEP);
		}

cleanup:
		printf(""); /* otherwise the pthread_cleanup_pop won't compile */
		pthread_cleanup_pop(1);
		return NULL;
}
