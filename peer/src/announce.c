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
		struct sockaddr_in multicast, broadcast;
		hostdef_t *message = NULL;
		unsigned char ttl = 3;
		unsigned char one = 1;
		int use_broadcast, use_multicast;

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
		if ((fd=socket(AF_INET,SOCK_DGRAM,0)) < 0) {
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

		/* set up destination address for broadcasting */
		/* only broadcast to localhost on which multiple peers may be listening */
		memset(&broadcast,0,sizeof(broadcast));
		broadcast.sin_family      = AF_INET;
		broadcast.sin_addr.s_addr = inet_addr("127.0.0.1");
		broadcast.sin_port        = htons(ANNOUNCE_PORT);

		/* now just sendto() our destination */
		while (1) {
				/* the host port and status are variable */
				pthread_mutex_lock(&mutexhost);
				memcpy(message, host, sizeof(hostdef_t));
				use_broadcast = (strncasecmp(host->name, "localhost", STRLEN)==0);
				use_multicast = (strncasecmp(host->name, "localhost", STRLEN)==1);
				pthread_mutex_unlock(&mutexhost);

				if (verbose>1)
						fprintf(stderr, "announce: use_broadcast = %d, use_multicast = %d\n", use_broadcast, use_multicast);

				if (use_broadcast)
						/* note that this is a thread cancelation point */
						if (sendto(fd,message,sizeof(hostdef_t),0,(struct sockaddr *) &multicast,sizeof(multicast)) < 0) {
								perror("announce sendto");
								goto cleanup;
						}

				if (use_multicast)
						/* note that this is a thread cancelation point */
						if (sendto(fd,message,sizeof(hostdef_t),0,(struct sockaddr *) &multicast,sizeof(multicast)) < 0) {
								perror("announce sendto");
								goto cleanup;
						}

				if (verbose>1)
						fprintf(stderr, "announce\n");

				/* note that this is a thread cancelation point */
				pthread_testcancel();
				usleep(ANNOUNCESLEEP);
		}

cleanup:
		printf(""); /* otherwise the pthread_cleanup_pop won't compile */
		pthread_cleanup_pop(1);
		return NULL;
}
