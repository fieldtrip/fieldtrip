#include <arpa/inet.h>
#include <netinet/in.h>
#include <stdio.h>
#include <string.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h> /* for sleep */

#include "peer.h"
#include "extern.h"

void cleanup_expire(void *arg) {
		pthread_mutex_lock(&mutexstatus);
		expireStatus = 0;
		pthread_mutex_unlock(&mutexstatus);

		pthread_mutex_lock(&mutexthreadcount);
		threadcount--;
		pthread_mutex_unlock(&mutexthreadcount);
}

void *expire(void *arg) {
		int found, verbose = 0;
		peerlist_t *peer = NULL, *next = NULL;

		/* this is for debugging */
		pthread_mutex_lock(&mutexthreadcount);
		threadcount++;
		pthread_mutex_unlock(&mutexthreadcount);

		pthread_cleanup_push(cleanup_expire, NULL);

		/* the status contains the thread id when running, or zero when not running */
		pthread_mutex_lock(&mutexstatus);
		if (expireStatus==0) {
				expireStatus = 1;
				pthread_mutex_unlock(&mutexstatus);
		}
		else {
				pthread_mutex_unlock(&mutexstatus);
				goto cleanup;
		}

		/* now just enter a read-print loop */
		while (1) {
				pthread_mutex_lock(&mutexpeerlist);

				/* remove expired observations from the list */
				/* test the first item on the list */
				if (peerlist)
				{
						found = (difftime(time(NULL), peerlist->time)>EXPIRATION);
						found = found | !ismember_userlist (peerlist->host->user);
						found = found | !ismember_grouplist(peerlist->host->group);
						found = found | !ismember_hostlist (peerlist->host->name);
						if (found) {
								if (verbose>0) {
										fprintf(stderr, "expire: name = %s\n", peerlist->host->name);
										fprintf(stderr, "expire: addr = %s\n", peerlist->host->addr);
										fprintf(stderr, "expire: port = %d\n", peerlist->host->port);
										fprintf(stderr, "expire: time = %s", ctime(&(peerlist->time)));
								}
								/* delete the first item in the list */
								next = peerlist->next;
								FREE(peerlist->host);
								FREE(peerlist);
								peerlist = next;
						}
				}

				/* remove expired observations from the list */
				/* traverse the list */
				peer = peerlist;
				while(peer) {
						/* test the next item on the list, remember the current one */
						next = peer->next;
						if (next) {
								found = (difftime(time(NULL), next->time)>EXPIRATION);
								found = found | !ismember_userlist (next->host->user);
								found = found | !ismember_grouplist(next->host->group);
								found = found | !ismember_hostlist (next->host->name);
								if (found) {
										if (verbose>0) {
												fprintf(stderr, "expire: name = %s\n", next->host->name);
												fprintf(stderr, "expire: addr = %s\n", next->host->addr);
												fprintf(stderr, "expire: port = %d\n", next->host->port);
												fprintf(stderr, "expire: time = %s", ctime(&(next->time)));
										}
										/* delete the next item in the list */
										peer->next = next->next;
										FREE(next->host);
										FREE(next);
										break;
								}
						}
						peer = peer->next;
				}

				pthread_mutex_unlock(&mutexpeerlist);

				/* note that this is a thread cancelation point */
				pthread_testcancel();
				usleep(EXPIRESLEEP);

		} /* while (1) */

cleanup:
		printf(""); /* otherwise the pthread_cleanup_pop won't compile */
		pthread_cleanup_pop(1);
		return NULL;
}
