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

#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>

#include "peer.h"
#include "extern.h"
#include "platform_includes.h"

void cleanup_expire(void *arg) {
        syslog(LOG_DEBUG, "cleanup_expire()");

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

		syslog(LOG_NOTICE, "expire()");

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
								syslog(LOG_INFO, "expire: name = %s", peerlist->host->name);
								syslog(LOG_INFO, "expire: port = %d", peerlist->host->port);
								syslog(LOG_INFO, "expire: time = %s", ctime(&(peerlist->time)));
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
										syslog(LOG_INFO, "expire: name = %s\n", next->host->name);
										syslog(LOG_INFO, "expire: port = %d\n", next->host->port);
										syslog(LOG_INFO, "expire: time = %s", ctime(&(next->time)));
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
