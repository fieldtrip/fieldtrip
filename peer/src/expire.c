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
        DEBUG(LOG_DEBUG, "cleanup_expire()");

		if (expireStatus==0)
				return;

		pthread_mutex_lock(&mutexstatus);
		expireStatus = 0;
		pthread_mutex_unlock(&mutexstatus);

		pthread_mutex_lock(&mutexthreadcount);
		threadcount--;
		pthread_mutex_unlock(&mutexthreadcount);
}

void *expire(void *arg) {
		int found;
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
				/* signal that this thread has started */
				pthread_cond_signal(&condstatus);
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
								DEBUG(LOG_INFO, "expire: name = %s", peerlist->host->name);
								DEBUG(LOG_INFO, "expire: port = %u", peerlist->host->port);
								DEBUG(LOG_INFO, "expire: time = %s", ctime(&(peerlist->time)));
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
										DEBUG(LOG_INFO, "expire: name = %s", next->host->name);
										DEBUG(LOG_INFO, "expire: port = %u", next->host->port);
										DEBUG(LOG_INFO, "expire: time = %s", ctime(&(next->time)));
										/* delete the next item in the list */
										peer->next = next->next;
										FREE(next->host);
										FREE(next);
										continue;
								}
						}
						peer = peer->next;
				}

				pthread_mutex_lock(&mutexkillswitch);

/*
fprintf(stderr, "killswitch: time = %d, now = %d\n", killswitch.time, time(NULL));
fprintf(stderr, "killswitch: masterid = %d\n", killswitch.masterid);
*/

				/* check whether the kill switch should be triggered for the time */
				if (killswitch.enabled && killswitch.time) {
						if (difftime(time(NULL), killswitch.time)>0) {
								/* the maximum allowed time has elapsed */
								DEBUG(LOG_CRIT, "expire: killswitch triggered (time)");
								exit(0);
						}
				}

				/* check whether the kill switch should be triggered for the masterid */
				if (killswitch.enabled && killswitch.masterid) {
						found = 0;

						/* look whether the master is stil available */
						peer = peerlist;
						while(peer && !found) {
								found = 1;
								found = found && (peer->host->id == killswitch.masterid);
								found = found && (peer->host->status == STATUS_MASTER);
								peer = peer->next;
						}

						if (!found)
								killswitch.evidence++;
						else
								killswitch.evidence==0;

						if (killswitch.evidence>2) {
								/* the master is not available any more */
								DEBUG(LOG_CRIT, "expire: killswitch triggered (master)");
								exit(0);
						}

				} /* if killswitch enabled */

				pthread_mutex_unlock(&mutexkillswitch);
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

