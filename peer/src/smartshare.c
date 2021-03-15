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

#include <math.h>
#include <time.h>

#include "peer.h"
#include "extern.h"
#include "platform_includes.h"

/* reset the timer */
void smartshare_reset(void) {
		pthread_mutex_lock(&mutexsmartshare);
		smartshare.n             = 0;
		smartshare.time          = time(NULL);
		smartshare.prevhostid    = 0;
		smartshare.prevhostcount = 0;
		pthread_mutex_unlock(&mutexsmartshare);
}

/* use a probabilistic approach to determine whether the connection should be dropped */
int smartshare_check(float t, int hostid) {
		float p, r;
		UINT64_T mintimreq;
		smartsharelist_t *listitem;

		DEBUG(LOG_DEBUG, "smartshare_check()");

		pthread_mutex_lock(&mutexsmartshare);

		/* always accept jobs when smartshare is disabled */
		if (smartshare.enabled!=1) {
				pthread_mutex_unlock(&mutexsmartshare);
				return 1;
		}

		/* accept the job when running in controller mode */
		if (hoststatus()==STATUS_CONTROLLER) {
				pthread_mutex_unlock(&mutexsmartshare);
				return 1;
		}

		/* accept the job if it does not take any time, e.g. writing results back to the controller */
		if (t<=0) {
				pthread_mutex_unlock(&mutexsmartshare);
				return 1;
		}

        /* count the number of subsequent requests from the same host */
		if (smartshare.prevhostid==hostid)
				smartshare.prevhostcount++;
		else
				smartshare.prevhostcount = 0;

		smartshare.n++;
		smartshare.prevhostid = hostid;

		/* accept the job if all previous requests originated from the same host */
		if (smartshare.prevhostcount >= SMARTSHARE_PREVHOSTCOUNT) {
						DEBUG(LOG_DEBUG, "smartshare_check: prevhostcount exceeded");
				smartshare.prevhostcount = 0;
				pthread_mutex_unlock(&mutexsmartshare);
				return 1;
		}

		/* accept the job if all previous jobs were too long ago */
		if (difftime(time(NULL), smartshare.time) > SMARTSHARE_TIMEOUT) {
				DEBUG(LOG_DEBUG, "smartshare_check: timer has elapsed");
				pthread_mutex_unlock(&mutexsmartshare);
				return 1;
		}

		/* determine the baseline for the time, based on the recent job history */
		if (smartsharelist) {
				listitem = smartsharelist;
				mintimreq = (float)listitem->timreq;
				while (listitem) {
						if (listitem->timreq < mintimreq)
								mintimreq = listitem->timreq;
						listitem = listitem->next;
				}
				/* scale the time of this job request with the minimal time required */
				t = t/((float)mintimreq);
		}
		else {
				/* the scale factor cannot be determined from the list of known jobs */
				t = 1;
		}  

		/* compute the probability of accepting the job */
		if (t<=1)
				p = 1.0;
		else
				p = 1.0 / t;

		/* compute random number between 0 and 1 */
		r = ((float)rand()) / ((float)RAND_MAX);

		DEBUG(LOG_DEBUG, "smartshare_check: t = %f, p = %f, r = %f, n = %d", t, p, r, smartshare.n);

		pthread_mutex_unlock(&mutexsmartshare);

		/* return 1 if the connection should be accepted, 0 if it should not be accepted */
		p = (p>r);
		DEBUG(LOG_INFO, "smartshare_check: return value = %d", p);
		return p;
}


/* keep a short history of the jobs that are currently submidded */
void smartshare_history(jobdef_t *job) {
		int historycount = 0;
		int peercount = 0;
		smartsharelist_t *listitem;
		peerlist_t *peer;

		listitem = malloc(sizeof(smartsharelist_t));
		listitem->timreq = job->timreq;

		pthread_mutex_lock(&mutexsmartshare);
		if (smartsharelist==NULL) {
				listitem->next = NULL;
				smartsharelist = listitem;
		}
		else {
				listitem->next = smartsharelist;
				smartsharelist = listitem;
		}

		/* remember the time at which the last job request was observed */
		smartshare.time = time(NULL);

		/* count the number of history items */
		listitem = smartsharelist;
		while (listitem) {
				historycount++;
				listitem = listitem->next;
		}

		/* count the number of peers */
		pthread_mutex_lock(&mutexpeerlist);
		peer = peerlist;
		while (peer) {
				peercount++;
				peer = peer->next;
		}
		pthread_mutex_unlock(&mutexpeerlist);

		DEBUG(LOG_DEBUG, "smartshare_history: historycount = %d, peercount = %d", historycount, peercount);

		while (historycount > peercount*SMARTSHARE_HISTORY) {
				/* remove the oldest item from the history */
				listitem = smartsharelist;
				while (listitem && listitem->next && listitem->next->next)
						listitem = listitem->next;

				FREE(listitem->next);
				listitem->next = NULL;
				historycount--;
		}

		pthread_mutex_unlock(&mutexsmartshare);
}
