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
		smartshare.t0            = time(NULL);
		smartshare.prevhostid    = 0;
		smartshare.prevhostcount = 0;
		pthread_mutex_unlock(&mutexsmartshare);
}

/* use a probabilistic approach to determine whether the connection should be dropped */
int smartshare_check(float t, int hostid) {
		int verbose = 0;
		float p, r, baseline = 1;
		smartsharelist_t *listitem;

		if (verbose)
				fprintf(stderr, "smartshare_check\n");

		pthread_mutex_lock(&mutexsmartshare);

		/* always accept jobs when smartshare is disabled */
		if (smartshare.enabled!=1) {
				pthread_mutex_unlock(&mutexsmartshare);
				return 1;
		}

		/* always accept jobs when running in master mode */
		if (hoststatus()==STATUS_MASTER) {
				pthread_mutex_unlock(&mutexsmartshare);
				return 1;
		}

		/* always accept jobs that don't take any time, e.g. writing results back to the master */
		if (t<=0) {
				pthread_mutex_unlock(&mutexsmartshare);
				return 1;
		}

		if (smartshare.prevhostid==hostid)
				smartshare.prevhostcount++;
		else
				smartshare.prevhostcount = 0;

		smartshare.n++;
		smartshare.prevhostid = hostid;

		if (smartshare.prevhostcount >= SMARTSHARE_PREVHOSTCOUNT) {
				if (verbose)
						fprintf(stderr, "smartshare_check: prevhostcount exceeded\n");
				smartshare.prevhostcount = 0;
				pthread_mutex_unlock(&mutexsmartshare);
				return 1;
		}

		if ((time(NULL)-smartshare.t0) >= SMARTSHARE_TIMER) {
				if (verbose)
						fprintf(stderr, "smartshare_check: timer has elapsed\n");
				pthread_mutex_unlock(&mutexsmartshare);
				return 1;
		}

		/* determine the baseline for the time, based on the recent job history */
		if (smartsharelist) {
				listitem = smartsharelist;
				baseline = listitem->timreq;
				while (listitem) {
						if (listitem->timreq < baseline)
								baseline = listitem->timreq;
						listitem = listitem->next;
				}
		}
		else {
				baseline = t;
		}  

		/* scale the time of this job request with the baseline */
		t = t/baseline;

		if (t<=1)
				p = 1;
		else
				p = 1.0 / t;

		r = (float)rand() / (float)INT32_MAX;

		if (verbose)
				fprintf(stderr, "smartshare_check: t = %f, p = %f, r = %f, n = %d\n", t, p, r, smartshare.n);

		pthread_mutex_unlock(&mutexsmartshare);

		/* return 1 if the connection should be accepted, 0 if it should not be accepted */
		p = (p>r);
		if (verbose>0)
				fprintf(stderr, "smartshare_check: return value = %d\n", p);
		return p;
}


/* keep a short history of the jobs that are currently submidded */
void smartshare_history(jobdef_t *job) {
		int verbose = 0;
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

		if (verbose)
				fprintf(stderr, "smartshare_history: historycount = %d, peercount = %d\n", historycount, peercount);

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
