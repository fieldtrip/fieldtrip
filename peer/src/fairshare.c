#include <math.h>
#include <time.h>

#include "peer.h"
#include "extern.h"
#include "platform_includes.h"

/* reset the timer */
void fairshare_reset(void) {
		pthread_mutex_lock(&mutexfairshare);
		fairshare.n             = 0;
		fairshare.t0            = time(NULL);
		fairshare.prevhostid    = 0;
		fairshare.prevhostcount = 0;
		pthread_mutex_unlock(&mutexfairshare);
}

/* use a probabilistic approach to determine whether the connection should be dropped */
int fairshare_check(float t, int hostid) {
		int verbose = 0;
		float p, r, baseline = 1;
		fairsharelist_t *listitem;

		pthread_mutex_lock(&mutexfairshare);

		/* always accept jobs when fairshare is disabled */
		if (fairshare.enabled!=1) {
				pthread_mutex_unlock(&mutexfairshare);
				return 1;
		}

		/* always accept jobs when running in master mode */
		if (hoststatus()==2) {
				pthread_mutex_unlock(&mutexfairshare);
				return 1;
		}

		/* always accept jobs that don't take any time, e.g. writing results back to the master */
		if (t<=0) {
				pthread_mutex_unlock(&mutexfairshare);
				return 1;
		}

		if (fairshare.prevhostid==hostid)
				fairshare.prevhostcount++;
		else
				fairshare.prevhostcount = 0;

		fairshare.n++;
		fairshare.prevhostid = hostid;

		if (fairshare.prevhostcount >= FAIRSHARE_PREVHOSTCOUNT) {
				if (verbose)
						fprintf(stderr, "fairshare: prevhostcount exceeded\n");
				fairshare.prevhostcount = 0;
				pthread_mutex_unlock(&mutexfairshare);
				return 1;
		}

		if ((time(NULL)-fairshare.t0) >= FAIRSHARE_TIMER) {
				if (verbose)
						fprintf(stderr, "fairshare: timer has elapsed\n");
				pthread_mutex_unlock(&mutexfairshare);
				return 1;
		}

		/* determine the baseline for the time, based on the recent job history */
		if (fairsharelist) {
				listitem = fairsharelist;
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
				p = powf(t, -1);

		r = (float)random() / (float)INT32_MAX;

		if (verbose)
				fprintf(stderr, "fairshare: t = %f, p = %f, r = %f, n = %d\n", t, p, r, fairshare.n);

		pthread_mutex_unlock(&mutexfairshare);

		/* return 1 if the connection should be accepted, 0 if it should not be accepted */
		return (p>r);
}


/* keep a short history of the jobs that are currently submidded */
void fairshare_history(jobdef_t *job) {
		int verbose = 0;
		int historycount = 0;
		int peercount = 0;
		fairsharelist_t *listitem;
		peerlist_t *peer;

		listitem = malloc(sizeof(fairsharelist_t));
		listitem->timreq = job->timreq;

		pthread_mutex_lock(&mutexfairshare);
		if (fairsharelist==NULL) {
				listitem->next = NULL;
				fairsharelist = listitem;
		}
		else {
				listitem->next = fairsharelist;
				fairsharelist = listitem;
		}

		/* count the number of history items */
		listitem = fairsharelist;
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
				fprintf(stderr, "historycount = %d, peercount = %d\n", historycount, peercount);

		while (historycount > peercount*FAIRSHARE_HISTORY) {
				/* remove the oldest item from the history */
				listitem = fairsharelist;
				while (listitem && listitem->next && listitem->next->next)
						listitem = listitem->next;

				FREE(listitem->next);
				listitem->next = NULL;
				historycount--;
		}

		pthread_mutex_unlock(&mutexfairshare);
}
