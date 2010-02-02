#include <math.h>
#include <time.h>

#include "peer.h"
#include "extern.h"

/* reset the timer */
void fairshare_timer(void) {

		pthread_mutex_lock(&mutexfairshare);
		param.t0 = time(NULL);
		param.prevhostid    = 0;
		param.prevhostcount = 0;
		param.n  = 0;
		pthread_mutex_unlock(&mutexfairshare);

}

/* use a probabilistic approach to determine whether the connection should be dropped */
int fairshare_check(float t, int hostid) {
		int verbose = 0, retval;
		float p, r;

		if (t<=0) {
				/* always accept jobs that don't take any time, e.g. writing results back to the master */
				return 1;
		}

		pthread_mutex_lock(&mutexfairshare);
		param.n++;

		if (param.prevhostid==hostid)
				param.prevhostcount++;
		else
				param.prevhostcount = 0;
		param.prevhostid = hostid;

		if (param.prevhostcount >= param.c) {
				if (verbose)
						fprintf(stderr, "fairshare: prevhostcount exceeded\n");
				pthread_mutex_unlock(&mutexfairshare);
				return 1;
		}

		if ((time(NULL)-param.t0) >= param.d) {
				if (verbose)
						fprintf(stderr, "fairshare: timer has elapsed\n");
				pthread_mutex_unlock(&mutexfairshare);
				return 1;
		}

		r = (float)random() / (float)INT32_MAX;
		t = t/param.a;

		if (t<=1)
				p = 1;
		else
				p = powf(t, param.b);

		if (verbose)
				fprintf(stderr, "fairshare: t = %f, p = %f, r = %f, n = %d\n", t, p, r, param.n);

		pthread_mutex_unlock(&mutexfairshare);

		/* return 1 if the connection should be accepted, 0 if it should not be accepted */
		return (p>r);
}

