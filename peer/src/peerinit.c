#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>     /* for getpwuid and geteuid */

#include "peer.h"
#include "extern.h"
#include "platform_includes.h"

/******************************************************************************/
void peerinit(void *arg) {
		struct passwd *pwd;
		int family, s, verbose = 0;

		if (verbose)
				fprintf(stderr, "peerinit\n");

		pthread_mutex_lock(&mutexhost);
		if (host) {
				/* it should only be initialized once */
				pthread_mutex_unlock(&mutexhost);
				return;
		}

		/* check that all datatypes are of the expected size */
		check_datatypes();

		/* initialize the random number generator */
		/* this is used for creating host and job IDs */
		srand(time(NULL));

		if ((host = malloc(sizeof(hostdef_t)))==NULL) {
				perror("announce malloc");
				exit(1);
		}

		/* specify the host parameters */
		host->version  = VERSION;
		host->port     = DEFAULT_PORT;
		host->status   = DEFAULT_STATUS;
		host->memavail = DEFAULT_MEMAVAIL;
		host->cpuavail = DEFAULT_TIMAVAIL;
		host->timavail = DEFAULT_CPUAVAIL;
		host->id       = rand();

#if defined (PLATFORM_LINUX) || defined (PLATFORM_OSX)

		/* get the user name */
		pwd = getpwuid(geteuid());
		strncpy(host->user, pwd->pw_name, STRLEN);

		/* set the default group name */
		strncpy(host->group, DEFAULT_GROUP, STRLEN);

		/* get the host name */
		if (gethostname(host->name, STRLEN)) {
				perror("announce gethostname");
				exit(1);
		}


#else

		/* set the default user name */
		strncpy(host->group, DEFAULT_USER, STRLEN);

		/* set the default group name */
		strncpy(host->group, DEFAULT_GROUP, STRLEN);

		/* set the default group name */
		strncpy(host->group, DEFAULT_HOST, STRLEN);

#endif

		if (verbose>0) {
				fprintf(stderr, "peerinit: host.name =  %s\n", host->name);
				fprintf(stderr, "peerinit: host.port =  %d\n", host->port);
				fprintf(stderr, "peerinit: host.id   =  %d\n", host->id);
		}

		pthread_mutex_unlock(&mutexhost);

		pthread_mutex_lock(&mutexfairshare);
		fairshare.n             = 0;
		fairshare.t0            = time(NULL);
		fairshare.prevhostid    = 0;
		fairshare.prevhostcount = 0;
		fairshare.enabled       = 1;
		pthread_mutex_unlock(&mutexfairshare);

		return;
}

/******************************************************************************/
/* free the dynamical memory that is shared between the threads               */
void peerexit(void *arg) {
		int verbose = 0;
		peerlist_t *peer = NULL;
		joblist_t *job = NULL;
		userlist_t *user = NULL;
		grouplist_t *group = NULL;
		fairsharelist_t *listitem = NULL;

		if (verbose)
				fprintf(stderr, "peerexit\n");

		pthread_mutex_lock(&mutexhost);
		FREE(host);
		pthread_mutex_unlock(&mutexhost);

		pthread_mutex_lock(&mutexpeerlist);
		peer = peerlist;
		while (peer) {
				peerlist = peer->next;
				FREE(peer->host);
				FREE(peer);
				peer = peerlist;
		}
		pthread_mutex_unlock(&mutexpeerlist);

		pthread_mutex_lock(&mutexjoblist);
		job = joblist;
		while (job) {
				joblist = job->next;
				FREE(job->job);
				FREE(job->host);
				FREE(job->arg);
				FREE(job->opt);
				FREE(job);
				job = joblist;
		}
		pthread_mutex_unlock(&mutexjoblist);

		pthread_mutex_lock(&mutexuserlist);
		user = userlist;
		while (user) {
				userlist = user->next;
				FREE(user->name);
				FREE(user);
				user = userlist;
		}
		pthread_mutex_unlock(&mutexuserlist);

		pthread_mutex_lock(&mutexgrouplist);
		group = grouplist;
		while (group) {
				grouplist = group->next;
				FREE(group->name);
				FREE(group);
				group = grouplist;
		}
		pthread_mutex_unlock(&mutexgrouplist);

		pthread_mutex_lock(&mutexfairshare);
		listitem = fairsharelist;
		while (listitem) {
				fairsharelist = listitem->next;
				FREE(listitem);
				listitem = fairsharelist;
		}
		pthread_mutex_unlock(&mutexfairshare);

		return;
}

