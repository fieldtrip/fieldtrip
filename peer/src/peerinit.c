
#include <arpa/inet.h>
#include <ifaddrs.h>
#include <netdb.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/socket.h>
#include <unistd.h>
#include <sys/types.h>     /* for getpwuid and geteuid */
#include <pwd.h>           /* for getpwuid and geteuid */

#include "peer.h"
#include "extern.h"

/******************************************************************************/
void peerinit(void *arg) {
		struct ifaddrs *ifaddr, *ifa;
		struct passwd *pwd;
		int family, s, verbose = 0;
		char str[NI_MAXHOST];

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
		srandom(time(NULL));

		if ((host = malloc(sizeof(hostdef_t)))==NULL) {
				perror("announce malloc");
				exit(1);
		}

		/* get the user name */
		pwd = getpwuid(geteuid());

		/* specify the host parameters */
		host->version  = VERSION;
		host->port     = DEFAULT_PORT;
		host->status   = DEFAULT_STATUS;
		host->memavail = DEFAULT_MEMAVAIL;
		host->cpuavail = DEFAULT_TIMAVAIL;
		host->timavail = DEFAULT_CPUAVAIL;
		host->id       = random();
		strncpy(host->user, pwd->pw_name, STRLEN);
		strncpy(host->group, DEFAULT_GROUP, STRLEN);

		if (gethostname(host->name, STRLEN)) {
				perror("announce gethostname");
				exit(1);
		}

		if (getifaddrs(&ifaddr) == -1) {
				perror("peerinit getifaddrs");
				exit(1);
		}

		/* Walk through linked list, maintaining head pointer so we
		   can free list later */

		for (ifa = ifaddr; ifa != NULL; ifa = ifa->ifa_next) {
				family = ifa->ifa_addr->sa_family;

				/* Display interface name and family (including symbolic
				   form of the latter for the common families) */

				if (verbose>1)
						printf("%s  address family: %d%s\n",
										ifa->ifa_name, family,
										(family == AF_INET) ?   " (AF_INET)" :
										(family == AF_INET6) ?  " (AF_INET6)" : "");

				/* For an AF_INET* interface address, display the address */

				if (family == AF_INET || family == AF_INET6) {

						s = getnameinfo(ifa->ifa_addr, (family == AF_INET) ? sizeof(struct sockaddr_in) :sizeof(struct sockaddr_in6),
										str, NI_MAXHOST, NULL, 0, NI_NUMERICHOST);
						if (s != 0) {
								fprintf(stderr, "peerinfo: getnameinfo failed (%s)\n", gai_strerror(s));
								exit(1);
						}

						if (verbose>1)
								printf("\taddress: <%s>\n", str);

						/* there is only one address that we are interested in */
						if (strcmp(ifa->ifa_name, "lo0")!=0 && family==AF_INET)
								strncpy(host->addr, str, STRLEN);

				}
		}

		freeifaddrs(ifaddr);

		if (verbose>0) {
				fprintf(stderr, "peerinit: host.name =  %s\n", host->name);
				fprintf(stderr, "peerinit: host.addr =  %s\n", host->addr);
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

		if (verbose)
				fprintf(stderr, "peerexit\n");

		pthread_mutex_lock(&mutexhost);
		FREE(host);
		pthread_mutex_unlock(&mutexhost);

		pthread_mutex_lock(&mutexpeerlist);
		peerlist_t *peer = peerlist;
		while (peer) {
				peerlist = peer->next;
				FREE(peer->host);
				FREE(peer);
				peer = peerlist;
		}
		pthread_mutex_unlock(&mutexpeerlist);

		pthread_mutex_lock(&mutexjoblist);
		joblist_t *job = joblist;
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
		userlist_t *user = userlist;
		while (user) {
				userlist = user->next;
				FREE(user->name);
				FREE(user);
				user = userlist;
		}
		pthread_mutex_unlock(&mutexuserlist);

		pthread_mutex_lock(&mutexgrouplist);
		grouplist_t *group = grouplist;
		while (group) {
				grouplist = group->next;
				FREE(group->name);
				FREE(group);
				group = grouplist;
		}
		pthread_mutex_unlock(&mutexgrouplist);

		pthread_mutex_lock(&mutexfairshare);
		fairsharelist_t *item = fairsharelist;
		while (item) {
				fairsharelist = item->next;
				FREE(item);
				item = fairsharelist;
		}
		pthread_mutex_unlock(&mutexfairshare);

		return;
}

