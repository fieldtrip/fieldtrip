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
#include <stdlib.h>
#include <sys/types.h>     /* for getpwuid and geteuid */

#include "peer.h"
#include "extern.h"
#include "platform_includes.h"


unsigned long hash(unsigned char *str)
{
		unsigned long hash = 5381;
		int c;

		while (c = *str++)
				hash = ((hash << 5) + hash) + c; /* hash * 33 + c */

		return hash;
}

/******************************************************************************/
void peerinit(void *arg) {
		struct passwd *pwd;
		int s, verbose = 0;

#if defined(PLATFORM_WIN32) || defined(PLATFORM_WIN64)
		DWORD nStrLen; 
#endif

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

		/* initialize the unix domain socket as empty */
		bzero(host->socket, STRLEN);

#if defined (PLATFORM_LINUX) || defined (PLATFORM_OSX)

		/* get the user name */
		pwd = getpwuid(geteuid());
		if (pwd==NULL) {
				fprintf(stderr, "error: announce getpwuid\n");
				exit(1);
		}
		strncpy(host->user, pwd->pw_name, STRLEN);

		/* set the default group name */
		strncpy(host->group, DEFAULT_GROUP, STRLEN);

		/* get the host name */
		if (gethostname(host->name, STRLEN)) {
				perror("announce gethostname");
				exit(1);
		}

#else

		/* set the user name */
		nStrLen = STRLEN;
		GetUserName(host->user, &nStrLen);

		/* set the default group name */
		strncpy(host->group, DEFAULT_GROUP, STRLEN);

		/* get the host name */
		nStrLen = STRLEN;
		GetComputerName(host->name, &nStrLen);

#endif

		/* initialize the random number generator */
		/* this is used for creating host and job IDs */
		srand(getpid());

		/* this turns out not to be sufficiently random in the case of many peers */
		host->id       = rand();

		/* increasse the host ID with some host specific details */
		host->id += hash(host->name);
		host->id += hash(host->user);
		host->id += hash(host->group);

		if (verbose>0) {
				fprintf(stderr, "peerinit: host.name =  %s\n", host->name);
				fprintf(stderr, "peerinit: host.port =  %d\n", host->port);
				fprintf(stderr, "peerinit: host.id   =  %d\n", host->id);
		}

		pthread_mutex_unlock(&mutexhost);

		pthread_mutex_lock(&mutexsmartmem);
		smartmem_enabled  = 0;
		pthread_mutex_unlock(&mutexsmartmem);

		pthread_mutex_lock(&mutexsmartcpu);
		smartcpu_enabled  = 0;
		pthread_mutex_unlock(&mutexsmartcpu);

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

		clear_peerlist();
		clear_joblist();
		clear_userlist();
		clear_grouplist();
		clear_hostlist();
		clear_fairsharelist();

		return;
}

