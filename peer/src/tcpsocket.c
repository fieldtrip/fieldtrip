/*
 * Copyright (C) 2008-2010, Robert Oostenveld
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

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "peer.h"
#include "extern.h"
#include "platform_includes.h"

typedef struct {
		void **message;
		int *fd;
} threadlocal_t;

void cleanup_tcpsocket(void *arg) {
		threadlocal_t *threadlocal;
		threadlocal = (threadlocal_t *)arg;

		if (threadlocal && *threadlocal->message) {
/*
                these pointers are added to the joblist and therefore should be left alone
				FREE(((message_t *)(*threadlocal->message))->host);
				FREE(((message_t *)(*threadlocal->message))->job);
*/
				FREE(*threadlocal->message);
		}

		if (threadlocal && (*threadlocal->fd)>0) {
				closesocket(*threadlocal->fd);
				*threadlocal->fd = 0;
		}

		pthread_mutex_lock(&mutexsocketcount);
		socketcount--;
		pthread_mutex_unlock(&mutexsocketcount);

		pthread_mutex_lock(&mutexthreadcount);
		threadcount--;
		pthread_mutex_unlock(&mutexthreadcount);
}

/* this function deals with the incoming message */
/* the return value is always NULL */
void *tcpsocket(void *arg) {
		int n, jobcount;
		int connect_accept = 1, connect_continue = 1, handshake;
		joblist_t *job;

		/* these are used for communication over the TCP socket */
		int fd = 0;
		message_t *message = NULL;

		threadlocal_t threadlocal;
		threadlocal.message = &message;
		threadlocal.fd      = &fd;

		/* the connection to the client has been made by the server */
		fd = ((int)arg);

		pthread_cleanup_push(cleanup_tcpsocket, &threadlocal);

		/* this is for debugging */
		pthread_mutex_lock(&mutexthreadcount);
		threadcount++;
		pthread_mutex_unlock(&mutexthreadcount);

		/* this is for debugging */
		pthread_mutex_lock(&mutexsocketcount);
		socketcount++;
		pthread_mutex_unlock(&mutexsocketcount);

		/* prevent smartcpu_update during the negociation and reading of the job */
		pthread_mutex_lock(&mutexsmartcpu);
		smartcpu.freeze = 1;
		pthread_mutex_unlock(&mutexsmartcpu);

		/* prevent smartmem_update during the negociation and reading of the job */
		pthread_mutex_lock(&mutexsmartmem);
		smartmem.freeze = 1;
		pthread_mutex_unlock(&mutexsmartmem);

		DEBUG(LOG_DEBUG, "tcpsocket: fd = %d, socketcount = %d, threadcount = %d", fd, socketcount, threadcount);

		pthread_mutex_lock(&mutexjoblist);
		jobcount = 0;
		job = joblist;
		while (job) {
				jobcount++;
				job = job->next;
		}
		pthread_mutex_unlock(&mutexjoblist);

		pthread_mutex_lock(&mutexhost);
		if (host->status==STATUS_MASTER) {
				connect_accept = 1;
		}
		else if (host->status==STATUS_IDLE && jobcount==0) {
				connect_accept = 1;
		}
		else {
				DEBUG(LOG_INFO, "tcpsocket: failed on status (%d) and/or jobcount (%d)", host->status, jobcount);
				connect_accept = 0;
		}
		pthread_mutex_unlock(&mutexhost);

		/* give a handshake */
		handshake = connect_accept || connect_continue;
		if ((n = bufwrite(fd, &handshake, sizeof(int))) != sizeof(int)) {
				DEBUG(LOG_ERR, "tcpsocket: could not write handshake, n = %d, should be %d", n, sizeof(int));
				goto cleanup;
		}

		if (!connect_continue) {
				DEBUG(LOG_INFO, "tcpsocket: dropping connection");
				goto cleanup;
		}

		message       = (message_t*)malloc(sizeof(message_t));
		message->host = (hostdef_t*)malloc(sizeof(hostdef_t));
		message->job  = (jobdef_t*)malloc(sizeof(jobdef_t));
		message->arg  = NULL;
		message->opt  = NULL;

		/* read the host details */
		if ((n = bufread(fd, message->host, sizeof(hostdef_t))) != sizeof(hostdef_t)) {
				DEBUG(LOG_DEBUG, "tcpsocket: read size = %d, should be %d", n, sizeof(hostdef_t));
				goto cleanup;
		}

		/* test whether the version is compatible */
		if (message->host->version!=VERSION) {
				DEBUG(LOG_ERR, "tcpsocket: incorrect host version (%u, %u)", message->host->version, VERSION);
				connect_accept   = 0;
				connect_continue = 0; /* prevent another read request */
		}

		/* determine whether the host, group and user is allowed to execute a job */
		if (message->host->version==VERSION && !security_check(message->host)) {
				DEBUG(LOG_INFO, "tcpsocket: failed security check");
				DEBUG(LOG_DEBUG, "tcpsocket: ismember_userlist  = %d", ismember_userlist(message->host->user));
				DEBUG(LOG_DEBUG, "tcpsocket: ismember_grouplist = %d", ismember_grouplist(message->host->group));
				DEBUG(LOG_DEBUG, "tcpsocket: ismember_hostlist  = %d", ismember_hostlist(message->host->name));
				connect_accept = 0;
		}

		/* give a handshake */
		handshake = connect_accept || connect_continue;
		if ((n = bufwrite(fd, &handshake, sizeof(int))) != sizeof(int)) {
				DEBUG(LOG_ERR, "tcpsocket: could not write handshake, n = %d, should be %d", n, sizeof(int));
				goto cleanup;
		}

		if (!connect_continue) {
				DEBUG(LOG_INFO, "tcpsocket: dropping connection");
				goto cleanup;
		}

		/* read the job details */
		if ((n = bufread(fd, message->job, sizeof(jobdef_t))) != sizeof(jobdef_t)) {
				DEBUG(LOG_ERR, "tcpsocket: packet size = %d, should be %d", n, sizeof(jobdef_t));
				goto cleanup;
		}

		/* test whether the request can be accepted based on the job characteristics */
		if (message->job->version!=VERSION) {
				DEBUG(LOG_ERR, "tcpsocket: incorrect job version (%u, %u)", message->job->version, VERSION);
				connect_accept   = 0;
				connect_continue = 0; /* prevent another read request */
		}

		pthread_mutex_lock(&mutexhost);
		if (message->job->memreq > host->memavail) {
				DEBUG(LOG_INFO, "tcpsocket: memory request too large");
				connect_accept = 0;
		}
		if (message->job->cpureq > host->cpuavail) {
				DEBUG(LOG_INFO, "tcpsocket: cpu request too large");
				connect_accept = 0;
		}
		if (message->job->timreq > host->timavail) {
				DEBUG(LOG_INFO, "tcpsocket: time request too large");
				connect_accept = 0;
		}
		pthread_mutex_unlock(&mutexhost);

		if (message->job->argsize>MAXARGSIZE) {
				DEBUG(LOG_ERR, "tcpsocket: argsize too large");
				connect_accept = 0;
		}

		if (message->job->optsize>MAXARGSIZE) {
				DEBUG(LOG_ERR, "tcpsocket: optsize too large");
				connect_accept = 0;
		}

		/* remember the job characteristics for the smartshare algorithm */
		smartshare_history(message->job);

		/* use a probabilistic approach to determine whether the connection should be dropped */
		if (!smartshare_check(message->job->timreq, message->host->id)) {
				DEBUG(LOG_INFO, "tcpsocket: failed smartshare_check");
				connect_accept = 0;
		}

		/* don't continue reading the content of the job, drop the connection before the job arguments are sent */
		if (!connect_accept)
				connect_continue = 0; /* prevent another read request */

		/* give a handshake */
		handshake = connect_accept || connect_continue;
		if ((n = bufwrite(fd, &handshake, sizeof(int))) != sizeof(int)) {
				DEBUG(LOG_ERR, "tcpsocket: could not write handshake, n = %d, should be %d", n, sizeof(int));
				goto cleanup;
		}

		if (!connect_continue) {
				DEBUG(LOG_INFO, "tcpsocket: dropping connection");
				goto cleanup;
		}

		/* read the job request arguments */
		if (message->job->argsize>0) {
				message->arg = malloc(message->job->argsize);
				if ((n = bufread(fd, message->arg, message->job->argsize)) != message->job->argsize) {
						DEBUG(LOG_ERR, "tcpsocket: read size = %d, should be %d", n, message->job->argsize);
						goto cleanup;
				}
		}

		/* give a handshake */
		handshake = connect_accept || connect_continue;
		if ((n = bufwrite(fd, &handshake, sizeof(int))) != sizeof(int)) {
				DEBUG(LOG_ERR, "tcpsocket: could not write handshake, n = %d, should be %d", n, sizeof(int));
				goto cleanup;
		}

		if (!connect_continue) {
				DEBUG(LOG_INFO, "tcpsocket: dropping connection");
				goto cleanup;
		}

		/* read the job request options */
		if (message->job->optsize>0) {
				message->opt = malloc(message->job->optsize);
				if ((n = bufread(fd, message->opt, message->job->optsize)) != message->job->optsize) {
						DEBUG(LOG_ERR, "tcpsocket: read size = %d, should be %d", n, message->job->optsize);
						goto cleanup;
				}
		}

		/* give a handshake */
		handshake = connect_accept || connect_continue;
		if ((n = bufwrite(fd, &handshake, sizeof(int))) != sizeof(int)) {
				DEBUG(LOG_ERR, "tcpsocket: could not write handshake, n = %d, should be %d", n, sizeof(int));
				goto cleanup;
		}

		if (!connect_continue) {
				DEBUG(LOG_INFO, "tcpsocket: dropping connection");
				goto cleanup;
		}

		/* create a new list item */
		job = (joblist_t *)malloc(sizeof(joblist_t));
		job->job  = message->job;
		job->host = message->host;
		job->arg  = message->arg;
		job->opt  = message->opt;

		pthread_mutex_lock(&mutexjoblist);
		/* add the item to the beginning of the list */
		job->next = joblist;
		joblist = job;

		DEBUG(LOG_DEBUG, "tcpsocket: job.version  = %u", job->job->version);
		DEBUG(LOG_DEBUG, "tcpsocket: job.id       = %u", job->job->id);
		DEBUG(LOG_DEBUG, "tcpsocket: job.argsize  = %u", job->job->argsize);
		DEBUG(LOG_DEBUG, "tcpsocket: job.optsize  = %u", job->job->optsize);
		DEBUG(LOG_DEBUG, "tcpsocket: host.name    = %s", job->host->name);
		DEBUG(LOG_DEBUG, "tcpsocket: host.port    = %u", job->host->port);
		DEBUG(LOG_DEBUG, "tcpsocket: host.id      = %u", job->host->id);

		pthread_mutex_unlock(&mutexjoblist);

cleanup:

		/* from now on it is again allowed to use smartcpu_update */
		pthread_mutex_lock(&mutexsmartcpu);
		smartcpu.freeze = 0;
		pthread_mutex_unlock(&mutexsmartcpu);

		/* from now on it is again allowed to use smartcpu_update */
		pthread_mutex_lock(&mutexsmartmem);
		smartmem.freeze = 0;
		pthread_mutex_unlock(&mutexsmartmem);

		printf(""); /* otherwise the pthread_cleanup_pop won't compile */
		pthread_cleanup_pop(1);
		return NULL;
}

