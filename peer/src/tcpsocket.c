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

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "peer.h"
#include "extern.h"
#include "platform_includes.h"

typedef struct {
		void *message;
		int fd;
} threadlocal_t;

void cleanup_tcpsocket(void *arg) {
		threadlocal_t *threadlocal;
		threadlocal = (threadlocal_t *)arg;
		if (threadlocal && threadlocal->message) {
				FREE(threadlocal->message);
		}
		if (threadlocal && threadlocal->fd>0) {
				close(threadlocal->fd);
				threadlocal->fd = -1;
		}

		pthread_mutex_lock(&mutexsocketcount);
		socketcount--;
		pthread_mutex_unlock(&mutexsocketcount);

		pthread_mutex_lock(&mutexthreadcount);
		threadcount--;
		pthread_mutex_unlock(&mutexthreadcount);
}

/* this function deals with the incoming message */
void *tcpsocket(void *arg) {
		int n, verbose = 0;
		int connect_accept = 1, connect_continue = 1, handshake;
		joblist_t *job;

		/* these are used for communication over the TCP socket */
		int fd = 0;
		message_t *message = NULL;

		threadlocal_t threadlocal;
		threadlocal.message = NULL;
		threadlocal.fd = -1;

		/* the connection to the client has been made by the server */
		fd = ((int)arg);

		/* this will be closed at cleanup */
		threadlocal.fd = fd;

		pthread_cleanup_push(cleanup_tcpsocket, &threadlocal);

		/* this is for debugging */
		pthread_mutex_lock(&mutexthreadcount);
		threadcount++;
		pthread_mutex_unlock(&mutexthreadcount);

		/* this is for debugging */
		pthread_mutex_lock(&mutexsocketcount);
		socketcount++;
		/*
		   if (socketcount>100) {
		   pthread_mutex_unlock(&mutexsocketcount);
		   fprintf(stderr, "tcpsocket: too many open sockets\n");
		   goto cleanup;
		   }
		   else
		 */
		pthread_mutex_unlock(&mutexsocketcount);

		if (verbose>1)
				fprintf(stderr, "tcpsocket: hoststatus = %d, jobcount = %d\n", hoststatus(), jobcount());

		if (verbose>1)
				fprintf(stderr, "tcpsocket: fd = %d, socketcount = %d, threadcount = %d\n", fd, socketcount, threadcount);

		/* status = 0 means zombie mode, don't accept anything   */
		/* status = 1 means slave mode, accept only a single job */
		/* status = 2 means master mode, accept everything       */
		/* any other status is interpreted as zombie mode        */

		if (hoststatus()==1 && jobcount()==0) {
				connect_accept   = 1;
		}
		else if (hoststatus()==2) {
				connect_accept   = 1;
		}
		else {
				connect_accept   = 0;
		}

		/* give a handshake */
		handshake = connect_accept | connect_continue;
		if ((n = bufwrite(fd, &handshake, sizeof(int))) != sizeof(int)) {
				if (verbose>0)
						fprintf(stderr, "tcpsocket: could not write handshake, n = %d, should be %d\n", n, sizeof(int));
				goto cleanup;
		}

		if (!connect_accept)
				if (verbose>1)
						fprintf(stderr, "tcpsocket: dropping connection based on status\n");
		if (!connect_continue)
				goto cleanup;

		message       = (message_t*)malloc(sizeof(message_t));
		message->host = (hostdef_t*)malloc(sizeof(hostdef_t));
		message->job  = (jobdef_t*)malloc(sizeof(jobdef_t));
		message->arg  = NULL;
		message->opt  = NULL;

		/* this will be deallocated at cleanup */
		threadlocal.message = message;

		/* read the host details */
		if ((n = bufread(fd, message->host, sizeof(hostdef_t))) != sizeof(hostdef_t)) {
				if (verbose>0)
						fprintf(stderr, "tcpsocket: read size = %d, should be %d\n", n, sizeof(hostdef_t));
				goto cleanup;
		}

		if (verbose>1) {
				fprintf(stderr, "tcpsocket: ismember_userlist  = %d\n", ismember_userlist(message->host->user));
				fprintf(stderr, "tcpsocket: ismember_grouplist = %d\n", ismember_grouplist(message->host->group));
				fprintf(stderr, "tcpsocket: ismember_hostlist  = %d\n", ismember_hostlist(message->host->name));
		}

		/* test whether the version is compatible */
		if (message->host->version!=VERSION) {
				if (verbose>0)
						fprintf(stderr, "tcpsocket: incorrect version (%d, %d)\n", message->host->version, VERSION);
				connect_accept   = 0;
				connect_continue = 0;
		}

		/* determine whether the host, group and user is allowed to execute a job */
		if (!security_check(message->host)) {
				if (verbose>0)
						fprintf(stderr, "tcpsocket: security_check returned zero\n");
				connect_accept = 0;
		}

		/* give a handshake */
		handshake = connect_accept | connect_continue;
		if ((n = bufwrite(fd, &handshake, sizeof(int))) != sizeof(int)) {
				if (verbose>0)
						fprintf(stderr, "tcpsocket: could not write handshake, n = %d, should be %d\n", n, sizeof(int));
				goto cleanup;
		}

		if (!connect_accept)
				if (verbose>1)
						fprintf(stderr, "tcpsocket: dropping connection based on host details\n");
		if (!connect_continue)
				goto cleanup;

		/* read the job details */
		if ((n = bufread(fd, message->job, sizeof(jobdef_t))) != sizeof(jobdef_t)) {
				if (verbose>0)
						fprintf(stderr, "tcpsocket: packet size = %d, should be %d\n", n, sizeof(jobdef_t));
				goto cleanup;
		}

		/* test whether the request can be accepted based on the job characteristics */
		if (message->job->version!=VERSION) {
				if (verbose>0)
						fprintf(stderr, "tcpsocket: incorrect version\n");
				connect_accept   = 0;
				connect_continue = 0;
		}

		pthread_mutex_lock(&mutexhost);
		if (message->job->memreq > host->memavail) {
				if (verbose>0)
						fprintf(stderr, "tcpsocket: memory request too large\n");
				connect_accept = 0;
		}
		if (message->job->cpureq > host->cpuavail) {
				if (verbose>0)
						fprintf(stderr, "tcpsocket: cpu request too large\n");
				connect_accept = 0;
		}
		if (message->job->timreq > host->timavail) {
				if (verbose>0)
						fprintf(stderr, "tcpsocket: time request too large\n");
				connect_accept = 0;
		}
		pthread_mutex_unlock(&mutexhost);

		if (message->job->argsize>MAXARGSIZE) {
				if (verbose>0)
						fprintf(stderr, "tcpsocket: argsize too large\n");
				connect_accept = 0;
		}

		if (message->job->optsize>MAXARGSIZE) {
				if (verbose>0)
						fprintf(stderr, "tcpsocket: optsize too large\n");
				connect_accept = 0;
		}

		/* remember the job characteristics for the fairshare algorithm */
		fairshare_history(message->job);

		/* use a probabilistic approach to determine whether the connection should be dropped */
		if (!fairshare_check(message->job->timreq, message->host->id)) {
				if (verbose>0)
						fprintf(stderr, "tcpsocket: fairshare_check returned zero\n");
				connect_accept = 0;
		}

		/* don't continue reading the content of the job, drop the connection before the job arguments are sent */
		if (!connect_accept)
				connect_continue = 0;

		/* give a handshake */
		handshake = connect_accept | connect_continue;
		if ((n = bufwrite(fd, &handshake, sizeof(int))) != sizeof(int)) {
				if (verbose>0)
						fprintf(stderr, "tcpsocket: could not write handshake, n = %d, should be %d\n", n, sizeof(int));
				goto cleanup;
		}

		if (!connect_accept)
				if (verbose>1)
						fprintf(stderr, "tcpsocket: dropping connection based on job details\n");
		if (!connect_continue)
				goto cleanup;

		/* read the job request arguments */
		if (message->job->argsize>0) {
				message->arg = malloc(message->job->argsize);
				if ((n = bufread(fd, message->arg, message->job->argsize)) != message->job->argsize) {
						if (verbose>0)
								fprintf(stderr, "tcpsocket: read size = %d, should be %d\n", n, message->job->argsize);
						goto cleanup;
				}
		}

		/* give a handshake */
		handshake = connect_accept | connect_continue;
		if ((n = bufwrite(fd, &handshake, sizeof(int))) != sizeof(int)) {
				if (verbose>0)
						fprintf(stderr, "tcpsocket: could not write handshake, n = %d, should be %d\n", n, sizeof(int));
				goto cleanup;
		}

		if (!connect_accept)
				if (verbose>1)
						fprintf(stderr, "tcpsocket: dropping connection based on host details\n");
		if (!connect_continue)
				goto cleanup;

		/* read the job request options */
		if (message->job->optsize>0) {
				message->opt = malloc(message->job->optsize);
				if ((n = bufread(fd, message->opt, message->job->optsize)) != message->job->optsize) {
						if (verbose>0)
								fprintf(stderr, "tcpsocket: read size = %d, should be %d\n", n, message->job->optsize);
						goto cleanup;
				}
		}

		/* give a handshake */
		handshake = connect_accept | connect_continue;
		if ((n = bufwrite(fd, &handshake, sizeof(int))) != sizeof(int)) {
				if (verbose>0)
						fprintf(stderr, "tcpsocket: could not write handshake, n = %d, should be %d\n", n, sizeof(int));
				goto cleanup;
		}

		if (!connect_accept)
				if (verbose>1)
						fprintf(stderr, "tcpsocket: dropping connection based on host details\n");
		if (!connect_continue)
				goto cleanup;

		pthread_mutex_lock(&mutexjoblist);

		/* create a new list item */
		job = (joblist_t *)malloc(sizeof(joblist_t));
		job->job  = message->job;
		job->host = message->host;
		job->arg  = message->arg;
		job->opt  = message->opt;
		/* add the item to the beginning of the list */
		job->next = joblist;
		joblist = job;

		/* the message will be deallocated in the cleanup function */

		/* the responsibility for the dynamical memory of the message content has been
		   reassigned to the joblist and should not be deallocated here */

		if (verbose>2) {
				fprintf(stderr, "tcpsocket: job.version  = %d\n", job->job->version);
				fprintf(stderr, "tcpsocket: job.id       = %d\n", job->job->id);
				fprintf(stderr, "tcpsocket: job.argsize  = %d\n", job->job->argsize);
				fprintf(stderr, "tcpsocket: job.optsize  = %d\n", job->job->optsize);
				fprintf(stderr, "tcpsocket: host.name    = %s\n", job->host->name);
				fprintf(stderr, "tcpsocket: host.port    = %d\n", job->host->port);
				fprintf(stderr, "tcpsocket: host.id      = %d\n", job->host->id);
		}

		pthread_mutex_unlock(&mutexjoblist);

cleanup:
		printf(""); /* otherwise the pthread_cleanup_pop won't compile */
		pthread_cleanup_pop(1);
		return NULL;
}

