#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "peer.h"
#include "extern.h"
#include "socket_includes.h"
#include "unix_includes.h"

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
		int handshake = 0;
		joblist_t *job;

		threadlocal_t threadlocal;
		threadlocal.message = NULL;
		threadlocal.fd = -1;

		/* these are used for communication over the TCP socket */
		int fd = 0;
		message_t *message = NULL;

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

		if (hoststatus()==1 && jobcount()==0)
				handshake = 1;
		else if (hoststatus()==2)
				handshake = 1;
		else
				handshake = 0;

		/* give a handshake */
		if ((n = bufwrite(fd, &handshake, sizeof(int))) != sizeof(int)) {
				if (verbose>0) fprintf(stderr, "tcpsocket: could not write handshake, n = %d, should be %d\n", n, sizeof(int));
				goto cleanup;
		}
		else if (!handshake) {
				if (verbose>1)
						fprintf(stderr, "tcpsocket: dropping connection based on status\n");
				goto cleanup;
		}

		message       = (message_t*)malloc(sizeof(message_t));
		message->host = (hostdef_t*)malloc(sizeof(hostdef_t));
		message->job  = (jobdef_t*)malloc(sizeof(jobdef_t));
		message->arg  = NULL;
		message->opt  = NULL;

		/* this will be deallocated at cleanup */
		threadlocal.message = message;

		/* read the host details */
		if ((n = bufread(fd, message->host, sizeof(hostdef_t))) != sizeof(hostdef_t)) {
				if (verbose>0) fprintf(stderr, "tcpsocket: read size = %d, should be %d\n", n, sizeof(hostdef_t));
				goto cleanup;
		}

		if (verbose>1) {
				fprintf(stderr, "tcpsocket: ismember_userlist  = %d\n", ismember_userlist(message->host->user));
				fprintf(stderr, "tcpsocket: ismember_grouplist = %d\n", ismember_grouplist(message->host->group));
				fprintf(stderr, "tcpsocket: ismember_hostlist  = %d\n", ismember_hostlist(message->host->name));
		}

		/* test whether the request can be accepted based on the host characteristics */
		handshake = 1;
		handshake = (handshake & ismember_userlist (message->host->user));
		handshake = (handshake & ismember_grouplist(message->host->group));
		handshake = (handshake & ismember_hostlist (message->host->name));

		/* give a handshake */
		if ((n = bufwrite(fd, &handshake, sizeof(int))) != sizeof(int)) {
				if (verbose>0) fprintf(stderr, "tcpsocket: could not write handshake, n = %d, should be %d\n", n, sizeof(int));
				goto cleanup;
		}
		else if (!handshake) {
				if (verbose>1)
						fprintf(stderr, "tcpsocket: dropping connection based on host details\n");
				goto cleanup;
		}

		/* read the job details */
		if ((n = bufread(fd, message->job, sizeof(jobdef_t))) != sizeof(jobdef_t)) {
				if (verbose>0) fprintf(stderr, "tcpsocket: packet size = %d, should be %d\n", n, sizeof(jobdef_t));
				goto cleanup;
		}

		/* test whether the request can be accepted based on the job characteristics */
		handshake = 1;

		if (message->job->version!=VERSION) {
				if (verbose>0) fprintf(stderr, "tcpsocket: incorrect version\n");
				handshake = 0;
		}

		pthread_mutex_lock(&mutexhost);
		if (message->job->memreq > host->memavail) {
				if (verbose>0) fprintf(stderr, "tcpsocket: memory request too large\n");
				handshake = 0;
		}
		if (message->job->cpureq > host->cpuavail) {
				if (verbose>0) fprintf(stderr, "tcpsocket: cpu request too large\n");
				handshake = 0;
		}
		if (message->job->timreq > host->timavail) {
				if (verbose>0) fprintf(stderr, "tcpsocket: time request too large\n");
				handshake = 0;
		}
		pthread_mutex_unlock(&mutexhost);

		if (message->job->argsize>MAXARGSIZE) {
				if (verbose>0) fprintf(stderr, "tcpsocket: argsize too large\n");
				handshake = 0;
		}

		if (message->job->optsize>MAXARGSIZE) {
				if (verbose>0) fprintf(stderr, "tcpsocket: optsize too large\n");
				handshake = 0;
		}

		/* give a handshake */
		if ((n = bufwrite(fd, &handshake, sizeof(int))) != sizeof(int)) {
				if (verbose>0) fprintf(stderr, "tcpsocket: could not write handshake, n = %d, should be %d\n", n, sizeof(int));
				goto cleanup;
		}
		else if (!handshake) {
				if (verbose>1)
						fprintf(stderr, "tcpsocket: dropping connection based on job details\n");
				goto cleanup;
		}

		/* read the job request arguments */
		if (message->job->argsize>0) {
				message->arg = malloc(message->job->argsize);
				if ((n = bufread(fd, message->arg, message->job->argsize)) != message->job->argsize) {
						if (verbose>0) fprintf(stderr, "tcpsocket: read size = %d, should be %d\n", n, message->job->argsize);
						goto cleanup;
				}
		}

		/* give a handshake */
		if ((n = bufwrite(fd, &handshake, sizeof(int))) != sizeof(int)) {
				if (verbose>0) fprintf(stderr, "tcpsocket: could not write handshake, n = %d, should be %d\n", n, sizeof(int));
				goto cleanup;
		}
		else if (!handshake) {
				if (verbose>1)
						fprintf(stderr, "tcpsocket: dropping connection based on host details\n");
				goto cleanup;
		}

		/* read the job request options */
		if (message->job->optsize>0) {
				message->opt = malloc(message->job->optsize);
				if ((n = bufread(fd, message->opt, message->job->optsize)) != message->job->optsize) {
						if (verbose>0) fprintf(stderr, "tcpsocket: read size = %d, should be %d\n", n, message->job->optsize);
						goto cleanup;
				}
		}

		/* give a handshake */
		if ((n = bufwrite(fd, &handshake, sizeof(int))) != sizeof(int)) {
				if (verbose>0) fprintf(stderr, "tcpsocket: could not write handshake, n = %d, should be %d\n", n, sizeof(int));
				goto cleanup;
		}
		else if (!handshake) {
				if (verbose>1)
						fprintf(stderr, "tcpsocket: dropping connection based on host details\n");
				goto cleanup;
		}

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
				fprintf(stderr, "tcpsocket: host.addr    = %s\n", job->host->addr);
				fprintf(stderr, "tcpsocket: host.port    = %d\n", job->host->port);
		}

		pthread_mutex_unlock(&mutexjoblist);

cleanup:
		printf(""); /* otherwise the pthread_cleanup_pop won't compile */
		pthread_cleanup_pop(1);
		return NULL;
}

