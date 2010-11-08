
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <getopt.h>
#include <string.h>

#include "engine.h"
#include "matrix.h"

#include "peer.h"
#include "extern.h"
#include "platform_includes.h"

int main(int argc, char *argv[]) {
		int rc;
		joblist_t *job;

		/* the thread IDs are needed for cancelation at cleanup */
		pthread_t udsserverThread;
		pthread_t tcpserverThread;
		pthread_t announceThread;
		pthread_t discoverThread;
		pthread_t expireThread;

		openlog("peerslave", LOG_PID | LOG_PERROR, LOG_USER);
		setlogmask(LOG_MASK(LOG_EMERG) | LOG_MASK(LOG_ALERT) | LOG_MASK(LOG_CRIT));
		peerinit(NULL);

		if ((rc = pthread_create(&tcpserverThread, NULL, tcpserver, (void *)NULL))>0) {
				PANIC("failed to start tcpserver thread\n");
		}
		else {
				DEBUG(LOG_NOTICE, "started tcpserver thread");
		}

		if ((rc = pthread_create(&announceThread, NULL, announce, (void *)NULL))>0) {
				PANIC("failed to start announce thread\n");
		}
		else {
				DEBUG(LOG_NOTICE, "started announce thread");
		}

		if ((rc = pthread_create(&discoverThread, NULL, discover, (void *)NULL))>0) {
				PANIC("failed to start discover thread\n");
		}
		else {
				DEBUG(LOG_NOTICE, "started discover thread");
		}

		if ((rc = pthread_create(&expireThread, NULL, expire, (void *)NULL))>0) {
				DEBUG(LOG_NOTICE, "failed to start expire thread");
				exit(1);
		}
		else {
				DEBUG(LOG_NOTICE, "started expire thread");
		}

		pthread_mutex_lock(&mutexhost);
		/* switch the peer to idle slave */
		host->status = STATUS_IDLE;
		pthread_mutex_unlock(&mutexhost);

		while (1) {
				fprintf(stderr, "test: --------------------------\n");
				pthread_mutex_lock(&mutexjoblist);
				if (joblist) {
						fprintf(stderr, "test: job.version  = %u\n", joblist->job->version);
						fprintf(stderr, "test: job.id       = %u\n", joblist->job->id);
						fprintf(stderr, "test: job.argsize  = %u\n", joblist->job->argsize);
						fprintf(stderr, "test: job.optsize  = %u\n", joblist->job->optsize);
						fprintf(stderr, "test: host.name    = %s\n", joblist->host->name);
						fprintf(stderr, "test: host.port    = %u\n", joblist->host->port);
						fprintf(stderr, "test: host.id      = %u\n", joblist->host->id);
				}
				pthread_mutex_unlock(&mutexjoblist);
				usleep(1000000);
		}

		return 0;
}
