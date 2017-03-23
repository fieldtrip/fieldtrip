/*
 * demo - a multithreaded and networked buffer for EEG/MEG data and events
 *
 * Copyright (C) 2008, Christian Hesse & Robert Oostenveld
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "buffer.h"
#include <pthread.h>

/***********************************************************************
  this is the main thread
 ***********************************************************************/
int main(int argc, char *argv[]) {
	host_t host;

	/* these variables are for the threading */
	int rc;
	pthread_t tid;

  check_datatypes();

	if (argc>2) {
		sprintf(host.name, argv[1]);
		host.port = atoi(argv[2]);
	}
	else {
		sprintf(host.name, DEFAULT_HOSTNAME);
		host.port = DEFAULT_PORT;
	}

	/* start the buffer in a seperate thread */
	rc = pthread_create(&tid, NULL, tcpserver, (void *)(&host));
	if (rc) {
		fprintf(stderr, "Error: return code from pthread_create() is %d\n", rc);
		exit(-1);
	}
  pthread_detach(tid);
	usleep(1000000);

	/* start the acquisition */
    sinewave_thread((void *)(&host));

	return 0;
}

