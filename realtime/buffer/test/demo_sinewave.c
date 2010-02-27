/*
 * Copyright (C) 2008, Christian Hesse & Robert Oostenveld
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>
#include "buffer.h"

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

    fprintf(stderr, "demo_sinewave: host.name =  %s\n", host.name);
    fprintf(stderr, "demo_sinewave: host.port =  %d\n", host.port);

	/* start the acquisition */
	sinewave_thread((void *)(&host));

	return 0;
}

