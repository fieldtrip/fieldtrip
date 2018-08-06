/*
 * Copyright (C) 2008, Christian Hesse & Robert Oostenveld
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "buffer.h"

int main(int argc, char *argv[]) {
	host_t host;

    /* verify that all datatypes have the expected syze in bytes */
    check_datatypes();

	sprintf(host.name, DEFAULT_HOSTNAME);
	if (argc>1) {
		host.port = atoi(argv[1]);
	}
	else {
	    printf("Using default port, recommended usage 'buffer [port]'. \n");
		host.port = DEFAULT_PORT;
	}

	/* start the buffer */
	printf("Starting FieldTrip buffer on port %d... \n", host.port);
	tcpserver((void *)(&host));

	return 0;
}

