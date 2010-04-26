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
#include "rdaserver.h"

int main(int argc, char *argv[]) {
	host_t host;
	rda_server_ctrl_t *rdac;
	int errval;

    /* verify that all datatypes have the expected syze in bytes */
    check_datatypes();

	if (argc>2) {
		sprintf(host.name, argv[1]);
		host.port = atoi(argv[2]);
	}
	else {
		sprintf(host.name, DEFAULT_HOSTNAME);
		host.port = DEFAULT_PORT;
	}
	
	/* 0,0,0 = dma-connection to fieldtrip, default = float, default port */
	rdac = rda_start_server(0, 0, 0, &errval);
	if (errval != 0) {
		printf("Could not start rdaserver: %i\n", errval);
		return errval;
	}
	
	/* start the buffer */
	tcpserver((void *)(&host));

	return 0;
}

