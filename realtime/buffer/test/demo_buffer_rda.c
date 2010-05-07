/*
 * Copyright (C) 2008, Christian Hesse & Robert Oostenveld
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <time.h>
#include "buffer.h"
#include "rdaserver.h"

rda_server_ctrl_t *rdac = NULL;

void abortHandler(int sig) {
	printf("Ctrl-C pressed -- stopping RDA server...\n");
	rda_stop_server(rdac);
	printf("Done.\n");
	/* TODO: do the same with the tcp server if run from a thread */
	exit(0);
}

int main(int argc, char *argv[]) {
	host_t host;
	sigset_t sigInt;
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
	
	/* Make sure the RDA server thread does not receive CTRL-C, instead it will
		get terminated properly using rda_stop_server from a handler. We need this
		because otherwise it is unspecified in which thread the signal will be
		received - the thread created in rda_start_server will inherit the signal
		mask we temporarily modify here.
	*/
	sigemptyset(&sigInt);
	sigaddset(&sigInt, SIGINT);
	sigprocmask(SIG_BLOCK, &sigInt, NULL);
	
	/* 0,0,0 = dma-connection to fieldtrip, default = float, default port */
	rdac = rda_start_server(0, 0, 0, &errval);
	if (errval != 0) {
		printf("Could not start rdaserver: %i\n", errval);
		return errval;
	}
	
	/* We want CTRL-C in this thread */
	sigprocmask(SIG_UNBLOCK, &sigInt, NULL);
	signal(SIGINT, abortHandler);
	
	/* start the buffer */
	tcpserver((void *)(&host));

	return 0;
}

