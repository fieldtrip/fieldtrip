/*
 * Copyright (C) 2010, Stefan Klanke
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "buffer.h"

int main(int argc, char *argv[]) {
	int server;
	message_t request, *response;
	messagedef_t requestdef;
	waitdef_t waitdef;

	if (argc != 6) {
		fprintf(stderr, "USAGE: application <server_ip> <port> <sample nr> <event nr> <ms>\n");
		exit(1);
	}

	/* open the TCP socket */
	if ((server = open_connection(argv[1], atoi(argv[2]))) < 0) {
		fprintf(stderr, "ERROR; failed to create socket\n");
		exit(1);
	}
	
	waitdef.threshold.nsamples = atoi(argv[3]);
	waitdef.threshold.nevents  = atoi(argv[4]);
	waitdef.milliseconds = atoi(argv[5]);

	requestdef.version = VERSION;
	requestdef.command = WAIT_DAT;
	requestdef.bufsize = sizeof(waitdef);
	
	request.def = &requestdef;
	request.buf = &waitdef;
	
	printf("Want to wait up to %i milliseconds for nsamples > %i || nevents > %i\n", waitdef.milliseconds, waitdef.threshold.nsamples, waitdef.threshold.nevents);
	
	printf("Sending request...\n");
	tcprequest(server, &request, &response);
	printf("Got it\n");
	
	if (response == NULL || response->def == NULL) {
		fprintf(stderr, "Unknown error in response\n");
	} else if (response->def->command != WAIT_OK) {
		fprintf(stderr, "Response != WAIT_OK : %X\n", response->def->command);
	} else if (response->def->bufsize != sizeof(samples_events_t) || response->buf == NULL) {
		fprintf(stderr, "Response->buf != samples_event_t\n");
	} else {
		samples_events_t *ns = (samples_events_t *) response->buf;
		printf("Number of samples: %i\n",ns->nsamples);
		printf("Number of events : %i\n",ns->nevents);
	}
	
	if (response!=NULL) {
		FREE(response->buf);
		FREE(response->def);
		free(response);
	}
		
	close(server);
	exit(0);
}

