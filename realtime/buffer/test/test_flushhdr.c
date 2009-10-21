/*
 * Copyright (C) 2008, Robert Oostenveld
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 * $Log: test_flushhdr.c,v $
 * Revision 1.4  2008/10/29 20:12:27  roboos
 * renamed open_remotehost into open_connection and added function close_connection
 * added some fprintf statements, no functional changes
 *
 * Revision 1.3  2008/07/08 20:25:13  roboos
 * print a separator line between the output sections
 *
 * Revision 1.2  2008/07/08 16:12:05  roboos
 * A lot of small changes related to cleaning up the release version. No functional changes in the actual core code.
 *
 * Revision 1.1  2008/05/29 13:47:31  roboos
 * added command line executables for flushing the buffer
 *
 */

#include <netdb.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <unistd.h>
#include <pthread.h>
#include <math.h>

#include "buffer.h"
#include "message.h"

int main(int argc, char *argv[]) {
	int server;
	int n;
	messagedef_t request, response;
	void *buf = NULL;

	if (argc != 3) {
		fprintf(stderr, "USAGE: application <server_ip> <port>\n");
		exit(1);
	}

	/* open the TCP socket */
	if ((server = open_connection(argv[1], atoi(argv[2]))) < 0) {
		fprintf(stderr, "ERROR; failed to create socket\n");
		exit(1);
	}

	request.version = VERSION;
	request.command = FLUSH_HDR;
	request.bufsize = 0;

	fprintf(stderr, "------------------------------\n");
	print_request(&request);
	write(server, &request, sizeof(messagedef_t));

	read(server, &response, sizeof(messagedef_t));
	fprintf(stderr, "------------------------------\n");
	print_response(&response);
	fprintf(stderr, "------------------------------\n");

	if (response.command==GET_OK) {
		buf = malloc(response.bufsize);
		if ((n = bufread(server, buf, response.bufsize)) < response.bufsize) {
			fprintf(stderr, "problem reading enough bytes (%d)\n", n);
		}
		FREE(buf);
	}

	close(server);
	exit(0);
}

