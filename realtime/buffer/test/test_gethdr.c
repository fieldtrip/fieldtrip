/*
 * Copyright (C) 2008, Robert Oostenveld
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
	int n;
	messagedef_t request, response;
	header_t header;
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
	request.command = GET_HDR;
	request.bufsize = 0;

	fprintf(stderr, "------------------------------\n");
	print_request(&request);
	bufwrite(server, &request, sizeof(messagedef_t));

	bufread(server, &response, sizeof(messagedef_t));
	fprintf(stderr, "------------------------------\n");
	print_response(&response);
	fprintf(stderr, "------------------------------\n");

	if (response.command==GET_OK) {
		buf = malloc(response.bufsize);
		if ((n = bufread(server, buf, response.bufsize)) < response.bufsize) {
			fprintf(stderr, "problem reading enough bytes (%d)\n", n);
		}
		else {
			header.def = (headerdef_t *)buf;
			header.buf = (char *) buf+sizeof(headerdef_t);
			print_headerdef(header.def);
		}
		FREE(buf);
	}

	close(server);
	exit(0);
}

