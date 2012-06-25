/*
 * Copyright (C) 2008, Robert Oostenveld
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/* #include <unistd.h> */
#include <math.h>
#include "buffer.h"

int main(int argc, char *argv[]) {
	int server;
	int i, j, k, n;
	messagedef_t request, response;
	data_t    data;
	datasel_t datasel;
	void *buf = NULL;

	if (argc < 3) {
		fprintf(stderr, "USAGE: application <server_ip> <port>\n");
		exit(1);
	}

	/* open the TCP socket */
	if ((server = open_connection(argv[1], atoi(argv[2]))) < 0) {
		fprintf(stderr, "ERROR; failed to create socket\n");
		exit(1);
	}

	request.version = VERSION;
	request.command = GET_DAT;

	if (argc>4) {
		request.bufsize = sizeof(datasel_t);
		datasel.begsample = atoi(argv[3]);
		datasel.endsample = atoi(argv[4]);
	}
	else {
		request.bufsize = 0;
		datasel.begsample = 0;
		datasel.endsample = 0;
	}

	fprintf(stderr, "------------------------------\n");
	print_request(&request);
	bufwrite(server, &request, sizeof(messagedef_t));
	if (argc>4)
		bufwrite(server, &datasel, sizeof(datasel_t));

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
			data.def = (datadef_t *)buf;
			data.buf = (char *) buf+sizeof(datadef_t);
			print_datadef(data.def);
			if (data.def->nsamples * data.def->nchans > 200) {
				printf("Lots of data - not printing...\n");
			} 
			else for (j=0; j<data.def->nsamples; j++) {
				for (i=0; i<data.def->nchans; i++) {
					//k = i*data.def->nsamples + j;
					k = j*data.def->nchans + i;
					switch (data.def->data_type) {
						case DATATYPE_INT8:
							fprintf(stdout, "%9d ", ((INT8_T*)data.buf)[k]);
							break;
						case DATATYPE_INT16:
							fprintf(stdout, "%9d ", ((INT8_T*)data.buf)[k]);
							break;
						case DATATYPE_INT32:
							fprintf(stdout, "%9d ", ((INT8_T*)data.buf)[k]);
							break;
						case DATATYPE_INT64:
							fprintf(stdout, "%9d ", ((INT8_T*)data.buf)[k]);
							break;
						case DATATYPE_FLOAT32:
							fprintf(stdout, "%9.3f ", ((FLOAT32_T*)data.buf)[k]);
							break;
						case DATATYPE_FLOAT64:
							fprintf(stdout, "%9.3f ", ((FLOAT64_T*)data.buf)[k]);
							break;
					}
				}
				fprintf(stdout, "\n");
			}
		}
		FREE(buf);
	}

	close(server);
	exit(0);
}

