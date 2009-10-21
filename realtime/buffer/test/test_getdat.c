/*
 * Copyright (C) 2008, Robert Oostenveld
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 * $Log: test_getdat.c,v $
 * Revision 1.7  2008/10/29 20:12:27  roboos
 * renamed open_remotehost into open_connection and added function close_connection
 * added some fprintf statements, no functional changes
 *
 * Revision 1.6  2008/07/08 20:25:13  roboos
 * print a separator line between the output sections
 *
 * Revision 1.5  2008/07/08 16:12:05  roboos
 * A lot of small changes related to cleaning up the release version. No functional changes in the actual core code.
 *
 * Revision 1.4  2008/03/26 14:38:59  thohar
 * portability again
 *
 * Revision 1.3  2008/02/26 21:43:26  roboos
 * renamed packet_t structure definition into messagedef_t, added message_t structure (contains def+buf)
 *
 * Revision 1.2  2008/02/20 13:48:49  roboos
 * specification of begsample and endsample on command line
 * correct printing of other datatypes
 *
 * Revision 1.1  2008/02/19 10:25:19  roboos
 * initial implementation, seems functional
 *
 *
 */

#ifdef WIN32
	#include <winsock2.h>
	//#include <win32compat.h>
#else
	#include <netdb.h>
#endif
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
	write(server, &request, sizeof(messagedef_t));
	if (argc>4)
		write(server, &datasel, sizeof(datasel_t));

	read(server, &response, sizeof(messagedef_t));
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
			data.buf = buf+sizeof(datadef_t);
			print_datadef(data.def);
			for (j=0; j<data.def->nsamples; j++) {
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

