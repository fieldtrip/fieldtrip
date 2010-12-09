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
#include "socketserver.h"
#include <signal.h>

ft_buffer_server_t *S;
volatile int keepRunning = 1;
double timePutHeader = 0.0;

int my_request_handler(const message_t *request, message_t **response, void *user_data) {
	struct timeval tv;
	double tAbs, tRel;
	int res;
	
	gettimeofday(&tv, NULL);
	tAbs = tv.tv_sec + tv.tv_usec*1e-6;
	tRel = tAbs - timePutHeader;
	
	printf("t=%8.3f ", tRel);
	switch(request->def->command) {
		case PUT_HDR:
			printf("Put header, bufsize = %i ... ", request->def->bufsize);
			timePutHeader = tAbs;
			break;
		case GET_HDR:
			printf("Get header ...");
			break;
		case FLUSH_HDR:
			printf("Flush header ...");
			break;
		case PUT_DAT:
			{
				const datadef_t *ddef = (const datadef_t *) request->buf;
				printf("Put data, %i channels, %i samples, type=%i, size=%i ... ", ddef->nchans, ddef->nsamples, ddef->data_type, ddef->bufsize);
			}
			break;
		case GET_DAT:
			if (request->def->bufsize >= sizeof(datasel_t)) {
				const datasel_t *ds = (const datasel_t *) request->buf;
				printf("Get data, start=%i, end=%i ... ", ds->begsample, ds->endsample);
			} else {
				printf("Get all data ... ");
			}
			break;
		case FLUSH_DAT:
			printf("Flush data ... ");
			break;
		case PUT_EVT:
			printf("Put events, bufsize = %i ... ", request->def->bufsize);
			break;
		case GET_EVT:
			if (request->def->bufsize >= sizeof(eventsel_t)) {
				const eventsel_t *es = (const eventsel_t *) request->buf;
				printf("Get events, start=%i, end=%i ... ", es->begevent, es->endevent);
			} else {
				printf("Get all events ... ");
			}
			break;
		case FLUSH_EVT:
			printf("Flush events ... ");
			break;
		case WAIT_DAT:
			if (request->def->bufsize >= sizeof(waitdef_t)) {
				const waitdef_t *wd = (const waitdef_t *) request->buf;
				printf("Wait data, nsamples=%i, nevents=%i, timeout=%i ... \n", wd->threshold.nsamples, wd->threshold.nevents, wd->milliseconds);
			} else {
				printf("Wait data, malformed! ... \n");
			}
			break;
	}
	
	res = dmarequest(request, response);
	if (res != 0) {
		printf("ERROR\n");
	} else {
		switch((*response)->def->command) {
			case WAIT_OK:
				gettimeofday(&tv, NULL);
				tAbs = tv.tv_sec + tv.tv_usec*1e-6;
				tRel = tAbs - timePutHeader;
				printf("t=%8.3f WAIT_OK\n", tRel);
				break;
			case WAIT_ERR:				
				printf("WAIT_ERR\n");
				break;			
			case PUT_OK:
			case GET_OK:
			case FLUSH_OK:
				printf("OK\n");
				break;
			case PUT_ERR:
			case GET_ERR:
			case FLUSH_ERR:
				printf("FAILED\n");
				break;
			default:
				printf("UNRECOGNIZED\n");
		}
	}
	return res;
}

void abortHandler(int sig) {
	keepRunning = 0;
}

int main(int argc, char *argv[]) {
	int port;
	char *name = NULL;
	
    /* verify that all datatypes have the expected syze in bytes */
    check_datatypes();

	if (argc>1) {
		port = atoi(argv[1]);
		if (port == 0) {
			name = argv[1];
		}	
	} else {
		port = 1972;
	}
	
	/* with enabled debug output */
	S = ft_start_buffer_server(port, name, my_request_handler, NULL);
	/* plain server without extra output */
	/*
	S = ft_start_buffer_server(port, name, NULL, NULL);
	*/
	if (S==NULL) return 1;
	signal(SIGINT, abortHandler);
	while (keepRunning) {
		usleep(1000000);
	}
	printf("Ctrl-C pressed -- stopping buffer server...\n");
	ft_stop_buffer_server(S);
	printf("Done.\n");
	return 0;
}

