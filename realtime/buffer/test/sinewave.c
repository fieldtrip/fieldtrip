/*
 * This piece of code demonstrates how an acquisition device could be
 * writing header information (e.g. number of channels and sampling
 * frequency) and streaming EEG data to the buffer.
 *
 * Copyright (C) 2008, Robert Oostenveld
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "buffer.h"
#include <pthread.h>

#define FREQ       1
#define PI         3.1415926

void *sinewave_thread(void *arg) {
	int i, j, k, sample = 0, status = 0, verbose = 0;
	long tdif;
	host_t *host = (host_t *)arg;

	/* these represent the acquisition system properties */
	int nchans         = 16;
	int fsample        = 250;
	int blocksize      = 125;

	/* these are used in the communication and represent statefull information */
	int server             = -1;
	message_t    *request  = NULL;
	message_t    *response = NULL;
	header_t     *header   = NULL;
	data_t       *data     = NULL;
	event_t      *event    = NULL;

	/* this determines the hostname and port on which the client will write */
	if (!arg)
		exit(1);

	if (verbose>0) fprintf(stderr, "sinewave: host.name =  %s\n", host->name);
	if (verbose>0) fprintf(stderr, "sinewave: host.port =  %d\n", host->port);

	/* allocate the elements that will be used in the communication */
	request      = malloc(sizeof(message_t));
	request->def = malloc(sizeof(messagedef_t));
	request->buf = NULL;
	request->def->version = VERSION;
	request->def->bufsize = 0;

	header      = malloc(sizeof(header_t));
	header->def = malloc(sizeof(headerdef_t));
	header->buf = NULL;

	data      = malloc(sizeof(data_t));
	data->def = malloc(sizeof(datadef_t));
	data->buf = NULL;

	event      = malloc(sizeof(event_t));
	event->def = malloc(sizeof(eventdef_t));
	event->buf = NULL;

	/* define the header */
	header->def->nchans    = nchans;
	header->def->nsamples  = 0;
	header->def->nevents   = 0;
	header->def->fsample   = fsample;
	header->def->data_type = DATATYPE_FLOAT32;
	header->def->bufsize   = 0;
	FREE(header->buf);

	/* define the constant part of the data and allocate space for the variable part */
	data->def->nchans    = nchans;
	data->def->nsamples  = blocksize;
	data->def->data_type = DATATYPE_FLOAT32;
	data->def->bufsize   = WORDSIZE_FLOAT32*nchans*blocksize;
	FREE(data->buf);
	data->buf            = malloc(WORDSIZE_FLOAT32*nchans*blocksize);

	/* initialization phase, send the header */
	request->def->command = PUT_HDR;
	request->def->bufsize = append(&request->buf, request->def->bufsize, header->def, sizeof(headerdef_t));
	request->def->bufsize = append(&request->buf, request->def->bufsize, header->buf, header->def->bufsize);

	server = open_connection(host->name, host->port);
	status = clientrequest(server, request, &response);
    if (verbose>0) fprintf(stderr, "sinewave: clientrequest returned %d\n", status);
    if (status) {
      fprintf(stderr, "sinewave: err1\n");
      goto cleanup;
    }
	status = close_connection(server);
    if (status) {
      fprintf(stderr, "sinewave: err2\n");
      goto cleanup;
    }

	/* FIXME do someting with the response, i.e. check that it is OK */
	cleanup_message(&request);
	cleanup_message(&response);
    request = NULL;
    response = NULL;

	while (1) {
	//for (k=0; k<10; k++) {

		for (j=0; j<blocksize; j++) {
			if (nchans>0)
				((FLOAT32_T *)(data->buf))[j*nchans+0] = sample;
			if (nchans>1)
				((FLOAT32_T *)(data->buf))[j*nchans+1] = sin(2*PI*FREQ*sample/fsample);
			if (nchans>2)
				((FLOAT32_T *)(data->buf))[j*nchans+2] = 2.0*((FLOAT32_T)rand())/RAND_MAX - 1.0;
			if (nchans>3)
				for (i=3; i<nchans; i++)
					((FLOAT32_T *)(data->buf))[j*nchans+i] = sin(2*PI*FREQ*sample/fsample) + 2.0*((FLOAT32_T)rand())/RAND_MAX - 1.0;
			sample++;
		}

		/* create the request */
		request      = malloc(sizeof(message_t));
		request->def = malloc(sizeof(messagedef_t));
		request->buf = NULL;
		request->def->version = VERSION;
		request->def->bufsize = 0;
		request->def->command = PUT_DAT;
		request->def->bufsize = append(&request->buf, request->def->bufsize, data->def, sizeof(datadef_t));
		request->def->bufsize = append(&request->buf, request->def->bufsize, data->buf, data->def->bufsize);

		server = open_connection(host->name, host->port);
		status = clientrequest(server, request, &response);
        if (verbose>0) fprintf(stderr, "sinewave: clientrequest returned %d\n", status);
    	if (status) {
			fprintf(stderr, "sinewave: err3\n");
      		goto cleanup;
		}
		status = close_connection(server);
    	if (status) {
      		fprintf(stderr, "sinewave: err4\n");
     		goto cleanup;
		}

		/* FIXME do someting with the response, i.e. check that it is OK */
		cleanup_message(&request);
		cleanup_message(&response);
    	request = NULL;
    	response = NULL;

		/* approximate delay in microseconds */
		tdif = (long)(blocksize * 1000000 / fsample);
		usleep(tdif);

	} /* while(1) */

cleanup:
	cleanup_event(&event);
	cleanup_data(&data);
	cleanup_header(&header);
	cleanup_message(&request);
	cleanup_message(&response);

    pthread_exit(0);
	return 0;
}

