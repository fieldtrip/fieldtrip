/*
 * This piece of code demonstrates how an acquisition device could
 * be writing events (e.g. TTL triggers or battery status information)
 * to the buffer.
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

void *event_thread(void *arg) {
	long tdif;
	host_t *host = (host_t *)arg;

	/* these will be sent as event */
    int sample  = 0;
    int value   = 0;
    char *type  = "trigger";

	/* these are used in the communication and represent statefull information */
	int server             = -1;
	message_t    *request  = NULL;
	message_t    *response = NULL;
	header_t     *header   = NULL;
	data_t       *data     = NULL;
	event_t      *event    = NULL;

	/* these represent the acquisition system properties */
	int fsample        = 250;
	int blocksize      = 125;

	/* this is to prevent closing the thread at an unwanted moment and memory from leaking */
	int oldcancelstate, oldcanceltype;
	pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, &oldcancelstate);
	pthread_setcanceltype(PTHREAD_CANCEL_DEFERRED, &oldcanceltype);
	pthread_cleanup_push(cleanup_message, request);
	pthread_cleanup_push(cleanup_message, response);
	pthread_cleanup_push(cleanup_header,  header);
	pthread_cleanup_push(cleanup_data,    data);
	pthread_cleanup_push(cleanup_event,   event);
	pthread_cleanup_push(cleanup_socket,  &server);

	/* this determines the hostname and port on which the client will write */
	if (!arg)
		exit(1);

	fprintf(stderr, "event: host.name =  %s\n", host->name);
	fprintf(stderr, "event: host.port =  %d\n", host->port);

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

	while (1) {

		/* increment the trigger value on each iteration */
		value++;
		sample += blocksize;

		/* construct the event definition */
		event->def->type_type   = DATATYPE_CHAR; 
		event->def->type_numel  = strlen(type);
		event->def->value_type  = DATATYPE_INT32;
		event->def->value_numel = 1;
		event->def->sample      = sample;
		event->def->offset      = 0;
		event->def->duration    = 0;
		event->def->bufsize     = 0; /* see below */
		event->buf              = NULL;

		/* add the variable information to the event buffer */
		event->def->bufsize = append(&event->buf, event->def->bufsize, type,   event->def->type_numel);	/* type is a char-array */
		event->def->bufsize = append(&event->buf, event->def->bufsize, &value, WORDSIZE_INT32);			/* value is an integer  */

		/* construct the request */
		request->def->command = PUT_EVT;
		request->def->bufsize = append(&request->buf, request->def->bufsize, event->def, sizeof(eventdef_t));
		request->def->bufsize = append(&request->buf, request->def->bufsize, event->buf, event->def->bufsize);

		server = open_connection(host->name, host->port);
		clientrequest(server, request, &response);
		if (server>=0)
			closesocket(server);

		/* FIXME do someting with the response, i.e. check that it is OK */
		request->def->bufsize = 0;
		FREE(request->buf);
		event->def->bufsize = 0;
		FREE(event->buf);
		if (response) {
			FREE(response->def);
			FREE(response->buf);
			FREE(response);
		}

		/* approximate delay in microseconds */
		tdif = (long)(blocksize * 1000000 / fsample);
		usleep(tdif);

	} /* while(1) */

cleanup:
	/* from now on it is save to cancel the thread */
	pthread_setcancelstate(oldcancelstate, NULL);
	pthread_setcanceltype(oldcanceltype, NULL);

	pthread_cleanup_pop(1);  /* server */
	pthread_cleanup_pop(1);  /* event */
	pthread_cleanup_pop(1);  /* data */
	pthread_cleanup_pop(1);  /* header */
	pthread_cleanup_pop(1);  /* response */
	pthread_cleanup_pop(1);  /* request */

	pthread_exit(NULL);
	return NULL;
}

