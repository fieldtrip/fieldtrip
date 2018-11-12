/*
 * Copyright (C) 2008, Robert Oostenveld & Christian Hesse
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 */

#include <stdio.h>
#include <stdlib.h>

#include "message.h"

void print_request(messagedef_t *request) {
	fprintf(stderr, "request.version = 0x%04x\n", request->version);
	fprintf(stderr, "request.command = 0x%04x\n", request->command);
	fprintf(stderr, "request.bufsize = %u\n",     request->bufsize);
}

void print_response(messagedef_t *response) {
	fprintf(stderr, "response.version = 0x%04x\n", response->version);
	fprintf(stderr, "response.command = 0x%04x\n", response->command);
	fprintf(stderr, "response.bufsize = %u\n",     response->bufsize);
}

void print_headerdef(headerdef_t *headerdef) {
	if (headerdef==NULL)
		fprintf(stderr, "headerdef==NULL\n");
	else {
		fprintf(stderr, "headerdef.nchans    = %u\n", headerdef->nchans);
		fprintf(stderr, "headerdef.nsamples  = %u\n", headerdef->nsamples);
		fprintf(stderr, "headerdef.nevents   = %u\n", headerdef->nevents);
		fprintf(stderr, "headerdef.fsample   = %f\n", headerdef->fsample);
		fprintf(stderr, "headerdef.data_type = %u\n", headerdef->data_type);
		fprintf(stderr, "headerdef.bufsize   = %u\n", headerdef->bufsize);
	}
}

void print_datadef(datadef_t *datadef) {
	if (datadef==NULL)
		fprintf(stderr, "datadef==NULL\n");
	else {
		fprintf(stderr, "datadef.nchans    = %u\n", datadef->nchans);
		fprintf(stderr, "datadef.nsamples  = %u\n", datadef->nsamples);
		fprintf(stderr, "datadef.data_type = %u\n", datadef->data_type);
		fprintf(stderr, "datadef.bufsize   = %u\n", datadef->bufsize);
	}
}

void print_eventdef(eventdef_t *eventdef) {
	if (eventdef==NULL)
		fprintf(stderr, "eventdef==NULL\n");
	else {
		fprintf(stderr, "eventdef.type_type       = %u\n", eventdef->type_type      );
		fprintf(stderr, "eventdef.type_numel      = %u\n", eventdef->type_numel     );
		fprintf(stderr, "eventdef.value_type      = %u\n", eventdef->value_type     );
		fprintf(stderr, "eventdef.value_numel     = %u\n", eventdef->value_numel    );
		fprintf(stderr, "eventdef.sample          = %d\n", eventdef->sample         );
		fprintf(stderr, "eventdef.offset          = %d\n", eventdef->offset         );
		fprintf(stderr, "eventdef.duration        = %d\n", eventdef->duration       );
		fprintf(stderr, "eventdef.bufsize         = %u\n", eventdef->bufsize        );
	}
}

void print_datasel(datasel_t *datasel) {
	fprintf(stderr, "datasel.begsample  = %u\n", datasel->begsample);
	fprintf(stderr, "datasel.endsample  = %u\n", datasel->endsample);
}

void print_eventsel(eventsel_t *eventsel) {
	fprintf(stderr, "eventsel.begevent  = %u\n", eventsel->begevent);
	fprintf(stderr, "eventsel.endevent  = %u\n", eventsel->endevent);
}

void print_buf(void *buf, int bufsize) {
	int i;
	fprintf(stderr, "buf     =");
	if (buf==NULL)
		fprintf(stderr, " NULL");
	else
		for (i=0; i<bufsize; i++)
			fprintf(stderr, " %02x", ((unsigned char*)buf)[i]);
	/* FIXME: fprint expects uint. */
	fprintf(stderr, "\n");
}
