/*
 * Simple application that sets up a buffer and writes all incoming messages (with timestamps) to disk for detailled debugging.
 * The data stream can be played back with the same timing using the 'playback' application.
 *
 * Copyright (C) 2010, Stefan Klanke
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <signal.h>

#include "buffer.h"
#include "socketserver.h"
#include <ft_storage.h>

#if !defined(PLATFORM_WINDOWS) && defined(COMPILER_MINGW)
#include <direct.h>
#else
#include <sys/stat.h>
#include <sys/time.h>
#endif

#define QUEUE_SIZE  100

char content_descr[]=
"This directory contains FieldTrip buffer data in V1 format.\n" \
"\n" \
"The numbered subdirectories correspond to a complete cycle of\n" \
"filling the buffer by first writing the header, and then\n" \
"writing samples and events. All of these quantities are stored\n" \
"in the same binary form as in the buffer itself. Additionally,\n" \
"a plain text file 'header.txt' is written to allow manual inspection.\n" \
"Finally, the file 'timing' is a plain text file that describes how the\n" \
"buffer was filled over time. In that file, a line\n" \
"S 200 0.030\n" \
"means that 0.03 seconds after writing the header, a block of 200 samples\n" \
"has been written. Similarly,\n" \
"E 2 0.124\n" \
"means that 124ms after writing the header, 2 events were written.\n";

typedef struct {
	volatile int command;
	int quantity;
	double t;
} QueueElement;

QueueElement queue[QUEUE_SIZE];
int qWritePos = 0, qReadPos = 0;
pthread_mutex_t qMutex = PTHREAD_MUTEX_INITIALIZER;

ft_buffer_server_t *S;
ft_storage_t *OS = NULL;
volatile int keepRunning = 1;
double timePutHeader = 0.0;
char baseDirectory[512];
int setCounter = 0;
int sampleCounter = 0;
int eventCounter = 0;
char endianness[10];

double getCurrentTime() {
#if defined(PLATFORM_WINDOWS) && !defined(COMPILER_MINGW)
	return timeGetTime() * 0.001;
#else
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + tv.tv_usec*1e-6;
#endif
}

int my_request_handler(const message_t *request, message_t **response, void *user_data) {

	double tAbs, tRel;
	int res;
	int quantity;

	tAbs = getCurrentTime();
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
				quantity = ddef->nsamples;
				printf("Put data, %i channels, %i samples, type=%i ... ", ddef->nchans, ddef->nsamples, ddef->data_type);
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
			quantity = check_event_array(request->def->bufsize, request->buf);
			printf("Put events, number = %i, bufsize = %i ... ", quantity, request->def->bufsize);
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
				tAbs = getCurrentTime();
				tRel = tAbs - timePutHeader;
				printf("t=%8.3f WAIT_OK\n", tRel);
				break;
			case WAIT_ERR:
				printf("WAIT_ERR\n");
				break;
			case PUT_OK:
				printf("OK\n");
				pthread_mutex_lock(&qMutex);
				if (queue[qWritePos].command == 0) {
					queue[qWritePos].command = request->def->command;
					queue[qWritePos].quantity = quantity;
					queue[qWritePos].t = tRel;
					if (++qWritePos == QUEUE_SIZE) qWritePos = 0;
				} else {
					printf("WARNING - saving thread does not keep up!\n");
				}
				pthread_mutex_unlock(&qMutex);
				break;
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


int write_contents() {
	char name[512];
	FILE *f;
	int r;

	#if defined(PLATFORM_WINDOWS) && !defined(COMPILER_MINGW)
	r = GetFileAttributes(baseDirectory);
	#else
	r = access(baseDirectory, F_OK);
	#endif
	if (r!=-1) {
		fprintf(stderr, "ERROR: %s already exists in path\n", baseDirectory);
		return 0;
	}

	#if defined(WIN32)
	r = mkdir(baseDirectory);
	#else
	r = mkdir(baseDirectory, 0700);
	#endif
	if (r==-1) {
		fprintf(stderr, "ERROR: cannot create directory %s\n", name);
		return 0;
	}

	snprintf(name, sizeof(name), "%s/contents.txt", baseDirectory);
	f = fopen(name, "w");
	if (f==NULL) {
		fprintf(stderr, "Contents file can not be written - please check base directory.\n");
		return 0;
	}
	fputs(content_descr, f);
	fclose(f);
	return 1;
}

int write_header_to_disk() {
	char name[512];
	messagedef_t reqdef = {VERSION, GET_HDR, 0};
	message_t request;
	message_t *response;
	headerdef_t *hdef;
	int r;

	request.def = &reqdef;
	request.buf = NULL;

	r = dmarequest(&request, &response);
	if (r!=0 || response == NULL || response->def == NULL || response->buf == NULL) {
		fprintf(stderr, "ERROR: Cannot retrieve header for writing to disk\n");
		goto cleanup;
	}

	hdef = (headerdef_t *) response->buf;

	setCounter++;
	sampleCounter = eventCounter = 0;
	snprintf(name, sizeof(name), "%s/%04i", baseDirectory, setCounter);

	OS = ft_storage_create(name, hdef, hdef+1, &r);

cleanup:
	if (response!=NULL) {
		if (response->def != NULL) free(response->def);
		if (response->buf != NULL) free(response->buf);
		free(response);
	}
	return r;
}


int write_samples_to_disk(int nsamps, double t) {
	messagedef_t reqdef;
	datasel_t ds;
	message_t request;
	message_t *response;
	datadef_t *ddef;
	int r;

	ds.begsample = sampleCounter;
	ds.endsample = sampleCounter + nsamps - 1;
	reqdef.version = VERSION;
	reqdef.command = GET_DAT;
	reqdef.bufsize = sizeof(ds);
	request.def = &reqdef;
	request.buf = &ds;

	r = dmarequest(&request, &response);
	if (r!=0 || response == NULL || response->def == NULL || response->buf == NULL) {
		fprintf(stderr, "ERROR: Cannot retrieve samples for writing to disk\n");
		goto cleanup;
	}

	sampleCounter += nsamps;

	ddef = (datadef_t *) response->buf;
	r = ft_storage_add_samples(OS, nsamps, ddef+1);
	if (r==0) {
		ft_timing_element_t te;
		te.numSamples = nsamps;
		te.numEvents = 0;
		te.time = t;
		r = ft_storage_add_timing(OS, &te);
	}

cleanup:
	if (response!=NULL) {
		if (response->def != NULL) free(response->def);
		if (response->buf != NULL) free(response->buf);
		free(response);
	}
	return r;
}

int write_events_to_disk(int nevs, double t) {
	messagedef_t reqdef;
	eventsel_t es;
	message_t request;
	message_t *response;
	int r;

	es.begevent = eventCounter;
	es.endevent = eventCounter + nevs - 1;
	reqdef.version = VERSION;
	reqdef.command = GET_EVT;
	reqdef.bufsize = sizeof(es);
	request.def = &reqdef;
	request.buf = &es;

	r = dmarequest(&request, &response);
	if (r!=0 || response == NULL || response->def == NULL || response->buf == NULL) {
		fprintf(stderr, "ERROR: Cannot retrieve events for writing to disk\n");
		goto cleanup;
	}

	eventCounter += nevs;

	r = ft_storage_add_events(OS, response->def->bufsize, response->buf);
	if (r==0) {
		ft_timing_element_t te;
		te.numSamples = 0;
		te.numEvents = nevs;
		te.time = t;
		r = ft_storage_add_timing(OS, &te);
	}
cleanup:
	if (response!=NULL) {
		if (response->def != NULL) free(response->def);
		if (response->buf != NULL) free(response->buf);
		free(response);
	}
	return r;
}



int main(int argc, char *argv[]) {
	int port;
	char *name = NULL;
	union {
		short word;
		char bytes[2];
	} endianTest;

	endianTest.bytes[0] = 1;
	endianTest.bytes[1] = 0;

	if (endianTest.word == 1) {
		strcpy(endianness, "little");
	} else {
		strcpy(endianness, "big");
	}

	/* verify that all datatypes have the expected syze in bytes */
	check_datatypes();

	#if defined(PLATFORM_WINDOWS) && !defined(COMPILER_MINGW)
	timeBeginPeriod(1);
	#endif

	if (argc<2) {
		fprintf(stderr, "Usage: recording <directory> [port/unix socket]\n");
		return 1;
	}

	strncpy(baseDirectory, argv[1], sizeof(baseDirectory));

	if (argc>2) {
		port = atoi(argv[2]);
		if (port == 0) {
			name = argv[2];
		}
	} else {
		port = 1972;
	}

	memset(queue, sizeof(queue), 0);

	if (!write_contents()) goto cleanup;

	S = ft_start_buffer_server(port, name, my_request_handler, NULL);
	if (S==NULL) return 1;

	signal(SIGINT, abortHandler);
	while (keepRunning) {
		if (queue[qReadPos].command == 0) {
			usleep(1000);
		} else {
			switch(queue[qReadPos].command) {
				case PUT_HDR:
					if (setCounter>0) ft_storage_close(OS);
					write_header_to_disk();
					if (OS==NULL) {
						fprintf(stderr, "!!!!!!!!!!!!! WARNING !!!!!!!!!!\n");
						fprintf(stderr, "!! will NOT save this dataset !!\n");
						fprintf(stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
					}
					break;
				case PUT_DAT:
					if (OS) write_samples_to_disk(queue[qReadPos].quantity, queue[qReadPos].t);
					break;
				case PUT_EVT:
					if (OS) write_events_to_disk(queue[qReadPos].quantity, queue[qReadPos].t);
					break;
			}
			queue[qReadPos].command = 0;
			if (++qReadPos == QUEUE_SIZE) qReadPos = 0;
		}
	}
	printf("Ctrl-C pressed -- stopping buffer server...\n");
	if (setCounter>0 && OS!=NULL) ft_storage_close(OS);
	ft_stop_buffer_server(S);
	printf("Done.\n");

cleanup:
	#if defined(PLATFORM_WINDOWS) && !defined(COMPILER_MINGW)
	timeEndPeriod(1);
	#endif

	return 0;
}
