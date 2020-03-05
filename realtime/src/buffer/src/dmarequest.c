/*
 * Copyright (C) 2008, Robert Oostenveld
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <time.h>  /* for clock */
#include <sys/time.h>

#include "buffer.h"
#include "platform_includes.h"

/* FIXME should these be static? */
static header_t   *header   = NULL;
static data_t     *data     = NULL;
static event_t    *event    = NULL;

/* these are used for fine-tuning the sample number of incoming events */
struct timespec putdat_clock;
struct timespec putevt_clock;

static unsigned int current_max_num_sample = 0;

static int thissample = 0;    /* points at the buffer */
static int thisevent = 0;     /* points at the buffer */

/* Note that there have been problems with the order of the mutexes (e.g.
 * http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=933).
 * I have attempted to make the order of locking consistent, but can't give
 * guarantees. A more long term solution could be:
 * - find the dependencies between modifications of volatile data (e.g. events
 * depend on header),
 * - keep locks as shortly as possible (get info, release again).
 *
 * This could results in a global lock (robust, probably not optimal in terms
 * of speed), or a series of locks sandwiching modification code in dedicated
 * functions.
 *
 * -- Boris
 */

pthread_mutex_t mutexheader   = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutexdata     = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutexevent    = PTHREAD_MUTEX_INITIALIZER;

pthread_cond_t getData_cond   = PTHREAD_COND_INITIALIZER;
pthread_mutex_t getData_mutex = PTHREAD_MUTEX_INITIALIZER;

/*****************************************************************************/

void free_header() {
	int verbose = 0;
	if (verbose>0) fprintf(stderr, "free_header: freeing header buffer\n");
	if (header) {
		FREE(header->def);
		FREE(header->buf);
		FREE(header);
	}
}

void free_data() {
	int verbose = 0;
	if (verbose>0) fprintf(stderr, "free_data: freeing data buffer\n");
	if (data) {
		FREE(data->def);
		FREE(data->buf);
		FREE(data);
	}
	thissample = 0;
	if (header) header->def->nsamples = 0;
}

void free_event() {
	int verbose = 0;
	int i;
	if (verbose>0) fprintf(stderr, "free_event: freeing event buffer\n");
	if (event) {
		for (i=0; i<MAXNUMEVENT; i++) {
			FREE(event[i].def);
			FREE(event[i].buf);
		}
		FREE(event);
	}
	thisevent = 0;
	if (header) header->def->nevents = 0;
}

/*****************************************************************************/

void init_data(void) {
	int verbose = 0;
	if (verbose>0) fprintf(stderr, "init_data: creating data buffer\n");
	if (header) {
		unsigned int wordsize = wordsize_from_type(header->def->data_type);

		if (wordsize==0) {
			fprintf(stderr, "init_data: unsupported data type (%u)\n", header->def->data_type);
			return;
		}
		/* heuristic of choosing size of buffer:
			 set current_max_num_sample to MAXNUMSAMPLE if nchans <= 256
			 otherwise, allocate about MAXNUMBYTE and calculate current_max_num_sample from nchans + wordsize
		 */
		if (header->def->nchans <= 256) {
			current_max_num_sample = MAXNUMSAMPLE;
		} else {
			current_max_num_sample = MAXNUMBYTE / (wordsize * header->def->nchans);
		}
		data = (data_t*)malloc(sizeof(data_t));

		DIE_BAD_MALLOC(data);

		data->def = (datadef_t*)malloc(sizeof(datadef_t));

		DIE_BAD_MALLOC(data->def);

		data->def->nchans    = header->def->nchans;
		data->def->nsamples  = current_max_num_sample;
		data->def->data_type = header->def->data_type;
		data->buf = malloc(header->def->nchans*current_max_num_sample*wordsize);

		DIE_BAD_MALLOC(data->buf);
	}
}

void init_event(void) {
	int verbose = 0;
	int i;
	if (verbose>0) fprintf(stderr, "init_event: creating event buffer\n");
	if (header) {
		event = (event_t*)malloc(MAXNUMEVENT*sizeof(event_t));
		DIE_BAD_MALLOC(event);
		for (i=0; i<MAXNUMEVENT; i++) {
			event[i].def = NULL;
			event[i].buf = NULL;
		}
	}
}


/*****************************************************************************
 * this function handles the direct memory access to the buffer
 * and copies objects to and from memory
 *****************************************************************************/
int dmarequest(const message_t *request, message_t **response_ptr) {
	unsigned int offset;
	/*
		 int blockrequest = 0;
	 */
	int verbose = 0;

	/* these are used for blocking the read requests */
	struct timeval tp;
	struct timespec ts;

	/* use a local variable for datasel (in GET_DAT) */
	datasel_t datasel;

	/* these are for typecasting */
	headerdef_t    *headerdef;
	datadef_t      *datadef;
	eventdef_t     *eventdef;
	eventsel_t     *eventsel;

	/* this will hold the response */
	message_t *response;
	response      = (message_t*)malloc(sizeof(message_t));

	/* check for "out of memory" problems */
	if (response == NULL) {
		*response_ptr = NULL;
		return -1;
	}
	response->def = (messagedef_t*)malloc(sizeof(messagedef_t));

	/* check for "out of memory" problems */
	if (response->def == NULL) {
		*response_ptr = NULL;
		free(response);
		return -1;
	}
	response->buf = NULL;
	/* the response should be passed to the calling function, where it should be freed */
	*response_ptr = response;

	if (verbose>1) print_request(request->def);

	switch (request->def->command) {

		case PUT_HDR:
		case PUT_HDR_NORESPONSE:
			if (verbose>1) fprintf(stderr, "dmarequest: PUT_HDR\n");
			pthread_mutex_lock(&mutexheader);
			pthread_mutex_lock(&mutexdata);
			pthread_mutex_lock(&mutexevent);

			headerdef = (headerdef_t*)request->buf;
			if (verbose>1) print_headerdef(headerdef);

			/* delete the old header, data and events */
			free_header();
			free_data();
			free_event();

			/* store the header and re-initialize */
			header      = (header_t*)malloc(sizeof(header_t));
			DIE_BAD_MALLOC(header);
			header->def = (headerdef_t*)malloc(sizeof(headerdef_t));
			DIE_BAD_MALLOC(header->def);
			header->buf = malloc(headerdef->bufsize);
			DIE_BAD_MALLOC(header->buf);
			memcpy(header->def, request->buf, sizeof(headerdef_t));
			memcpy(header->buf, (char*)request->buf+sizeof(headerdef_t), headerdef->bufsize);
			header->def->nsamples = 0;
			header->def->nevents  = 0;

			init_data();
			init_event();

			response->def->version = VERSION;
			response->def->bufsize = 0;
			/* check whether memory could indeed be allocated */
			if (data!= NULL && data->buf != NULL && data->def != NULL) {
				response->def->command = PUT_OK;
			} else {
				/* let's at least tell the client that something's wrong */
				response->def->command = PUT_ERR;
			}

			pthread_mutex_unlock(&mutexevent);
			pthread_mutex_unlock(&mutexdata);
			pthread_mutex_unlock(&mutexheader);
			break;

		case PUT_DAT:
		case PUT_DAT_NORESPONSE:
			if (verbose>1) fprintf(stderr, "dmarequest: PUT_DAT\n");
			pthread_mutex_lock(&mutexheader);
			pthread_mutex_lock(&mutexdata);

			datadef = (datadef_t*)request->buf;
			if (verbose>1) print_datadef(datadef);
			if (verbose>2) print_buf(request->buf, request->def->bufsize);

			response->def->version = VERSION;
			response->def->bufsize = 0;
			if (request->def->bufsize < sizeof(datadef_t))
				response->def->command = PUT_ERR;
			else if (header==NULL || data==NULL)
				response->def->command = PUT_ERR;
			else if (header->def->nchans != datadef->nchans)
				response->def->command = PUT_ERR;
			else if (header->def->data_type != datadef->data_type)
				response->def->command = PUT_ERR;
			else if (datadef->nsamples > current_max_num_sample)
				response->def->command = PUT_ERR;
			else {
				unsigned int i;
				unsigned int wordsize = wordsize_from_type(header->def->data_type);
				unsigned int datasize = wordsize * datadef->nsamples * datadef->nchans;

				response->def->command = PUT_OK;

				if (wordsize == 0) {
					fprintf(stderr, "dmarequest: unsupported data type (%d)\n", datadef->data_type);
					response->def->command = PUT_ERR;
				} else if (datasize > datadef->bufsize || (datadef->bufsize + sizeof(datadef_t)) > request->def->bufsize) {
					fprintf(stderr, "dmarequest: invalid size definitions in PUT_DAT request\n");
					response->def->command = PUT_ERR;
				} else {

					/* record the time at which the data was received */
					if (clock_gettime(CLOCK_REALTIME, &putdat_clock) != 0) {
						perror("clock_gettime");
						return -1;
					}

					/* number of bytes per sample (all channels) is given by wordsize x number of channels */
					unsigned int chansize = wordsize * data->def->nchans;
					/* request_data points to actual data samples within the request, use char* for convenience */
					const char *request_data = (const char *) request->buf + sizeof(datadef_t);
					char *buffer_data = (char *)data->buf;

					for (i=0; i<datadef->nsamples; i++) {
						memcpy(buffer_data+(thissample*chansize), request_data+(i*chansize), chansize);
						header->def->nsamples++;
						thissample++;
						thissample = WRAP(thissample, current_max_num_sample);
					}
					/* Signal possibly waiting threads that we have received data */
					pthread_cond_broadcast(&getData_cond);
				}
			}

			pthread_mutex_unlock(&mutexdata);
			pthread_mutex_unlock(&mutexheader);
			break;

		case PUT_EVT:
		case PUT_EVT_NORESPONSE:
			if (verbose>1) fprintf(stderr, "dmarequest: PUT_EVT\n");
			pthread_mutex_lock(&mutexheader);
			pthread_mutex_lock(&mutexevent);

			/* record the time at which the event was received */
			if (clock_gettime(CLOCK_REALTIME, &putevt_clock) != 0) {
				perror("clock_gettime");
				return -1;
			}

			/* Give an error message if there is no header, or if the given event array is defined badly */
			if (header==NULL || event==NULL || check_event_array(request->def->bufsize, request->buf) < 0) {
				response->def->version = VERSION;
				response->def->command = PUT_ERR;
				response->def->bufsize = 0;
			}
			else {	/* go over all events and store them one by one */
				response->def->version = VERSION;
				response->def->command = PUT_OK;
				response->def->bufsize = 0;

				offset = 0; /* this represents the offset of the event in the buffer */
				while (offset<request->def->bufsize) {
					FREE(event[thisevent].def);
					FREE(event[thisevent].buf);

					eventdef = (eventdef_t*)((char*)request->buf+offset);
					if (verbose>1) print_eventdef(eventdef);

					event[thisevent].def = (eventdef_t*)malloc(sizeof(eventdef_t));
					DIE_BAD_MALLOC(event[thisevent].def);
					memcpy(event[thisevent].def, (char*)request->buf+offset, sizeof(eventdef_t));

					if (event[thisevent].def->sample == EVENT_AUTO_SAMPLE) {
						/* automatically convert event->def->sample to current sample number */
						/* make some fine adjustment of the assigned sample number */
						double adjust = (putevt_clock.tv_sec - putdat_clock.tv_sec) + (double)(putevt_clock.tv_nsec - putdat_clock.tv_nsec) / 1000000000L;
						event[thisevent].def->sample = header->def->nsamples + (int)(header->def->fsample*adjust);
					}

					offset += sizeof(eventdef_t);
					event[thisevent].buf = malloc(eventdef->bufsize);
					DIE_BAD_MALLOC(event[thisevent].buf);
					memcpy(event[thisevent].buf, (char*)request->buf+offset, eventdef->bufsize);
					offset += eventdef->bufsize;
					if (verbose>1) print_eventdef(event[thisevent].def);
					thisevent++;
					thisevent = WRAP(thisevent, MAXNUMEVENT);
					header->def->nevents++;
				}
			}

			pthread_mutex_unlock(&mutexevent);
			pthread_mutex_unlock(&mutexheader);
			break;

		case GET_HDR:
			if (verbose>1) fprintf(stderr, "dmarequest: GET_HDR\n");
			if (header==NULL) {
				response->def->version = VERSION;
				response->def->command = GET_ERR;
				response->def->bufsize = 0;
				break;
			}

			pthread_mutex_lock(&mutexheader);

			response->def->version = VERSION;
			response->def->command = GET_OK;
			response->def->bufsize = 0;
			response->def->bufsize = append(&response->buf, response->def->bufsize, header->def, sizeof(headerdef_t));
			response->def->bufsize = append(&response->buf, response->def->bufsize, header->buf, header->def->bufsize);

			pthread_mutex_unlock(&mutexheader);
			break;

		case GET_DAT:
			if (verbose>1) fprintf(stderr, "dmarequest: GET_DAT\n");
			if (header==NULL || data==NULL) {
				response->def->version = VERSION;
				response->def->command = GET_ERR;
				response->def->bufsize = 0;
				break;
			}

			pthread_mutex_lock(&mutexheader);
			pthread_mutex_lock(&mutexdata);

			if (request->def->bufsize) {
				/* the selection has been specified */
				memcpy(&datasel, request->buf, sizeof(datasel_t));
				/* If endsample is -1 read the buffer to the end */
				if(datasel.endsample == -1)
				{
					datasel.endsample = header->def->nsamples - 1;
				}
			}
			else {
				/* determine a valid selection */
				if (header->def->nsamples>current_max_num_sample) {
					/* the ringbuffer is completely full */
					datasel.begsample = header->def->nsamples - current_max_num_sample;
					datasel.endsample = header->def->nsamples - 1;
				}
				else {
					/* the ringbuffer is not yet completely full */
					datasel.begsample = 0;
					datasel.endsample = header->def->nsamples - 1;
				}
			}

			/*

			// if the read should block...
			if(blockrequest == 1)
			{
			// check whether data is available
			while((datasel.begsample >= (datasel.endsample+1)) || (datasel.endsample > header->def->nsamples - 1))
			{
			// if not unlock all mutexes
			pthread_mutex_unlock(&mutexdata);
			pthread_mutex_unlock(&mutexheader);

			// wait for the condition to be signaled
			pthread_mutex_lock(&getData_mutex);
			gettimeofday(&tp, NULL);
			ts.tv_sec = tp.tv_sec;
			ts.tv_nsec = tp.tv_usec * 1000;
			ts.tv_sec += 1;
			pthread_cond_timedwait(&getData_cond, &getData_mutex, &ts);
			pthread_mutex_unlock(&getData_mutex);

			// Lock the mutexes again
			pthread_mutex_lock(&mutexheader);
			pthread_mutex_lock(&mutexdata);
			if(datasel.begsample == (datasel.endsample+1))
			datasel.endsample = header->def->nsamples - 1;
			}
			}
			 */

			if (verbose>1) print_headerdef(header->def);
			if (verbose>1) print_datasel(&datasel);

			if (datasel.begsample < 0 || datasel.endsample < 0) {
				fprintf(stderr, "dmarequest: err1\n");
				response->def->version = VERSION;
				response->def->command = GET_ERR;
				response->def->bufsize = 0;
			}
			else if (datasel.begsample >= header->def->nsamples || datasel.endsample >= header->def->nsamples) {
				fprintf(stderr, "dmarequest: err2\n");
				response->def->version = VERSION;
				response->def->command = GET_ERR;
				response->def->bufsize = 0;
			}
			else if ((header->def->nsamples - datasel.begsample) > current_max_num_sample) {
				fprintf(stderr, "dmarequest: err3\n");
				response->def->version = VERSION;
				response->def->command = GET_ERR;
				response->def->bufsize = 0;
			}
			else {
				unsigned int wordsize = wordsize_from_type(data->def->data_type);
				if (wordsize==0) {
					fprintf(stderr, "dmarequest: unsupported data type (%d)\n", data->def->data_type);
					response->def->version = VERSION;
					response->def->command = GET_ERR;
					response->def->bufsize = 0;
				}  else {
					unsigned int n;
					response->def->version = VERSION;
					response->def->command = GET_OK;
					response->def->bufsize = 0;

					/* determine the number of samples to return */
					n = datasel.endsample - datasel.begsample + 1;

					response->buf = malloc(sizeof(datadef_t) + n*data->def->nchans*wordsize);
					if (response->buf == NULL) {
						/* not enough space for copying data into response */
						fprintf(stderr, "dmarequest: out of memory\n");
						response->def->command = GET_ERR;
					}
					else {
						/* number of bytes per sample (all channels) */
						unsigned int chansize = data->def->nchans * wordsize;

						/* convenience pointer to start of actual data in response */
						char *resp_data = ((char *) response->buf) + sizeof(datadef_t);

						/* this is the location of begsample within the ringbuffer */
						unsigned int start_index = 	WRAP(datasel.begsample, current_max_num_sample);

						/* have datadef point into the freshly allocated response buffer and directly
							 fill in the information */
						datadef = (datadef_t *) response->buf;
						datadef->nchans    = data->def->nchans;
						datadef->data_type = data->def->data_type;
						datadef->nsamples  = n;
						datadef->bufsize   = n*chansize;

						response->def->bufsize = sizeof(datadef_t) + datadef->bufsize;

						if (start_index + n <= current_max_num_sample) {
							/* we can copy everything in one go */
							memcpy(resp_data, (char*)(data->buf) + start_index*chansize, n*chansize);
						} else {
							/* need to wrap around at current_max_num_sample */
							unsigned int na = current_max_num_sample - start_index;
							unsigned int nb = n - na;

							memcpy(resp_data, (char*)(data->buf) + start_index*chansize, na*chansize);
							memcpy(resp_data + na*chansize, (char*)(data->buf), nb*chansize);

							/* printf("Wrapped around!\n"); */
						}
					}
				}
			}

			pthread_mutex_unlock(&mutexdata);
			pthread_mutex_unlock(&mutexheader);
			break;

		case GET_EVT:
			if (verbose>1) fprintf(stderr, "dmarequest: GET_EVT\n");
			if (header==NULL || event==NULL || header->def->nevents==0) {
				response->def->version = VERSION;
				response->def->command = GET_ERR;
				response->def->bufsize = 0;
				break;
			}

			pthread_mutex_lock(&mutexheader);
			pthread_mutex_lock(&mutexevent);

			eventsel = (eventsel_t*)malloc(sizeof(eventsel_t));
			DIE_BAD_MALLOC(eventsel);

			/* determine the selection */
			if (request->def->bufsize) {
				/* the selection has been specified */
				memcpy(eventsel, request->buf, sizeof(eventsel_t));
			}
			else {
				/* determine a valid selection */
				if (header->def->nevents>MAXNUMEVENT) {
					/* the ringbuffer is completely full */
					eventsel->begevent = header->def->nevents - MAXNUMEVENT;
					eventsel->endevent = header->def->nevents - 1;
				}
				else {
					/* the ringbuffer is not yet completely full */
					eventsel->begevent = 0;
					eventsel->endevent = header->def->nevents - 1;
				}
			}

			if (verbose>1) print_headerdef(header->def);
			if (verbose>1) print_eventsel(eventsel);

			if (eventsel==NULL) {
				response->def->version = VERSION;
				response->def->command = GET_ERR;
				response->def->bufsize = 0;
			}
			else if (eventsel->begevent < 0 || eventsel->endevent < 0) {
				fprintf(stderr, "dmarequest: err4\n");
				response->def->version = VERSION;
				response->def->command = GET_ERR;
				response->def->bufsize = 0;
			}
			else if (eventsel->begevent >= header->def->nevents || eventsel->endevent >= header->def->nevents) {
				fprintf(stderr, "dmarequest: err5\n");
				response->def->version = VERSION;
				response->def->command = GET_ERR;
				response->def->bufsize = 0;
			}
			else if ((header->def->nevents-eventsel->begevent) > MAXNUMEVENT) {
				fprintf(stderr, "dmarequest: err6\n");
				response->def->version = VERSION;
				response->def->command = GET_ERR;
				response->def->bufsize = 0;
			}
			else {
				unsigned int j,n;

				response->def->version = VERSION;
				response->def->command = GET_OK;
				response->def->bufsize = 0;

				/* determine the number of events to return */
				n = eventsel->endevent - eventsel->begevent + 1;

				for (j=0; j<n; j++) {
					if (verbose>1) print_eventdef(event[WRAP(eventsel->begevent+j, MAXNUMEVENT)].def);
					response->def->bufsize = append(&response->buf, response->def->bufsize, event[WRAP(eventsel->begevent+j, MAXNUMEVENT)].def, sizeof(eventdef_t));
					response->def->bufsize = append(&response->buf, response->def->bufsize, event[WRAP(eventsel->begevent+j, MAXNUMEVENT)].buf, event[WRAP(eventsel->begevent+j, MAXNUMEVENT)].def->bufsize);
				}
			}

			FREE(eventsel);
			pthread_mutex_unlock(&mutexevent);
			pthread_mutex_unlock(&mutexheader);
			break;

		case FLUSH_HDR:
			pthread_mutex_lock(&mutexheader);
			pthread_mutex_lock(&mutexdata);
			pthread_mutex_lock(&mutexevent);
			if (header) {
				free_header();
				free_data();
				free_event();
				response->def->version = VERSION;
				response->def->command = FLUSH_OK;
				response->def->bufsize = 0;
			}
			else {
				response->def->version = VERSION;
				response->def->command = FLUSH_ERR;
				response->def->bufsize = 0;
			}
			pthread_mutex_unlock(&mutexevent);
			pthread_mutex_unlock(&mutexdata);
			pthread_mutex_unlock(&mutexheader);
			break;

		case FLUSH_DAT:
			pthread_mutex_lock(&mutexheader);
			pthread_mutex_lock(&mutexdata);
			if (header && data) {
				header->def->nsamples = thissample = 0;
				response->def->version = VERSION;
				response->def->command = FLUSH_OK;
				response->def->bufsize = 0;
			}
			else {
				response->def->version = VERSION;
				response->def->command = FLUSH_ERR;
				response->def->bufsize = 0;
			}
			pthread_mutex_unlock(&mutexdata);
			pthread_mutex_unlock(&mutexheader);
			break;

		case FLUSH_EVT:
			pthread_mutex_lock(&mutexheader);
			pthread_mutex_lock(&mutexevent);
			if (header && event) {
				unsigned int i;

				header->def->nevents = thisevent = 0;
				for (i=0; i<MAXNUMEVENT; i++) {
					FREE(event[i].def);
					FREE(event[i].buf);
				}
				response->def->version = VERSION;
				response->def->command = FLUSH_OK;
				response->def->bufsize = 0;
			}
			else {
				response->def->version = VERSION;
				response->def->command = FLUSH_ERR;
				response->def->bufsize = 0;
			}
			pthread_mutex_unlock(&mutexevent);
			pthread_mutex_unlock(&mutexheader);
			break;

		case WAIT_DAT:
			/* SK: This request means that the client wants to wait until
				 MORE than waitdef_t.threshold.nsamples samples OR
				 MORE THAN waitdef_t.threshold.nevents events
				 are in the buffer, BUT
				 only for the time given in waitdef_t.milliseconds.
				 The response is just the number of samples and events
				 in the buffer as described by samples_events_t.
			 */
			response->def->version = VERSION;
			if (header==NULL || request->def->bufsize!=sizeof(waitdef_t)) {
				response->def->command = WAIT_ERR;
				response->def->bufsize = 0;
			} else {
				int waiterr;
				waitdef_t *wd = (waitdef_t *) request->buf;
				samples_events_t *nret = malloc(sizeof(samples_events_t));
				UINT32_T nsmp, nevt;

				if (nret == NULL) {
					/* highly unlikely, but we cannot allocate a sample_event_t - return an error */
					response->def->command = WAIT_ERR;
					response->def->bufsize = 0;
					break;
				}
				/* Let response->buf point to the new sample_event_t structure */
				response->def->command = WAIT_OK;
				response->def->bufsize = sizeof(samples_events_t);
				response->buf = nret;

				/* get current number of samples */
				pthread_mutex_lock(&mutexheader);
				nsmp = header->def->nsamples;
				nevt = header->def->nevents;
				pthread_mutex_unlock(&mutexheader);

				if (wd->milliseconds == 0 || nsmp > wd->threshold.nsamples || nevt > wd->threshold.nevents) {
					/* the client doesn't want to wait, or
						 we're already above the threshold:
						 return immediately */
					nret->nsamples = nsmp;
					nret->nevents = nevt;
					break;
				}
				gettimeofday(&tp, NULL);
				ts.tv_sec = tp.tv_sec + (wd->milliseconds/1000);
				ts.tv_nsec = 1000 * (tp.tv_usec + (wd->milliseconds % 1000)*1000);
				while (ts.tv_nsec >= 1000000000) {
					ts.tv_sec++;
					ts.tv_nsec-=1000000000;
				}

				/* FIXME: The getData condition variable is only triggered by incoming data, not events */
				do {
					pthread_mutex_lock(&getData_mutex);
					waiterr = pthread_cond_timedwait(&getData_cond, &getData_mutex, &ts);
					pthread_mutex_unlock(&getData_mutex);

					/* get current number of samples */
					pthread_mutex_lock(&mutexheader);
					nsmp = header->def->nsamples;
					nevt = header->def->nevents;
					pthread_mutex_unlock(&mutexheader);
				} while (nsmp <= wd->threshold.nsamples && nevt <= wd->threshold.nevents && waiterr==0);
				nret->nsamples = nsmp;
				nret->nevents = nevt;
			}
			break;

		default:
			fprintf(stderr, "dmarequest: unknown command\n");
	}

	if (verbose>0) fprintf(stderr, "dmarequest: thissample = %u, thisevent = %u\n", thissample, thisevent);

	/* everything went fine */
	return 0;
}
