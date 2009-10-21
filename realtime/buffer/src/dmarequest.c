/*
 * Copyright (C) 2008, Robert Oostenveld
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 * $Log: dmarequest.c,v $
 * Revision 1.20  2009/06/17 13:41:36  roboos
 * only print if verbose==1
 *
 * Revision 1.19  2009/05/13 20:45:11  roboos
 * switched off verbose mode
 *
 * Revision 1.18  2009/01/23 19:47:50  roboos
 * changed verbosity
 *
 * Revision 1.17  2009/01/23 08:26:44  roboos
 * fixed a serious bug that caused a lot of memory to leak (in fact all packets that were sent over the socket would eventually leack away), both on the client and server side
 *
 * Revision 1.16  2009/01/21 20:54:35  roboos
 * added some debugging code to keep track of number of sockets and threads (both seem ok)
 * cleaned up the debugging: more consistent use of verbose flag and function name in fprintf
 *
 * Revision 1.15  2008/07/09 11:16:17  roboos
 * removed winsock2 header, moved declaration of timespec and timeval to top
 *
 * Revision 1.14  2008/07/01 17:09:00  thohar
 * added cleanuo_buffer function
 * made pthread_cond_wait a timed_wait to avoid deadlocks
 *
 * Revision 1.13  2008/06/20 07:51:30  roboos
 * renamed property GET_DAT_Block into dmaBlockRequest for consistency with other properties
 *
 * Revision 1.12  2008/06/19 20:46:40  roboos
 * fixed bug in GET_EVT due to unsigned ints in combination with inappropriate wrapping
 * made GET_DAT consistent with the improved GET_EVT implementatoin (thanks to better wrapping)
 *
 * Revision 1.11  2008/06/02 15:31:21  roboos
 * fixed typo in fprintf feedback
 *
 * Revision 1.10  2008/05/22 09:55:15  roboos
 * some changes for compatibility wioth Borland, thanks to Jurgen
 *
 * Revision 1.9  2008/04/24 15:53:15  roboos
 * changed verbosity
 *
 * Revision 1.8  2008/04/15 14:08:07  thohar
 * added possibility do do a blocking GET_DAT by setting property "GET_DAT_Block" to 1
 * if datasel->ensample == -1 all data until the end of the buffer is read out
 *
 * Revision 1.7  2008/03/23 13:18:17  roboos
 * removed and disabled some fprintf debug info
 *
 * Revision 1.6  2008/03/19 09:27:46  thohar
 * added possibility of debug output on GET_DAT and PUT_DAT for DATATYPE_FLOAT32
 *
 * Revision 1.5  2008/03/17 13:47:12  roboos
 * fixed various problems with get/put property:
 * - unlock mutex
 * - init propertybuf at other place in code
 * - find property memcpy for various types
 * - allow get without propertysel, returns all properties
 *
 * Revision 1.4  2008/03/13 13:39:07  roboos
 * made the global variables static, not sure though
 *
 * Revision 1.3  2008/03/13 12:59:46  thohar
 * initialised thissample, thisevent and thisproperty to 0. Needed to build on MAC
 *
 * Revision 1.2  2008/03/10 09:39:32  roboos
 * added some property handling
 *
 * Revision 1.1  2008/03/08 10:31:50  roboos
 * renamed handle_request to dmarequest, this function can now be called by tcpsocket() in the buffer-server, and/or by clientrequest() in the client
 *
 * Revision 1.2  2008/03/07 14:51:54  roboos
 * added comment, no functional change
 *
 * Revision 1.1  2008/03/02 13:29:23  roboos
 * seperated socket and buffer-handing code
 * adde polling to the socket section
 * allow connections to remain open (i.e. multiple request/response pairs)
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "unix_includes.h"
#include <pthread.h>

#include "buffer.h"
#include "message.h"

// FIXME should these be static?
static header_t   *header   = NULL;
static data_t     *data     = NULL;
static event_t    *event    = NULL;
static property_t *property = NULL;

static int thissample = 0;    // points at the buffer
static int thisevent = 0;     // points at the buffer
static int thisproperty = 0;  // points at the buffer

pthread_mutex_t mutexheader   = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutexdata     = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutexevent    = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutexproperty = PTHREAD_MUTEX_INITIALIZER;

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

void free_property() {
	int verbose = 0;
	int i;
	if (verbose>0) fprintf(stderr, "free_property: freeing property buffer\n");
	if (property) {
		for (i=0; i<MAXNUMPROPERTY; i++) {
			FREE(property[i].def);
			FREE(property[i].buf);
		}
		FREE(property);
	}
	thisproperty = 0;
}

/*****************************************************************************/

void init_data(void) {
	int verbose = 0;
	if (verbose>0) fprintf(stderr, "init_data: creating data buffer\n");
	if (header) {
		data = (data_t*)malloc(sizeof(data_t));
		data->def = (datadef_t*)malloc(sizeof(datadef_t));
		data->def->nchans    = header->def->nchans;
		data->def->nsamples  = MAXNUMSAMPLE;
		data->def->data_type = header->def->data_type;
		switch (header->def->data_type) {
			case DATATYPE_INT8:
				data->buf = malloc(header->def->nchans*MAXNUMSAMPLE*WORDSIZE_INT8);
				break;

			case DATATYPE_INT16:
				data->buf = malloc(header->def->nchans*MAXNUMSAMPLE*WORDSIZE_INT16);
				break;

			case DATATYPE_INT32:
				data->buf = malloc(header->def->nchans*MAXNUMSAMPLE*WORDSIZE_INT32);
				break;

			case DATATYPE_INT64:
				data->buf = malloc(header->def->nchans*MAXNUMSAMPLE*WORDSIZE_INT64);
				break;

			case DATATYPE_FLOAT32:
				data->buf = malloc(header->def->nchans*MAXNUMSAMPLE*WORDSIZE_FLOAT32);
				break;

			case DATATYPE_FLOAT64:
				data->buf = malloc(header->def->nchans*MAXNUMSAMPLE*WORDSIZE_FLOAT64);
				break;

			default:
				fprintf(stderr, "init_data: unsupported data type (%u)\n", header->def->data_type);
				free_data();
		}
	}
}

void init_event(void) {
	int verbose = 0;
	int i;
	if (verbose>0) fprintf(stderr, "init_event: creating event buffer\n");
	if (header) {
		event = (event_t*)malloc(MAXNUMEVENT*sizeof(event_t));
		for (i=0; i<MAXNUMEVENT; i++) {
			event[i].def = NULL;
			event[i].buf = NULL;
		}
	}
}

void init_property(void) {
	int verbose = 0;
	int i;
	if (verbose>0) fprintf(stderr, "init_event: creating property buffer\n");
	property = (property_t*)malloc(MAXNUMPROPERTY*sizeof(property_t));
	for (i=0; i<MAXNUMPROPERTY; i++) {
		property[i].def = NULL;
		property[i].buf = NULL;
	}
}

/*****************************************************************************/

int find_property(property_t *desired) {
	int i, n = -1;
	if (property)
		for (i=0; i<MAXNUMPROPERTY; i++) {
			if (property[i].def==NULL || property[i].buf==NULL)
				continue;
			if ((property[i].def->type_type==desired->def->type_type) &&
				(property[i].def->type_numel==desired->def->type_numel)) {
				switch (desired->def->type_type) {
					case DATATYPE_CHAR:
						n = (memcmp(property[i].buf, desired->buf, desired->def->type_numel*WORDSIZE_CHAR)==0 ? i : -1);
						break;
					case DATATYPE_INT8:
					case DATATYPE_INT16:
					case DATATYPE_INT32:
					case DATATYPE_INT64:
					case DATATYPE_UINT8:
					case DATATYPE_UINT16:
					case DATATYPE_UINT32:
					case DATATYPE_UINT64:
					case DATATYPE_FLOAT32:
					case DATATYPE_FLOAT64:
					default:
						fprintf(stderr, "find_property: unsupported type\n");
				}
			}
			if (n>=0)
				// the desired property has been found
				break;
		}
	return n;
}

/***************************************************************************** 
 * this function handles the direct memory access to the buffer
 * and copies objects to and from memory
 *****************************************************************************/
int dmarequest(message_t *request, message_t **response_ptr) {
	int i, j, n, offset;
    int blockrequest = 0;
	int verbose = 0;

    // these are used for blocking the read requests
    struct timeval tp;
	struct timespec ts;

	// these are for typecasting
	headerdef_t    *headerdef;
	datadef_t      *datadef;
	datasel_t      *datasel;
	eventdef_t     *eventdef;
	eventsel_t     *eventsel;
	propertydef_t  *propertydef;

	// this is used to find a given property in the buffer
	property_t *desired;

	// this will hold the response
	message_t *response;
	response      = (message_t*)malloc(sizeof(message_t));
	response->def = (messagedef_t*)malloc(sizeof(messagedef_t));
	response->buf = NULL;
	// the response should be passed to the calling function, where it should be freed
	*response_ptr = response;

    if (verbose>1) print_request(request->def);

	switch (request->def->command) {

		case PUT_HDR:
			if (verbose>1) fprintf(stderr, "dmarequest: PUT_HDR\n");
			pthread_mutex_lock(&mutexheader);
			pthread_mutex_lock(&mutexdata);
			pthread_mutex_lock(&mutexevent);
			pthread_mutex_lock(&mutexproperty);

			headerdef = (headerdef_t*)request->buf;
			if (verbose>1) print_headerdef(headerdef);

			// delete the old header, data and events
			free_header();
			free_data();
			free_event();

			// store the header and re-initialize
			header      = (header_t*)malloc(sizeof(header_t));
			header->def = (headerdef_t*)malloc(sizeof(headerdef_t));
			header->buf = malloc(headerdef->bufsize);
			memcpy(header->def, request->buf, sizeof(headerdef_t));
			memcpy(header->buf, (char*)request->buf+sizeof(headerdef_t), headerdef->bufsize);
			header->def->nsamples = 0;
			header->def->nevents  = 0;

			init_data();
			init_event();

			response->def->version = VERSION;
			response->def->command = PUT_OK;
			response->def->bufsize = 0;

			pthread_mutex_unlock(&mutexheader);
			pthread_mutex_unlock(&mutexdata);
			pthread_mutex_unlock(&mutexevent);
			pthread_mutex_unlock(&mutexproperty);
			break;

		case PUT_DAT:
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
			else if (datadef->nsamples > MAXNUMSAMPLE)
				response->def->command = PUT_ERR;
			else {
				response->def->command = PUT_OK;
				for (i=0; i<datadef->nsamples; i++) {
					switch (datadef->data_type) {
						case DATATYPE_INT8:
							memcpy((char*)(data->buf)+(thissample*data->def->nchans)*sizeof(INT8_T), (char*)request->buf+sizeof(datadef_t)+(i*data->def->nchans)*sizeof(INT8_T), sizeof(INT8_T)*data->def->nchans);
							break;

						case DATATYPE_INT16:
							memcpy((char*)(data->buf)+(thissample*data->def->nchans)*sizeof(INT16_T), (char*)request->buf+sizeof(datadef_t)+(i*data->def->nchans)*sizeof(INT16_T), sizeof(INT16_T)*data->def->nchans);
							break;

						case DATATYPE_INT32:
							memcpy((char*)(data->buf)+(thissample*data->def->nchans)*sizeof(INT32_T), (char*)request->buf+sizeof(datadef_t)+(i*data->def->nchans)*sizeof(INT32_T), sizeof(INT32_T)*data->def->nchans);
							break;

						case DATATYPE_INT64:
							memcpy((char*)(data->buf)+(thissample*data->def->nchans)*sizeof(INT64_T), (char*)request->buf+sizeof(datadef_t)+(i*data->def->nchans)*sizeof(INT64_T), sizeof(INT64_T)*data->def->nchans);
							break;

						case DATATYPE_FLOAT32:
							memcpy((char*)(data->buf)+(thissample*data->def->nchans)*sizeof(FLOAT32_T), (char*)request->buf+sizeof(datadef_t)+(i*data->def->nchans)*sizeof(FLOAT32_T), sizeof(FLOAT32_T)*data->def->nchans);
							break;

						case DATATYPE_FLOAT64:
							memcpy((char*)(data->buf)+(thissample*data->def->nchans)*sizeof(FLOAT64_T), (char*)request->buf+sizeof(datadef_t)+(i*data->def->nchans)*sizeof(FLOAT64_T), sizeof(FLOAT64_T)*data->def->nchans);
							break;

						default:
							fprintf(stderr, "dmarequest: unsupported data type (%d)\n", datadef->data_type);
							response->def->command = PUT_ERR;
							continue; // thissample and header->def->nsamples will not be incremented
					}
					header->def->nsamples++;
					thissample++;
					thissample = WRAP(thissample, MAXNUMSAMPLE);
				}
				// Signal possibly waiting threads that we have received data
				pthread_cond_broadcast(&getData_cond);
			}

			pthread_mutex_unlock(&mutexdata);
			pthread_mutex_unlock(&mutexheader);
			break;

		case PUT_EVT:
			if (verbose>1) fprintf(stderr, "dmarequest: PUT_EVT\n");
			pthread_mutex_lock(&mutexheader);
			pthread_mutex_lock(&mutexevent);

			// go over all events and store them one by one
			if (header==NULL || event==NULL) {
				response->def->version = VERSION;
				response->def->command = PUT_ERR;
				response->def->bufsize = 0;
			}
			else {
				response->def->version = VERSION;
				response->def->command = PUT_OK;
				response->def->bufsize = 0;

				offset = 0; // this represents the offset of the event in the buffer
				while (offset<request->def->bufsize) {
					FREE(event[thisevent].def);
					FREE(event[thisevent].buf);

					eventdef = (eventdef_t*)((char*)request->buf+offset);
					if (verbose>1) print_eventdef(eventdef);

					event[thisevent].def = (eventdef_t*)malloc(sizeof(eventdef_t));
					memcpy(event[thisevent].def, (char*)request->buf+offset, sizeof(eventdef_t));
					offset += sizeof(eventdef_t);
					event[thisevent].buf = malloc(eventdef->bufsize);
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

		case PUT_PRP:
			if (verbose>1) fprintf(stderr, "dmarequest: PUT_PRP\n");
			if (request->def->bufsize==0) {
				response->def->version = VERSION;
				response->def->command = PUT_ERR;
				response->def->bufsize = 0;
				break;
			}

			pthread_mutex_lock(&mutexproperty);

			// the property buffer should exist seperate from header, data and events
			if (!property)
				init_property();

			propertydef = (propertydef_t*)request->buf;
			if (verbose>1) print_propertydef(propertydef);

			// determine whether the property already exists in the buffer
			// if it already exists, then it should be updated
			// otherwise it should be inserted as new property
			desired = (property_t*)malloc(sizeof(property_t));
			desired->def = (propertydef_t*)request->buf;
			if (desired->def->bufsize)
				desired->buf = (char*)request->buf+sizeof(propertydef_t);
			else
				desired->buf = NULL;
			n = find_property(desired);
			FREE(desired);

			if (n<0 && (thisproperty<MAXNUMPROPERTY)) {
				// insert as new property
				n = thisproperty;
				thisproperty++;
			}

			if (n>=0) {
				// clear the old property information (if any)
				FREE(property[n].def);
				FREE(property[n].buf);

				// insert the new property information
				property[n].def = (propertydef_t*)malloc(sizeof(propertydef_t));
				memcpy(property[n].def, request->buf, sizeof(propertydef_t));
				property[n].buf = malloc(property[n].def->bufsize);
				memcpy(property[n].buf, (char*)request->buf+sizeof(propertydef_t), property[n].def->bufsize);

				response->def->version = VERSION;
				response->def->command = PUT_OK;
				response->def->bufsize = 0;
			}
			else {
				response->def->version = VERSION;
				response->def->command = PUT_ERR;
				response->def->bufsize = 0;
			}

			pthread_mutex_unlock(&mutexproperty);
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
			
			// Check whether the read-request should block...
			blockrequest = 0;
			get_property(0, "dmaBlockRequest", &blockrequest);
			if (verbose>1) fprintf(stderr, "dmarequest: blockrequest = %d\n", blockrequest);

			pthread_mutex_lock(&mutexdata);
			pthread_mutex_lock(&mutexheader);

			datasel = (datasel_t*)malloc(sizeof(datasel_t));
			if (request->def->bufsize) {
				// the selection has been specified
				memcpy(datasel, request->buf, sizeof(datasel_t));
				// If endsample is -1 read the buffer to the end
				if(datasel->endsample == -1)
				{
					datasel->endsample = header->def->nsamples - 1;
				}
			}
			else {
				// determine a valid selection
				if (header->def->nsamples>MAXNUMSAMPLE) {
					// the ringbuffer is completely full
					datasel->begsample = header->def->nsamples - MAXNUMSAMPLE;
					datasel->endsample = header->def->nsamples - 1;
				}
				else {
					// the ringbuffer is not yet completely full
					datasel->begsample = 0;
					datasel->endsample = header->def->nsamples - 1;
				}
			}
			
			// if the read should block...
			if(blockrequest == 1)
			{
				// check whether data is available
				while((datasel->begsample >= (datasel->endsample+1)) || (datasel->endsample > header->def->nsamples - 1))
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
					pthread_mutex_lock(&mutexdata);
					pthread_mutex_lock(&mutexheader);
					if(datasel->begsample == (datasel->endsample+1))
						datasel->endsample = header->def->nsamples - 1;
				}
			}

			if (verbose>1) print_headerdef(header->def);
			if (verbose>1) print_datasel(datasel);

			if (datasel==NULL) {
				fprintf(stderr, "dmarequest: err0\n");
				response->def->version = VERSION;
				response->def->command = GET_ERR;
				response->def->bufsize = 0;
			}
			else if (datasel->begsample < 0 || datasel->endsample < 0) {
				fprintf(stderr, "dmarequest: err1\n");
				response->def->version = VERSION;
				response->def->command = GET_ERR;
				response->def->bufsize = 0;
			}
			else if (datasel->begsample >= header->def->nsamples || datasel->endsample >= header->def->nsamples) {
				fprintf(stderr, "dmarequest: err2\n");
				response->def->version = VERSION;
				response->def->command = GET_ERR;
				response->def->bufsize = 0;
			}
			else if ((header->def->nsamples - datasel->begsample) > MAXNUMSAMPLE) {
				fprintf(stderr, "dmarequest: err3\n");
				response->def->version = VERSION;
				response->def->command = GET_ERR;
				response->def->bufsize = 0;
			}
			else {
				// assume for the moment that it will be ok
				// the only problem might be an unsupported datatype, which is dealt with below
				response->def->version = VERSION;
				response->def->command = GET_OK;
				response->def->bufsize = 0;

				// determine the number of samples to return
				n = datasel->endsample - datasel->begsample + 1;

				switch (data->def->data_type) {
					case DATATYPE_INT8:
						data->def->nsamples = n;
						data->def->bufsize  = n*data->def->nchans*sizeof(INT8_T);
						response->def->bufsize = append(&response->buf, response->def->bufsize, data->def, sizeof(datadef_t));
						data->def->nsamples = MAXNUMSAMPLE;
						data->def->bufsize  = 0;
						for (j=0; j<n; j++) {
							response->def->bufsize = append(&response->buf, response->def->bufsize, (char*)(data->buf)+WRAP(datasel->begsample+j,MAXNUMSAMPLE)*(data->def->nchans)*sizeof(INT8_T), (data->def->nchans)*sizeof(INT8_T));
						}
						break;

					case DATATYPE_INT16:
						data->def->nsamples = n;
						data->def->bufsize  = n*data->def->nchans*sizeof(INT16_T);
						response->def->bufsize = append(&response->buf, response->def->bufsize, data->def, sizeof(datadef_t));
						data->def->nsamples = MAXNUMSAMPLE;
						data->def->bufsize  = 0;
						for (j=0; j<n; j++) {
							response->def->bufsize = append(&response->buf, response->def->bufsize, (char*)(data->buf)+WRAP(datasel->begsample+j,MAXNUMSAMPLE)*(data->def->nchans)*sizeof(INT16_T), (data->def->nchans)*sizeof(INT16_T));
						}
						break;

					case DATATYPE_INT32:
						data->def->nsamples = n;
						data->def->bufsize  = n*data->def->nchans*sizeof(INT32_T);
						response->def->bufsize = append(&response->buf, response->def->bufsize, data->def, sizeof(datadef_t));
						data->def->nsamples = MAXNUMSAMPLE;
						data->def->bufsize  = 0;
						for (j=0; j<n; j++) {
							response->def->bufsize = append(&response->buf, response->def->bufsize, (char*)(data->buf)+WRAP(datasel->begsample+j,MAXNUMSAMPLE)*(data->def->nchans)*sizeof(INT32_T), (data->def->nchans)*sizeof(INT32_T));
						}
						break;

					case DATATYPE_INT64:
						data->def->nsamples = n;
						data->def->bufsize  = n*data->def->nchans*sizeof(INT64_T);
						response->def->bufsize = append(&response->buf, response->def->bufsize, data->def, sizeof(datadef_t));
						data->def->nsamples = MAXNUMSAMPLE;
						data->def->bufsize  = 0;
						for (j=0; j<n; j++) {
							response->def->bufsize = append(&response->buf, response->def->bufsize, (char*)(data->buf)+WRAP(datasel->begsample+j,MAXNUMSAMPLE)*(data->def->nchans)*sizeof(INT64_T), (data->def->nchans)*sizeof(INT64_T));
						}
						break;

					case DATATYPE_FLOAT32:
						data->def->nsamples = n;
						data->def->bufsize  = n*data->def->nchans*sizeof(FLOAT32_T);
						response->def->bufsize = append(&response->buf, response->def->bufsize, data->def, sizeof(datadef_t));
						data->def->nsamples = MAXNUMSAMPLE;
						data->def->bufsize  = 0;
						for (j=0; j<n; j++) {
							response->def->bufsize = append(&response->buf, response->def->bufsize, (char*)(data->buf)+WRAP(datasel->begsample+j,MAXNUMSAMPLE)*(data->def->nchans)*sizeof(FLOAT32_T), (data->def->nchans)*sizeof(FLOAT32_T));
						}
						break;

					case DATATYPE_FLOAT64:
						data->def->nsamples = n;
						data->def->bufsize  = n*data->def->nchans*sizeof(FLOAT64_T);
						response->def->bufsize = append(&response->buf, response->def->bufsize, data->def, sizeof(datadef_t));
						data->def->nsamples = MAXNUMSAMPLE;
						data->def->bufsize  = 0;
						for (j=0; j<n; j++) {
							response->def->bufsize = append(&response->buf, response->def->bufsize, (char*)(data->buf)+WRAP(datasel->begsample+j,MAXNUMSAMPLE)*(data->def->nchans)*sizeof(FLOAT64_T), (data->def->nchans)*sizeof(FLOAT64_T));
						}
						break;

					default:
						fprintf(stderr, "dmarequest: unsupported data type (%d)\n", data->def->data_type);
						response->def->version = VERSION;
						response->def->command = GET_ERR;
						response->def->bufsize = 0;
				}
			}

			FREE(datasel);
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
			// determine the selection
			if (request->def->bufsize) {
				// the selection has been specified
				memcpy(eventsel, request->buf, sizeof(eventsel_t));
			}
			else {
				// determine a valid selection
				if (header->def->nevents>MAXNUMEVENT) {
					// the ringbuffer is completely full
					eventsel->begevent = header->def->nevents - MAXNUMEVENT;
					eventsel->endevent = header->def->nevents - 1;
				}
				else {
					// the ringbuffer is not yet completely full
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
				fprintf(stderr, "dmarequest: err1\n");
				response->def->version = VERSION;
				response->def->command = GET_ERR;
				response->def->bufsize = 0;
			}
			else if (eventsel->begevent >= header->def->nevents || eventsel->endevent >= header->def->nevents) {
				fprintf(stderr, "dmarequest: err2\n");
				response->def->version = VERSION;
				response->def->command = GET_ERR;
				response->def->bufsize = 0;
			}
			else if ((header->def->nevents-eventsel->begevent) > MAXNUMEVENT) {
				fprintf(stderr, "dmarequest: err3\n");
				response->def->version = VERSION;
				response->def->command = GET_ERR;
				response->def->bufsize = 0;
			}
			else {
				response->def->version = VERSION;
				response->def->command = GET_OK;
				response->def->bufsize = 0;

                // determine the number of events to return
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

		case GET_PRP:
			if (verbose>1) fprintf(stderr, "dmarequest: GET_PRP\n");
			if (property==NULL) {
				response->def->version = VERSION;
				response->def->command = GET_ERR;
				response->def->bufsize = 0;
				break;
			}

			pthread_mutex_lock(&mutexproperty);

			if (request->def->bufsize) {
				// determine whether the property exists in the buffer
				// the message payload should include a full property def and buf 
				// although the value_type and value_value will not be used here
				desired = (property_t*)malloc(sizeof(property));
				desired->def = (propertydef_t*)request->buf;
				if (desired->def->bufsize)
					desired->buf = (char*)request->buf+sizeof(propertydef_t);
				else
					desired->buf = NULL;
				n = find_property(desired);
				FREE(desired);

				if (n>=0) {
					response->def->version = VERSION;
					response->def->command = GET_OK;
					response->def->bufsize = 0;
					response->def->bufsize = append(&response->buf, response->def->bufsize, property[n].def, sizeof(propertydef_t));
					response->def->bufsize = append(&response->buf, response->def->bufsize, property[n].buf, property[n].def->bufsize);
				}
				else {
					response->def->version = VERSION;
					response->def->command = GET_ERR;
					response->def->bufsize = 0;
				}
			}
			else {
				// send all properties
				response->def->version = VERSION;
				response->def->command = GET_OK;
				response->def->bufsize = 0;
			    if (verbose>1) fprintf(stderr, "dmarequest: sending %d properties\n", thisproperty);
				for (n=0; n<thisproperty; n++) {
					if (verbose>1) print_propertydef(property[n].def);
					response->def->bufsize = append(&response->buf, response->def->bufsize, property[n].def, sizeof(propertydef_t));
					response->def->bufsize = append(&response->buf, response->def->bufsize, property[n].buf, property[n].def->bufsize);
				}
			}

			pthread_mutex_unlock(&mutexproperty);
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
			pthread_mutex_unlock(&mutexheader);
			pthread_mutex_unlock(&mutexdata);
			pthread_mutex_unlock(&mutexevent);
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

		case FLUSH_PRP:
			pthread_mutex_lock(&mutexproperty);
			if (property) {
				thisproperty = 0;
				for (i=0; i<MAXNUMPROPERTY; i++) {
					FREE(property[i].def);
					FREE(property[i].buf);
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
			pthread_mutex_unlock(&mutexproperty);
			break;

		default:
			fprintf(stderr, "dmarequest: unknown command\n");
	}

	if (verbose>0) fprintf(stderr, "dmarequest: thissample = %u, thisevent = %u\n", thissample, thisevent);

	/* everything went fine */
	return 0;
}

