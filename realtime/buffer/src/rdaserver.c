/*
 * Copyright (C) 2010 S. Klanke
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 *
 */

#include <rdaserver.h>
#include <fcntl.h>
#include <errno.h>

#define RDA_MAX_NUM_CLIENTS  32			/* On Windows, this number must be smaller than 64 */

const UINT8_T _rda_guid[16]={
	0x8E,0x45,0x58,0x43,0x96,0xC9,0x86,0x4C,0xAF,0x4A,0x98,0xBB,0xF6,0xC9,0x14,0x50
};

/** RDA is defined as a little-endian protocol. If this code runs on a big-endian
	computer, than the following flag will be set to True[=1] in rda_start_server,
	and all packets will be converted accordingly.
*/
static int _i_am_big_endian_ = 0;

/* return 0 on success, -1 on error (out of memory), item may not be NULL */
int rda_aux_alloc_item(rda_buffer_item_t *item, size_t size) {
	if (size > item->sizeAlloc) {
		void *data = malloc(size);
		if (data == NULL) return -1;
		if (item->data != NULL) free(item->data);
		item->data = data;
		item->sizeAlloc = size;
	}
	item->size = size;
	return 0;
}


/** Wait for new samples/events in the FieldTrip buffer
	@param ft_buffer  	FieldTrip connection
	@param current		Previous (already handled!) number of samples/events (must not be NULL!)
	@param result		Current number of samples/events (must not be NULL!)
	@param ms			Timeout in milliseconds
	@return 	0 on success, -1 on error
*/
int rda_aux_wait_dat(int ft_buffer, const samples_events_t *previous, samples_events_t *result, int ms) {
	int r;
	message_t req, *resp = NULL;
	messagedef_t msg_def;
	waitdef_t wait_def;
	
	wait_def.threshold = *previous;
	wait_def.milliseconds = ms;
	
	req.def = &msg_def;
	req.buf = &wait_def;
	msg_def.version = VERSION;
	msg_def.command = WAIT_DAT;
	msg_def.bufsize = sizeof(waitdef_t);
	
	r = clientrequest(ft_buffer, &req, &resp);
	if (r<0) return r;
	
	if (resp == NULL || resp->def == NULL || resp->buf == NULL || 
		resp->def->command != WAIT_OK || 
		resp->def->bufsize != sizeof(samples_events_t)) {
		printf("Weird things...\n");
		r = -1;
		goto cleanup;
	}
	memcpy(result, resp->buf, sizeof(samples_events_t));
cleanup:
	if (resp) {
		if (resp->buf) free(resp->buf);
		if (resp->def) free(resp->def);
		free(resp);
	}
	return r;
}


/** Prepares "start" RDA packet with channel names determined from the corresponding chunk (or empty)
	Also converts to little-endian if this machine is big endian
	@param item	Receives the packet (usually "startItem")
	@param hdr  Points to headerdef_t structure, will be filled, may not be NULL
	@return 0 on success, <0 on connection errors, >0 if out of memory 
*/
int rda_aux_get_hdr_prep_start(int ft_buffer, headerdef_t *hdr, rda_buffer_item_t *item) {
	const ft_chunk_t *chunk;
	rda_msg_start_t *R;
	char *str;
	double *dRes;
	const void *dResSource;
	size_t bytesTotal;
	int i,r,numExtraZeros,sizeOrgNames;
	message_t req, *resp = NULL;
	messagedef_t msg_def;	
	
	req.def = &msg_def;
	req.buf = NULL;
	msg_def.version = VERSION;
	msg_def.command = GET_HDR;
	msg_def.bufsize = 0;
	
	r = clientrequest(ft_buffer, &req, &resp);
	if (r<0) return r;
	
	if (resp == NULL || resp->def == NULL || resp->buf == NULL || resp->def->command != GET_OK || resp->def->bufsize < sizeof(headerdef_t)) {
		r = -1;
		goto cleanup;
	}
	
	memcpy(hdr, resp->buf, sizeof(headerdef_t));
	
	/* Ok, we have the basic header, now look for proper FT_CHUNK_RESOLUTIONS */
	chunk = find_chunk(resp->buf, sizeof(headerdef_t), resp->def->bufsize, FT_CHUNK_RESOLUTIONS);
	if (chunk != NULL && chunk->def.size == hdr->nchans*sizeof(double)) {
		dResSource = chunk->data;
		/* fine - we just need a memcpy later on */
	} else {
		dResSource = NULL;
		/* no suitable chunk found - set defaults later on */
	}
	
	/* Now see if we can find channel names */
	chunk = find_chunk(resp->buf, sizeof(headerdef_t), resp->def->bufsize, FT_CHUNK_CHANNEL_NAMES);
	if (chunk != NULL && chunk->def.size >= hdr->nchans) {
		/* The chunk seems ok - check whether we really have N 0-terminated string */
		int k,nz = 0;
		for (k = 0; k<chunk->def.size && nz<=hdr->nchans; k++) {
			if (chunk->data[k] == 0) nz++;
		}
		/* Okay, either k is at the end and we have nz<=N, or we have reached N=nz before the
		   end of the chunk. In both cases, it's safe to transmit the first 'k' bytes and
		   add (N-nz) trailing zeros */
		numExtraZeros = hdr->nchans - nz;
		sizeOrgNames = k;
	} else {
		sizeOrgNames = 0;
		numExtraZeros = hdr->nchans;
	}
	
	bytesTotal = sizeof(rda_msg_start_t) + hdr->nchans*(sizeof(double)) + sizeOrgNames + numExtraZeros;
	
	if (rda_aux_alloc_item(item, bytesTotal)) {
		r = 1;
		goto cleanup;
	}
	
	R = (rda_msg_start_t *) item->data;
	memcpy(R->hdr.guid, _rda_guid, sizeof(_rda_guid));
	R->hdr.nSize = bytesTotal;
	R->hdr.nType = RDA_START_MSG;
	R->nChannels = hdr->nchans;
	R->dSamplingInterval = 1.0/hdr->fsample;
	
	if (_i_am_big_endian_) {
		/* take care of hdr.nSize, hdr.nType, nChannels */
		ft_swap32(3, &(R->hdr.nSize)); 
		ft_swap64(1, &(R->dSamplingInterval));
	}
	
	/* R+1 points to first byte after header info */
	dRes = (double *) ((void *)(R+1)); 
	if (dResSource == NULL) {
		/* Fill with resolution = 1.0 -- we have nothing better */
		for (i=0;i<hdr->nchans;i++) dRes[i]=1.0;
	} else {
		memcpy(dRes, dResSource, hdr->nchans * sizeof(double));
	}
	/* swap byte order if necessary */
	if (_i_am_big_endian_) {
		ft_swap64(hdr->nchans, dRes);
	}
	
	/* Let 'str' point to first byte after the resolution values */
	str = (char *) ((void *)(dRes + hdr->nchans));
	if (sizeOrgNames > 0) {
		memcpy(str, chunk->data, sizeOrgNames);
	}
	for (i=0;i<numExtraZeros;i++) str[sizeOrgNames + i] = 0;
	/* done */
	r = 0;
cleanup:
	if (resp) {
		if (resp->def) free(resp->def);
		if (resp->buf) free(resp->buf);
		free(resp);
	}
	return r;
}	


void rda_aux_convert_to_float(UINT32_T N, void *dest, UINT32_T data_type, const void *src) {
	UINT32_T n;
	float *d = (float *) dest;
	switch(data_type) {
		case DATATYPE_CHAR:
		case DATATYPE_UINT8:
			{
				const UINT8_T *s = (const UINT8_T *) src;
				for (n=0;n<N;n++) d[n] = (float) s[n];
			}
			break;
		case DATATYPE_UINT16:
			{
				const UINT16_T *s = (const UINT16_T *) src;
				for (n=0;n<N;n++) d[n] = (float) s[n];
			}
			break;
		case DATATYPE_UINT32:
			{
				const UINT32_T *s = (const UINT32_T *) src;
				for (n=0;n<N;n++) d[n] = (float) s[n];
			}
			break;
		case DATATYPE_UINT64:
			{
				const UINT64_T *s = (const UINT64_T *) src;
				for (n=0;n<N;n++) d[n] = (float) s[n];
			}
			break;
		case DATATYPE_INT8:
			{
				const INT8_T *s = (const INT8_T *) src;
				for (n=0;n<N;n++) d[n] = (float) s[n];
			}
			break;
		case DATATYPE_INT16:
			{
				const INT16_T *s = (const INT16_T *) src;
				for (n=0;n<N;n++) d[n] = (float) s[n];
			}
			break;
		case DATATYPE_INT32:
			{
				const INT32_T *s = (const INT32_T *) src;
				for (n=0;n<N;n++) d[n] = (float) s[n];
			}
			break;
		case DATATYPE_INT64:
			{
				const INT64_T *s = (const INT64_T *) src;
				for (n=0;n<N;n++) d[n] = (float) s[n];
			}
			break;
		case DATATYPE_FLOAT32:
			memcpy(dest,src,N*sizeof(float));
			break;
		case DATATYPE_FLOAT64:
			{
				const double *s = (const double *) src;
				for (n=0;n<N;n++) d[n] = (float) s[n];
			}
			break;
	}
}



/* returns -1 on errors, or the number of markers if succesfully allocated and transformed */
int rda_aux_get_markers(int ft_buffer, const samples_events_t *last, const samples_events_t *cur, unsigned int sizeSamples, rda_buffer_item_t *item) {
	int numEvt = 0;
	message_t req, *resp = NULL;
	messagedef_t msg_def;
	eventsel_t evt_sel;
	size_t bytesHdrData = sizeof(rda_msg_data_t) + sizeSamples;
	
	if (cur->nevents<=last->nevents) {
		/* actually there are no markers present - just allocate enough space for the data */
		return rda_aux_alloc_item(item, bytesHdrData);	/* 0 on success - that's what we need here */
	}
	
	msg_def.version = VERSION;
	msg_def.command = GET_EVT;
	msg_def.bufsize = sizeof(evt_sel);
	evt_sel.begevent = last->nevents;
	evt_sel.endevent = cur->nevents-1;
	req.def = &msg_def;
	req.buf = &evt_sel;
	
	if (clientrequest(ft_buffer, &req, &resp) < 0) return -1;
	if (resp == NULL) return -1;
	if (resp->def == NULL || resp->buf == NULL) {
		numEvt = -1;
		goto cleanup;
	}
	
	if (resp->def->command == GET_OK) {
		size_t bytesMarkers = 0;
		int offset = 0;
		char *ptr;
		
		/* count the number of events, increase bytesTotal as required */
		while (offset + sizeof(eventdef_t) <= resp->def->bufsize) {
			eventdef_t *evdef = (eventdef_t *) ((char *)resp->buf + offset);
			offset += sizeof(eventdef_t) + evdef->bufsize;
			
			if (evdef->bufsize < evdef->type_numel*wordsize_from_type(evdef->type_type) + evdef->value_numel*wordsize_from_type(evdef->value_type)) {
				fprintf(stderr,"Invalid event received: Buffer to small for given value/type description\n");
				continue;	/* Skip to next event */
			}
			
			bytesMarkers += sizeof(rda_marker_t);
			if (evdef->type_type == DATATYPE_CHAR) {
				if (evdef->value_type == DATATYPE_CHAR) {
					/* Transform into TYPE\0VALUE\0 */
					bytesMarkers += evdef->type_numel + evdef->value_numel + 2;
				} else {
					/* Transform into TYPE\0-\0 */
					bytesMarkers += evdef->type_numel + 3;
				}
			} else {
				if (evdef->value_type == DATATYPE_CHAR) {
					/* Transform into FT\0VALUE\0 */
					bytesMarkers += evdef->value_numel + 4;
				} else {
					/* Transform into FT\0-\0 */
					bytesMarkers += 5; 
				}
			}
			numEvt++;
		}
		
		if (rda_aux_alloc_item(item, bytesHdrData + bytesMarkers)) {
			numEvt = -1;
			goto cleanup;
		}
		/* Ok, we've reserved enough space for header, data, and markers.
			Here we just fill in the markers */
		ptr = (char *) item->data + bytesHdrData;
		offset = 0;
		
		/* count the number of events, increase bytesTotal as required */
		while (offset + sizeof(eventdef_t) <= resp->def->bufsize) {
			eventdef_t *evdef = (eventdef_t *) ((char *)resp->buf + offset);
			char *evbuf = (char *)resp->buf + offset + sizeof(eventdef_t);
			rda_marker_t *marker = (rda_marker_t *) ptr;
			int i;
			
			offset += sizeof(eventdef_t) + evdef->bufsize;			
			
			if (evdef->bufsize < evdef->type_numel*wordsize_from_type(evdef->type_type) + evdef->value_numel*wordsize_from_type(evdef->value_type)) {
				continue; /* skip to next event */
			}
			
			marker->nChannel = -1; /* All channels, FieldTrip doesn't have this*/
			marker->nPosition = evdef->sample - (last->nsamples +1); /* relative to first sample in this block */
			marker->nPoints = evdef->duration;
			
			ptr += sizeof(rda_marker_t);
			
			if (evdef->type_type == DATATYPE_CHAR) {
				/* copy type */
				for (i=0;i<evdef->type_numel;i++) {
					*ptr++ = *evbuf++;
				}
			} else {
				*ptr++ = 'F';
				*ptr++ = 'T';
				/* Skip bytes in event->buffer */
				evbuf += evdef->type_numel;
			}
			*ptr++ = 0;				
			if (evdef->value_type == DATATYPE_CHAR) {
				/* copy value */
				for (i=0;i<evdef->value_numel;i++) {
					*ptr++ = *evbuf++;
				}
			} else {
				*ptr++ = '-';
			}
			/* add trailing 0, see how big the complete marker got */
			*ptr++ = 0;
			marker->nSize = (ptr - (char *) marker);
			
			if (_i_am_big_endian_) {
				/* convert the 4 int32's in the marker definition */
				ft_swap32(4, (void *) marker);
			}
		}		
	} 
	
cleanup:
	if (resp) {
		FREE(resp->buf);
		FREE(resp->def);
		free(resp);
	}
	return numEvt;
}


/** Get data from FT buffer and convert it into a single precision RDA packet
	@return  0 on success, -1 on connection errors, -2 on out of memory
*/
int rda_aux_get_float_data(int ft_buffer, unsigned int start, unsigned int end, unsigned int nChans, void *data) {
	int r;
	message_t req, *resp = NULL;
	messagedef_t msg_def;
	datasel_t dat_sel;
	
	msg_def.version = VERSION;
	msg_def.command = GET_DAT;
	msg_def.bufsize = sizeof(dat_sel);
	dat_sel.begsample = start;
	dat_sel.endsample = end;
	req.def = &msg_def;
	req.buf = &dat_sel;
	
	r = clientrequest(ft_buffer, &req, &resp);
	if (r<0) return r;
	
	if (resp == NULL) return -1;
	if (resp->def == NULL || resp->buf == NULL) {
		r = -1;
		goto cleanup;
	}
	if (resp->def->command == GET_OK) {
		datadef_t *ddef = (datadef_t *) resp->buf;
		unsigned int numTotal = ddef->nchans * ddef->nsamples;
		
		if (ddef->nchans != nChans || ddef->nsamples != (end-start+1)) {
			r = -1;
			goto cleanup;
		}
				
		/* convert the samples, note that (R+1) points to the first byte after the RDA data header
		   and (ddef+1) points to the first byte after the FieldTrip data header */
		rda_aux_convert_to_float(numTotal, data, ddef->data_type, (void *)(ddef+1));
		if (_i_am_big_endian_) ft_swap32(numTotal, data);
		
		/* done */
		r = 0;
	} else {
		r = -1;
	}
cleanup:
	if (resp) {
		FREE(resp->buf);
		FREE(resp->def);
		free(resp);
	}
	return r;
}

/** Get 16-bit data from FT buffer and convert it into a single precision RDA packet
	@return  0 on success, -1 on connection errors, -2 on out of memory
*/
int rda_aux_get_int16_data(int ft_buffer, unsigned int start, unsigned int end, unsigned int nChans, void *data) {
	int r;
	message_t req, *resp = NULL;
	messagedef_t msg_def;
	datasel_t dat_sel;
	
	msg_def.version = VERSION;
	msg_def.command = GET_DAT;
	msg_def.bufsize = sizeof(dat_sel);
	dat_sel.begsample = start;
	dat_sel.endsample = end;
	req.def = &msg_def;
	req.buf = &dat_sel;
	
	r = clientrequest(ft_buffer, &req, &resp);
	if (r<0) return r;
	
	if (resp == NULL) return -1;
	if (resp->def == NULL || resp->buf == NULL) {
		r = -1;
		goto cleanup;
	}
	if (resp->def->command == GET_OK) {
		datadef_t *ddef = (datadef_t *) resp->buf;
		unsigned int numTotal = ddef->nchans * ddef->nsamples;
		
		if (((ddef->data_type != DATATYPE_INT16) && (ddef->data_type != DATATYPE_UINT16)) || (ddef->nchans != nChans)) {
			r = -1;
			goto cleanup;
		}
		
		if (_i_am_big_endian_) {
			/* copy + swap the 16 bit samples */
			unsigned char *source = (unsigned char *) (ddef+1);
			unsigned char *dest = data;
			unsigned int i;
			
			for (i=0;i<numTotal;i++) {
				dest[0] = source[1];
				dest[1] = source[0];
				dest+=2;
				source+=2;
			}
		} else {
			/* copy the samples, note that (ddef+1) points to the first byte after the FieldTrip data header */
			memcpy(data, (void *)(ddef+1), numTotal*sizeof(INT16_T));
		}
		/* done */
		r = 0;
	} else {
		r = -1;
	}
cleanup:
	if (resp) {
		FREE(resp->buf);
		FREE(resp->def);
		free(resp);
	}
	return r;
}


/** Thread function of the RDA server. Can handle multiple clients in parallel.
	TODO: Fieldtrip events are not transmitted as RDA markers, yet 
*/
void *_rdaserver_thread(void *arg) {
	/* 'select' cannot handle more than 64 elements on Windows, but this
		should really be enough for all practical purposes. Depending on
		the sampling rate and number of channels, you would probably hit
		other performance boundaries first.
	*/
	rda_client_job_t clients[RDA_MAX_NUM_CLIENTS];
	rda_buffer_item_t startItem = {NULL,0,0,0,NULL}; 
	rda_buffer_item_t *firstDataItem = NULL;
	
	headerdef_t ftHdr;
	int i,typeOk;
	int ftTimeout = 0;
	unsigned int numBlock = 0;
	
	samples_events_t lastNum = {0,0}, curNum;
	
	fd_set readSet, writeSet;
	int fdMax = 1;				 /* not really used on Windows */	
	struct timeval tv;
	rda_server_ctrl_t *SC = (rda_server_ctrl_t *) arg;
	int numClients = 0;

	
	if (SC==NULL) return NULL;
	
	SC->is_running = 1;
	/* Set typeOk flag to 1 for floats, 0 for int16 
		(in the latter case we need to check the FT header first)
	*/
	typeOk = SC->use16bit ? 0 : 1;
		
	/* wait up to 10 ms for new events */
	tv.tv_sec = 0; 
	tv.tv_usec = 10000;
	
	#ifndef PLATFORM_WIN32
	fdMax = SC->server_socket;
	#endif
	
	while (!SC->should_exit) {
		int sel;
		
		/* Prepare read (also error!) and write sets: clear them first */
		FD_ZERO(&readSet);
		FD_ZERO(&writeSet);
		/* Add server socket to read set, but only if we can actually handle more clients */
		if (numClients < RDA_MAX_NUM_CLIENTS) {
			FD_SET(SC->server_socket, &readSet);
		}
		for (i=0;i<numClients;i++) {
			/* Every client is listened to (for disconnection!) */
			FD_SET(clients[i].sock, &readSet);
			/* But only clients with a pending start/data are added to the write set */
			if (clients[i].item != NULL) FD_SET(clients[i].sock, &writeSet);
		}

		/* Check server and client sockets for possible read and write operations */
		sel = select(fdMax + 1, &readSet, &writeSet, NULL, &tv);
		if (sel == -1) {
			perror("rdaserver_thread -- select");
			continue; /* TODO: think about stopping operation */
		}
		
		/* In case of a timeout (=no network events): check the FieldTrip buffer for header + data */
		if (sel == 0) {
			if (startItem.data == NULL) {
				if (rda_aux_get_hdr_prep_start(SC->ft_buffer, &ftHdr, &startItem)) continue;
				/* Now the start item should have proper size/sizeAlloc/data fields */
				if (SC->use16bit) {
					typeOk = (ftHdr.data_type == DATATYPE_INT16) || (ftHdr.data_type == DATATYPE_UINT16);
				}
				if (SC->verbosity > 4) {
					printf("Picked up FT header for first time!\n"); 
				}
				if (typeOk) {
					/* Add the start item to every client */
					for (i=0;i<numClients;i++) {
						clients[i].item = &startItem;
						clients[i].written = 0;
					}
					lastNum.nsamples = ftHdr.nsamples;
					lastNum.nevents = ftHdr.nevents;
				}
				continue;
			}
			/* Check for new data (samples / events) */
			if (rda_aux_wait_dat(SC->ft_buffer, &lastNum, &curNum, ftTimeout)) continue;
			
			/* TODO: check if new numbers are smaller, which indicates a rewritten header */
			if (!typeOk) continue;
				
			if (curNum.nsamples > lastNum.nsamples || curNum.nevents > lastNum.nevents) {
				/* There's new data to stream out */
				rda_buffer_item_t *item, *last;
				int newItem = 0;
				int numEvts;
				size_t sizeSamples;
				
				/* Three cases:
					a) first data packet ever (firstDataItem = NULL)
						-> create new item and have firstDataItem point to it
					b) firstDataItem has refCount=0
						-> put data packet into that item (possible recycle memory)
					c) otherwise
						-> create new item and link it into the back of the list
				*/
				if (firstDataItem == NULL) {
					/* case a */
					newItem = 1;
				} else if (firstDataItem->refCount > 0) {
					/* case c */
					newItem = 1;
					last = firstDataItem;
					while (last->next != NULL) last=last->next;
				}
				
				if (newItem) {
					item = (rda_buffer_item_t *) malloc(sizeof(rda_buffer_item_t));
					if (item == NULL) {
						fprintf(stderr,"rdaserver_thread: Out of memory\n");
						break;
					}
					item->next = NULL;
					item->refCount = 0;
					item->sizeAlloc = 0;
					item->data = NULL;
				} else {
					item = firstDataItem;
				}
				
				if (SC->verbosity > 5) {
					printf("Trying to get data [%i ; %i], events [%i;%i]\n", lastNum.nsamples, curNum.nsamples, lastNum.nevents, curNum.nevents); 
				}
				
				sizeSamples = (SC->use16bit) ? sizeof(INT16_T) : sizeof(float);
				sizeSamples *= (curNum.nsamples - lastNum.nsamples) * ftHdr.nchans;
				
				numEvts = rda_aux_get_markers(SC->ft_buffer, &lastNum, &curNum, sizeSamples, item);
				if (numEvts >= 0) {
					rda_msg_data_t *R = (rda_msg_data_t *) item->data;
					void *sample_data = (void *)(R+1);
					memcpy(R->hdr.guid, _rda_guid, sizeof(_rda_guid));
					R->hdr.nSize = item->size;
					R->nBlock = numBlock; 
					R->nPoints = curNum.nsamples - lastNum.nsamples;
					R->nMarkers = numEvts;
					if (_i_am_big_endian_) {
						/* take care of hdr.nSize, hdr.nType, nBlocks, nPoints, nMarkers */
						ft_swap32(5, &(R->hdr.nSize)); 
						/* samples are swapped in rda_aux_get_***_data */
						/* markers were swapped in rda_aux_get_markers */
					}
					if (curNum.nsamples > lastNum.nsamples) {
						int err;
						if (SC->use16bit) {
							R->hdr.nType = RDA_INT_MSG;
							err = rda_aux_get_int16_data(SC->ft_buffer, lastNum.nsamples, curNum.nsamples-1, ftHdr.nchans, sample_data);
						} else {
							R->hdr.nType = RDA_FLOAT_MSG;
							err = rda_aux_get_float_data(SC->ft_buffer, lastNum.nsamples, curNum.nsamples-1, ftHdr.nchans, sample_data);
						}
						if (err!=0) {
							fprintf(stderr, "Error while getting data from FieldTrip buffer!\n");
							numEvts = -1;
						}
					}
				} else {
					fprintf(stderr, "Error in rda_aux_get_markers\n");
				}
				
				if (numEvts < 0) {
					/* get rid of the item again */
					if (newItem) {
						if (item->data != NULL) free(item->data);
						free(item);
					}
					continue;
				}
					
				
				if (newItem) {
					if (firstDataItem == NULL) {
						firstDataItem = item;
					} else {
						last->next = item;
					}
				}
				
				for (i=0;i<numClients;i++) {
					if (clients[i].item == NULL) {
						if (SC->verbosity>5) {
							printf("Adding new job for client %i\n", clients[i].sock);
						}
						clients[i].item = item;
						item->refCount++;
					}
				}
				lastNum = curNum;
				numBlock++;
			}
		}
		
		/* Check if we have a new client connection */
		if (sel>0 && FD_ISSET(SC->server_socket, &readSet)) {
			struct sockaddr_in sa;
			int size_sa = sizeof(sa);
			SOCKET newSock;
	
			newSock = accept(SC->server_socket, (struct sockaddr *)&sa, &size_sa);
			if (newSock == INVALID_SOCKET) {
				perror("rda_server_thread - accept");
			} else {
				if (SC->verbosity > 0) fprintf(stderr, "rdaserver_thread: opened connection to client (%i)\n",newSock);
				
				clients[numClients].sock = newSock;
				clients[numClients].item = startItem.data == NULL ? NULL : &startItem;
				clients[numClients].written = 0;
				numClients++;
				pthread_mutex_lock(&SC->mutex);
				SC->num_clients = numClients;
				pthread_mutex_unlock(&SC->mutex);
				#ifndef PLATFORM_WIN32
				if (newSock > fdMax) fdMax = newSock;
				#endif
			}
			--sel;
		}
		/* Check for sockets that are ready to be written to */
		for (i=0 ; i<numClients && sel>0 ; i++) {
			int sent;
			rda_client_job_t *C = &clients[i];

			if (!FD_ISSET(C->sock, &writeSet)) continue;
						
			/* printf("Sending out data (%i bytes)...\n", C->item->size - C->written); */
			sent = send(C->sock, (char *) C->item->data + C->written, C->item->size - C->written, 0);
			/* printf("%i bytes\n", sent); */
			if (sent > 0) {
				C->written += sent;
				if (C->written == C->item->size) {
					/* printf("Done with this packet on client %i\n", C->sock); */
					C->item->refCount--;
					C->written = 0;
					C->item = C->item->next;
					if (C->item != NULL) {
						C->item->refCount++;
					}
				}
			} else {
				if (SC->verbosity > 0) {
					fprintf(stderr, "rdaserver_thread: write error on socket %i\n",C->sock);
				}
				C->item = NULL;
			}
			--sel;
		}
		/* Check for sockets on which we can read (=> read 0 bytes means closure ) */
		for (i=0 ; i<numClients && sel>0 ; i++) {
			char dummy[1024];
			int len;
			if (!FD_ISSET(clients[i].sock, &readSet)) continue;
			
			len = recv(clients[i].sock, dummy, sizeof(dummy), 0);
			
			if (len>0) continue;	/* clients are not supposed to write, and if they do, we ignore it */
			if (len<0) {
				/* close this client */
				if (SC->verbosity > 0) fprintf(stderr, "rdaserver_thread: lost connection to client (%i)\n", clients[i].sock); 
			} else {
				if (SC->verbosity > 0) fprintf(stderr, "rdaserver_thread: client (%i) closed the connection\n", clients[i].sock); 
			}
			
			closesocket(clients[i].sock);
			
			if (clients[i].item != NULL) {
				/* apparently we lost this connection during transmission */
				clients[i].item->refCount--;
			}
			#ifndef PLATFORM_WIN32
			fdMax = SC->server_socket;
			for (i=0;i<numClients;i++) {
				if (clients[i].sock > fdMax) fdMax = clients[i].sock;
			}
			#endif
			/* Remove clients from list by moving the last one in its place */
			if (numClients > 1) {
				clients[i] = clients[numClients-1];
			}
			numClients--;
			pthread_mutex_lock(&SC->mutex);
			SC->num_clients = numClients;
			pthread_mutex_unlock(&SC->mutex);
			--sel;
		}
		
		/* just for debugging - this will be removed */
		if (0) {
			rda_buffer_item_t *item = firstDataItem;
			while (item!=NULL) {
				printf("Item at 0x%lX has refcount %i and points at 0x%lX\n", (long) (void *) item,  item->refCount, (long) (void *) item->next);
				item = item->next;
			}
		}
		
		/* Done with the network-specific stuff, now clean up the data items */
		if (firstDataItem!=NULL) {		
			/* We always keep one item to put the next data in, but also to be able to send
			   out data as soon as a client connects
			*/
			while (firstDataItem->refCount == 0 && firstDataItem->next != NULL) {
				if (firstDataItem->data != NULL) free(firstDataItem->data);
				/* printf("Shifting data item list...\n"); */
				firstDataItem = firstDataItem->next;
			}
			startItem.next = firstDataItem;
		}
	}		
	/* shutdown clients */
	for (i=0;i<numClients;i++) {
		#ifdef PLATFORM_WIN32
		shutdown(clients[i].sock, SD_BOTH);
		closesocket(clients[i].sock);
		#else
		shutdown(clients[i].sock, SHUT_RDWR);
		close(clients[i].sock);
		#endif
	}
	/* free memory pointed to by startItem ... */
	if (startItem.data) free(startItem.data);
	/* ... and firstDataItem (including following list elements) */
	while (firstDataItem!=NULL) {
		rda_buffer_item_t *item = firstDataItem;
		firstDataItem = item->next;
		if (item->data) free(item->data);
		free(item);
	}
	SC->is_running = 0;
	return NULL;
}
 
/* see header file for documentation */
rda_server_ctrl_t *rda_start_server(int ft_buffer, int use16bit, int port, int *errval) {
	rda_server_ctrl_t *SC = NULL;
	SOCKET s = INVALID_SOCKET;
	struct sockaddr_in sa;
	unsigned long optval;
	int interr = FT_ERR_SOCKET;  /* if things go wrong here, it's most often because of socket errors */
	UINT16_T testEndian = 0x0100;	/* this will be [0x00,0x01] on little-endian, [0x01,0x00] on big-endian */
		
#ifdef PLATFORM_WIN32
	WSADATA wsa;
 	if(WSAStartup(MAKEWORD(1, 1), &wsa))
	{
		fprintf(stderr, "tcpserver: cannot start sockets\n");
		goto cleanup;
	}
#endif

	if (*((UINT8_T *) &testEndian)) {
		_i_am_big_endian_ = 1;
		/*
		printf("Running on big-endian machine - conversion enabled\n");
		*/
	}

	/* allocate the control structure */
	SC = (rda_server_ctrl_t *) malloc(sizeof(rda_server_ctrl_t));
	if (SC == NULL) {
		fprintf(stderr,"start_rda_server: out of memory\n");
		interr = FT_ERR_OUT_OF_MEM;
		goto cleanup;
	}
	
	/* create TCP socket */
	interr = FT_ERR_SOCKET;
	s = socket(PF_INET, SOCK_STREAM, 0);
	if (s == INVALID_SOCKET) {
		perror("start_rda_server socket");
		goto cleanup;
	}
	
	/* prevent "bind: address already in use" */
	optval = 1;
	if (setsockopt(s, SOL_SOCKET, SO_REUSEADDR, (const char*)&optval, sizeof(optval)) < 0) {
		perror("start_rda_server setsockopt");
		/* not really critical - we go on */
	}
	
	/* check if user selected default port and select according to type */
	if (port == 0) {
		port = use16bit ? 51234 : 51244;
	}
	
	/* bind socket to the specified port (all interfaces) */
	bzero(&sa, sizeof(sa));
	sa.sin_family = AF_INET;
	sa.sin_port   = htons(port);
	sa.sin_addr.s_addr = htonl(INADDR_ANY);

	if (bind(s, (struct sockaddr *)&sa, sizeof(sa)) != 0) {
		perror("start_rda_server bind");
		goto cleanup;
	}
	
	/* place the socket in non-blocking mode */
#ifdef PLATFORM_WIN32
	optval = 1;
	if (ioctlsocket(s, FIONBIO, &optval) != 0) {
		fprintf(stderr,"start_rda_server: could not set non-blocking mode\n");
		goto cleanup;
	}
#else
	optval = fcntl(s, F_GETFL, NULL);
	optval = optval | O_NONBLOCK;
	if (fcntl(s, F_SETFL, optval)<0) {
		perror("start_rda_server fcntl");
		goto cleanup;
	}
#endif
	
	if (listen(s, BACKLOG)<0) {
		perror("tcpserver listen");
		goto cleanup;
	}
	
	/* set some control variables */
	SC->num_clients = 0;
	SC->should_exit = 0;
	SC->server_socket = s;
	SC->ft_buffer = ft_buffer;
	SC->use16bit = use16bit;
	SC->verbosity = 10; /* TODO: specify proper values */

	/* if things go wrong after this, it's because of pthread issues */
	interr = FT_ERR_THREADING;
	
	/* create the mutex */
	if (pthread_mutex_init(&SC->mutex, NULL) != 0) {
		fprintf(stderr,"start_rda_server: mutex could not be initialised\n");
		goto cleanup;
	}
	
	/* create thread with default attributes, select thread function depending on use16bit flag */
	if (pthread_create(&SC->thread, NULL, _rdaserver_thread, SC) == 0) {
		/* everything went fine - thread should be running now */
		if (errval != NULL) *errval = FT_NO_ERROR;
		return SC;
	}
	
	pthread_mutex_destroy(&SC->mutex);
cleanup:
	if (errval!=NULL) *errval = interr;
	if (SC != NULL) free(SC);
	if (s != INVALID_SOCKET) {
		#ifdef PLATFORM_WIN32
		shutdown(s, SD_BOTH);
		#else
		shutdown(s, SHUT_RDWR);
		#endif
		closesocket(s);
	}
	return NULL;
}	

/* stops the server and free's the given control structure */
/* see header file for documentation */
int rda_stop_server(rda_server_ctrl_t *SC) {
	if (SC==NULL) return -1;
	
	pthread_mutex_lock(&SC->mutex);
	SC->should_exit = 1;
	pthread_mutex_unlock(&SC->mutex);
	
	pthread_join(SC->thread, NULL);
	pthread_detach(SC->thread);
	pthread_mutex_destroy(&SC->mutex);
	
	closesocket(SC->server_socket);
	free(SC);
	return 0;
}

/* see header file for documentation */
int rda_get_num_clients(rda_server_ctrl_t *SC) {
	int nc;
	if (SC==NULL) return 0;
	pthread_mutex_lock(&SC->mutex);
	nc = SC->num_clients;
	pthread_mutex_unlock(&SC->mutex);
	return nc;
}
