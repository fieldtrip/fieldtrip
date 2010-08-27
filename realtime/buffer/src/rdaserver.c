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


const UINT8_T _rda_guid[16]={
	0x8E,0x45,0x58,0x43,0x96,0xC9,0x86,0x4C,0xAF,0x4A,0x98,0xBB,0xF6,0xC9,0x14,0x50
};

/** RDA is defined as a little-endian protocol. If this code runs on a big-endian
	computer, than the following flag will be set to True[=1] in rda_start_server,
	and all packets will be converted accordingly.
*/
static int _i_am_big_endian_ = 0;

/* returns new item or NULL on failure */
rda_buffer_item_t *rda_aux_alloc_item(size_t size) {
	rda_buffer_item_t *item = (rda_buffer_item_t *) malloc(sizeof(rda_buffer_item_t));
	if (item==NULL) return NULL;
	
	item->data = malloc(size);
	if (item->data == NULL) {
		free(item);
		return NULL;
	}
	item->size = size;
	item->next = NULL;
	item->refCount = 0;
	return item;
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
	
	if (r<0) {
		*result = *previous;
		return r;
	}
	
	if (resp == NULL || resp->def == NULL || resp->buf == NULL || 
		resp->def->command != WAIT_OK || 
		resp->def->bufsize != sizeof(samples_events_t)) {
		printf("Bad response from WAIT_DAT call\n");
		r = -1;
		*result = *previous;
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
	@param hdr  Points to headerdef_t structure, will be filled, may not be NULL
	@return created start item, or NULL on error (connection / out of memory)
*/
rda_buffer_item_t *rda_aux_get_hdr_prep_start(int ft_buffer, headerdef_t *hdr) {
	rda_buffer_item_t *item = NULL;
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
	if (r<0 || resp == NULL || resp->def == NULL) {
		goto cleanup;
	}
	
	if (resp->def->command != GET_OK || resp->def->bufsize < sizeof(headerdef_t) || resp->buf == NULL) {
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
		/* The chunk seems ok - check whether we really have N (0-terminated) strings */
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
		
	item = rda_aux_alloc_item(bytesTotal);
	if (item == NULL) goto cleanup;
		
	R = (rda_msg_start_t *) item->data;
	memcpy(R->hdr.guid, _rda_guid, sizeof(_rda_guid));
	R->hdr.nSize = bytesTotal;
	R->hdr.nType = RDA_START_MSG;
	R->nChannels = hdr->nchans;
	R->dSamplingInterval = 1.0e6/(double) hdr->fsample;	/* should be in microseconds */
	
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
cleanup:
	if (resp) {
		if (resp->def) free(resp->def);
		if (resp->buf) free(resp->buf);
		free(resp);
	}
	return item;
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



/** Retrieves samples and markers and returns them in as a new 'item', or NULL on errors.
*/
rda_buffer_item_t *rda_aux_get_samples_and_markers(int ft_buffer, const samples_events_t *last, const samples_events_t *cur, int numBlock, int use16bit) {
	rda_buffer_item_t *item = NULL;
	int numEvt = 0,numChans = 0,numSmp = 0;
	message_t req, *respSmp = NULL, *respEvt = NULL;
	messagedef_t msg_def;
	datadef_t *ddef;
	size_t bytesSamples = 0, bytesMarkers = 0;
	
	/* First, try to grab the samples */
	if (cur->nsamples > last->nsamples) {
		datasel_t dat_sel;
		
		msg_def.version = VERSION;
		msg_def.command = GET_DAT;
		msg_def.bufsize = sizeof(dat_sel);
		dat_sel.begsample = last->nsamples;
		dat_sel.endsample = cur->nsamples-1;
		req.def = &msg_def;
		req.buf = &dat_sel;
	
		if (clientrequest(ft_buffer, &req, &respSmp)<0) {
			goto cleanup;
		}
		
		numSmp = cur->nsamples - last->nsamples;
		
		if (respSmp == NULL || respSmp->def == NULL || respSmp->buf == NULL || respSmp->def->command != GET_OK) {
			goto cleanup;
		} else {
			ddef = (datadef_t *) respSmp->buf;
			
			if (ddef->nsamples != numSmp) goto cleanup;
			
			numChans = ddef->nchans;
			bytesSamples  = (use16bit ? sizeof(INT16_T) : sizeof(float)) * numSmp * numChans;
		}
	} 
	
	/* Now, try to grab the markers */
	if (cur->nevents > last->nevents) {
		eventsel_t evt_sel;
		int offset = 0;
		
		msg_def.version = VERSION;
		msg_def.command = GET_EVT;
		msg_def.bufsize = sizeof(evt_sel);
		evt_sel.begevent = last->nevents;
		evt_sel.endevent = cur->nevents-1;
		req.def = &msg_def;
		req.buf = &evt_sel;
	
		if (clientrequest(ft_buffer, &req, &respEvt) < 0) {
			goto cleanup;
		}
		if (respEvt == NULL || respEvt->def == NULL || respEvt->buf == NULL || respEvt->def->command != GET_OK) {
			goto cleanup;
		}
		
		/* count the number of events, increase bytesTotal as required */
		while (offset + sizeof(eventdef_t) <= respEvt->def->bufsize) {
			eventdef_t *evdef = (eventdef_t *) ((char *)respEvt->buf + offset);
			offset += sizeof(eventdef_t) + evdef->bufsize;
			
			if (evdef->bufsize < evdef->type_numel*wordsize_from_type(evdef->type_type) + evdef->value_numel*wordsize_from_type(evdef->value_type)) {
				fprintf(stderr,"Invalid event received: Buffer to small for given value/type description\n");
				continue;	/* Skip to next event */
			}
			
			bytesMarkers += sizeof(rda_marker_t);
			if (evdef->type_type == DATATYPE_CHAR) {
				if (evdef->value_type == DATATYPE_CHAR) {
					/* Transform into TYPE:VALUE\0 */
					bytesMarkers += evdef->type_numel + evdef->value_numel + 2;
				} else {
					/* Transform into TYPE:-\0 */
					bytesMarkers += evdef->type_numel + 3;
				}
			} else {
				if (evdef->value_type == DATATYPE_CHAR) {
					/* Transform into FT:VALUE\0 */
					bytesMarkers += evdef->value_numel + 4;
				} else {
					/* Transform into FT:-\0 */
					bytesMarkers += 5; 
				}
			}
			numEvt++;
		}
	}
	
	/* Now, allocate an item with enough space for both samples and markers */
	item = rda_aux_alloc_item(sizeof(rda_msg_data_t) + bytesSamples + bytesMarkers);
	if (item == NULL) {
		fprintf(stderr, "Out of memory\n");
		goto cleanup;
	}
	
	item->blockNumber = numBlock;
	
	/* Okay, we've got the samples in respSmp, events in respEvt, and a big enough 'item'.
		First fill in the header.
	*/
	{
		rda_msg_data_t *R = (rda_msg_data_t *) item->data;

		memcpy(R->hdr.guid, _rda_guid, sizeof(_rda_guid));
		R->hdr.nType = use16bit ? RDA_INT_MSG : RDA_FLOAT_MSG;
		R->hdr.nSize = item->size;
		R->nBlock = numBlock; 
		R->nPoints = numSmp;
		R->nMarkers = numEvt;
		if (_i_am_big_endian_) {
			/* take care of hdr.nSize, hdr.nType, nBlocks, nPoints, nMarkers */
			ft_swap32(5, &(R->hdr.nSize)); 
		}
	}
	
	/* Now, fill in the samples, possibly using conversion */
	if (numSmp > 0) {
		char *dataDest = ((char *) item->data + sizeof(rda_msg_data_t));
		char *dataSrc  = ((char *) respSmp->buf + sizeof(datadef_t));
		int numTotal = numSmp * numChans;
		
		if (use16bit) {
			if (_i_am_big_endian_) {
				/* copy + swap the 16 bit samples */
				int i;
				for (i=0;i<numTotal;i++) {
					dataDest[0] = dataSrc[1];
					dataDest[1] = dataSrc[0];
					dataDest+=2;
					dataSrc+=2;
				}
			} else {
				/* just copy the samples */
				memcpy(dataDest, dataSrc, bytesSamples);
			}
		} else {
			rda_aux_convert_to_float(numTotal, dataDest, ddef->data_type, dataSrc);
			if (_i_am_big_endian_) ft_swap32(numTotal, dataDest);
		}
	}
	
	/* Finally, fill in the events */
	if (numEvt>0) {
		char *ptr = (char *) item->data + sizeof(rda_msg_data_t) + bytesSamples;
		int offset = 0;
		
		/* count the number of events, increase bytesTotal as required */
		while (offset + sizeof(eventdef_t) <= respEvt->def->bufsize) {
			eventdef_t *evdef = (eventdef_t *) ((char *)respEvt->buf + offset);
			char *evbuf = (char *)respEvt->buf + offset + sizeof(eventdef_t);
			rda_marker_t *marker = (rda_marker_t *) ptr;
			int i, markerPos;
			
			offset += sizeof(eventdef_t) + evdef->bufsize;			
			
			if (evdef->bufsize < evdef->type_numel*wordsize_from_type(evdef->type_type) + evdef->value_numel*wordsize_from_type(evdef->value_type)) {
				continue; /* skip to next event */
			}
			
			markerPos = evdef->sample - (last->nsamples +1); /* relative to first sample in this block */
			marker->nPosition = (markerPos > 0) ? markerPos : 0;  /* needs to be unsigned! */
			marker->nChannel  = -1; /* All channels, FieldTrip doesn't have this*/
			marker->nPoints   = evdef->duration;
			
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
			*ptr++ = ':';				
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
	/* Done */
cleanup:
	if (respSmp) {
		FREE(respSmp->buf);
		FREE(respSmp->def);
		free(respSmp);
	}
	if (respEvt) {
		FREE(respEvt->buf);
		FREE(respEvt->def);
		free(respEvt);
	}
	return item;
}

/** Check for sockets that are ready to be written to, and write out as much of the corresponding job 
    as possible.
	@param 	remSelect	Number of remaining sockets to deal with
	@param	writeSet	The socket set to check for writable clients
	@param	numClients	The number of elements in the 'clients' array
	@param	clients		Array describing the clients and their jobs
	@param	verbosity	Determines how much status/error message to print
	@return	the remaining number of sockets to deal with (non-write operations)
*/
int rda_aux_check_writing(int remSelect, const fd_set *writeSet, int numClients, rda_client_job_t *clients, int verbosity) {
	int i;
	
	for (i=0 ; i<numClients && remSelect>0 ; i++) {
		int sent;
		rda_client_job_t *C = &clients[i];

		if (!FD_ISSET(C->sock, writeSet)) continue;	/* skip to next client in set */
						
		/* printf("Sending out data (%i bytes)...\n", C->item->size - C->written); */
		sent = send(C->sock, (char *) C->item->data + C->written, C->item->size - C->written, 0);
		
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
			if (verbosity>0) 
				fprintf(stderr, "rdaserver check_writing: write error on socket %i\n",C->sock);
			C->item = NULL;	/* do not attempt to send more data, will probably be closed later on */
		}
		--remSelect;
	}
	return remSelect;
}


/** Check for sockets that are have been closed from the other end, and remove the corresponding client.
	@param 	remSelect	Number of remaining sockets to deal with
	@param	readSet		The socket set to check for readable clients (read 0 bytes => close)
	@param	numClients	The number of elements in the 'clients' array
	@param	clients		Array describing the clients and their jobs
	@param	verbosity	Determines how much status/error message to print
	@return	the remaining number of clients
*/
int rda_aux_check_closure(int remSelect, const fd_set *readSet, int numClients, rda_client_job_t *clients, int verbosity) {
	int i = 0;	/* index of client we're looking at */
	
	while (i<numClients && remSelect>0) {
		char dummy[1024];
		int len;
		if (!FD_ISSET(clients[i].sock, readSet)) {
			++i;
			continue;
		}
			
		len = recv(clients[i].sock, dummy, sizeof(dummy), 0);
			
		if (len>0) {
			/* clients are not supposed to write, but if they do, we ignore it */
			++i;			
			continue;	
		}
		if (len<0) {
			/* close this client */
			if (verbosity > 0) fprintf(stderr, "rdaserver_thread: lost connection to client (%i)\n", clients[i].sock); 
		} else {	/* len == 0: smooth close from remote side */
			if (verbosity > 0) fprintf(stderr, "rdaserver_thread: client (%i) closed the connection\n", clients[i].sock); 
		}
			
		closesocket(clients[i].sock);
			
		if (clients[i].item != NULL) {
			/* apparently we lost this connection during transmission */
			clients[i].item->refCount--;
		}
		/* Remove clients from list by moving the last one in its place */
		if (i < numClients - 1) {
			clients[i] = clients[numClients-1];
		}
		--numClients;
		--remSelect;
		/* no increase of i here */
	}
	return numClients;
}


/** Check for sockets that are lagging behind significantly, and close them.
	@param 	minBlock    Threshold on block number of keeping the client alive
	@param	numClients	The number of elements in the 'clients' array
	@param	clients		Array describing the clients and their jobs
	@param	verbosity	Determines how much status/error message to print
	@return	the remaining number of clients
*/
int rda_aux_check_slow_clients(int minBlock, int numClients, rda_client_job_t *clients, int verbosity) {
	int i = 0;	/* index of client we're looking at */
	
	while (i<numClients) {
		if (clients[i].item != NULL && clients[i].item->blockNumber < minBlock) {
			if (verbosity > 0) {
				fprintf(stderr, "rdaserver_thread: disconnecting too slow client (%i)\n", clients[i].sock); 
			}
			/* Remove clients from list by moving the last one in its place */
			closesocket(clients[i].sock);
			if (i < numClients - 1) {
				clients[i] = clients[numClients-1];
			}
			clients[i].item->refCount--;
			--numClients;
			/* no increase of i here */
		} else {
			++i;
		}
	}
	return numClients;
}


/** Thread function of the RDA server. Can handle multiple clients in parallel.

	Incoming data (from the FieldTrip buffer) is first converted to "items" in a form that RDA clients expect,
	and a linked list of these items is kept. Each item has a reference count that determines on how many client
	sockets it's currently being streamed out. At the end of the main server loop, the linked list is inspected
	for unreferenced items at the start of the list, and those are removed.
	
	On top of the data items, also a "start item" is kept that contains the header information (RDA start packet).
	This is necessary for being able to write out the header information to newly connecting clients.
	
	Initially, both the "start item" and the "first data item" are empty, indicating that no header information
	and data/events	have been read yet.
*/
void *_rdaserver_thread(void *arg) {
	rda_server_ctrl_t *SC = (rda_server_ctrl_t *) arg;		/* our control structure */
	rda_client_job_t clients[RDA_MAX_NUM_CLIENTS];			/* list of clients and their current jobs */
	rda_buffer_item_t *startItem = NULL;				 	/* item containing start packet (header info) */
	rda_buffer_item_t *firstDataItem = NULL;				/* first item in list of (mostly) data packets */
	rda_buffer_item_t *latestItem = NULL;
	
	headerdef_t ftHdr;			/* contains header information */
	int i,typeOk;				/* typeOk is only interesting for 16-bit servers */
	int ftTimeout = 20;			/* in milliseconds, wait up to 20ms for new data/events */
	int selTimeout = 0;			/* in microseconds, for select */
	unsigned int numBlock = 0;	/* running count of data blocks received */
	samples_events_t lastNum = {0,0};	/* number of samples + events handled so far */
	samples_events_t curNum = {0,0};	/* ... currently available */
	fd_set readSet, writeSet;	/* for select on server + child sockets */
	int fdMax = 1;				/* for select, not really used on Windows */	
	struct timeval tv;			/* for select timeout */
	int numClients = 0;			/* current number of clients */
	int newNumClients = 0;      /* number of clients after error checks */
	int opState = 0; 			/* 	0 = waiting for clients or FT header, 
									1 = running, 
									2 = new header received, have all clients receive a stop message
								*/
	
	if (SC==NULL) return NULL;
	
	SC->is_running = 1;
	/* Set typeOk flag to 1 for floats, 0 for int16 
		(in the latter case we need to check the FT header first)
	*/
	typeOk = SC->use16bit ? 0 : 1;
			
	#ifndef PLATFORM_WIN32
	fdMax = SC->server_socket;
	#endif
	
	/* Loop this until errors occur or this flag is set from another thread */
	while (!SC->should_exit) {
		int sel;
		
		/* First, in case we have no header,  we need to read it from the FT buffer	*/
		if (startItem == NULL && opState == 0) {
			startItem = rda_aux_get_hdr_prep_start(SC->ft_buffer, &ftHdr);
			if (startItem == NULL) {
				/* no header yet, wait in select call for clients connecting */
				selTimeout = 100000;
			} else {
				/* yeah, we got it */
				if (SC->verbosity > 4) {
					printf("Picked up FieldTrip header: %i channels @ %.1f Hz, datatype=%i\n",ftHdr.nchans, ftHdr.fsample, ftHdr.data_type); 
				}
				/* don't wait in select call, but inside FT polling */
				selTimeout = 0;
				/* set 'typeOk' flag if we're running a 16 bit server */
				if (SC->use16bit) {
					typeOk = (ftHdr.data_type == DATATYPE_INT16) || (ftHdr.data_type == DATATYPE_UINT16);
				}
				/* start from the numbers of samples + events currently in the buffer */
				lastNum.nsamples = ftHdr.nsamples;
				lastNum.nevents  = ftHdr.nevents;
				/* if we already have clients, change operation state */
				if (numClients > 0) {
					opState = 1;
					/* add start packet to all clients */
					for (i=0;i<numClients;i++) {
						clients[i].item = startItem;
						clients[i].written = 0;
					}
				}
				/* reset block counter */
				numBlock = 0;
				startItem->blockNumber = -1;
			}
		}
		
		/* Second, we deal with the things specific to client sockets */
		
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
			/* But only clients with a pending start/data/stop item	are added to the write set */
			if (clients[i].item != NULL) FD_SET(clients[i].sock, &writeSet);
		}

		/* Check server and client sockets for possible read and write operations */
		tv.tv_sec  = 0; 
		tv.tv_usec = selTimeout;
		sel = select(fdMax + 1, &readSet, &writeSet, NULL, &tv);
		if (sel == -1) {
			perror("rdaserver_thread -- select");
			continue; /* TODO: think about stopping operation instead */
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
				if (SC->verbosity > 0) {
					fprintf(stderr, "rdaserver_thread: opened connection (%i) to client at %s\n",newSock,inet_ntoa(sa.sin_addr));
				}
				
				clients[numClients].sock = newSock;
				clients[numClients].item = startItem;
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
			if (numClients > 0 && startItem != NULL) {
				/* let's go running */
				opState = 1; 
			}
		}
		
		/* Check for sockets that are ready to be written to */
		if (sel>0) {
			sel = rda_aux_check_writing(sel, &writeSet, numClients, clients, SC->verbosity);
		}
		
		/* Check for sockets on which we can read (=> read 0 bytes means closure ) */
		if (sel>0) {
			newNumClients = rda_aux_check_closure(sel, &readSet, numClients, clients, SC->verbosity);
		} else {
			newNumClients = numClients;
		}
		/* Check for clients that lag behind */
		if (newNumClients > 0) {
			newNumClients = rda_aux_check_slow_clients(numBlock - RDA_MAX_LAG, newNumClients, clients, SC->verbosity);
		}
			
		if (newNumClients < numClients) {
			#ifndef PLATFORM_WIN32
			fdMax = SC->server_socket;
			for (i=0;i<newNumClients;i++) {
				if (clients[i].sock > fdMax) fdMax = clients[i].sock;
			}
			#endif
			pthread_mutex_lock(&SC->mutex);
			SC->num_clients = numClients = newNumClients;
			pthread_mutex_unlock(&SC->mutex);
		}
		if (numClients == 0) {
			/* no clients any more - wait */
			opState = 0;
		}

		/* Ok, (client) network stuff is done, let's see if there is more data */
		
		/* If we already have the header, check for new data (samples / events) */
		if (startItem != NULL && !rda_aux_wait_dat(SC->ft_buffer, &lastNum, &curNum, ftTimeout)) {
			int newBlock;
			
			if (SC->blocksize > 0) {
				/* only send out a DATA packet if there are enough new samples */
				newBlock = (curNum.nsamples - lastNum.nsamples) >= SC->blocksize;
				if (newBlock) {
					curNum.nsamples = lastNum.nsamples + SC->blocksize;
				}
			} else {
				/* blocksize = 0: send out a new block if there are new samples or new events */
				newBlock = curNum.nsamples > lastNum.nsamples || curNum.nevents > lastNum.nevents;
			}
			
			if (newBlock && SC->verbosity > 5) {
				printf("New samples/events: %i, %i (Total: %i, %i)\n", 	curNum.nsamples-lastNum.nsamples, 
																		curNum.nevents-lastNum.nevents,
																		curNum.nsamples, curNum.nevents); 
			}			
			
			/* 	In case the new number of samples is smaller than what we had so far,  
				we need to send a STOP packet to all clients and later re-read the header 
				information.
			*/
			if (curNum.nsamples < lastNum.nsamples) {
				rda_buffer_item_t *stopItem;
			
				free(startItem->data);
				free(startItem);
				startItem = NULL;
				
				if (opState == 1) {
					if (SC->verbosity > 4) {
						printf("Sample count in FieldTrip buffer decreased, sending STOP packet to all clients.\n");
					}

					/* if currently running with clients, send STOP packet */
					stopItem = rda_aux_alloc_item(sizeof(rda_msg_hdr_t));
					if (stopItem == NULL) {
						fprintf(stderr, "Out of memory\n");
						break;
					}			
					
					/* First add a STOP message to be sent to the list */
					memcpy(stopItem->data, _rda_guid, sizeof(_rda_guid));
					((rda_msg_hdr_t *) stopItem->data)->nSize = sizeof(rda_msg_hdr_t);
					((rda_msg_hdr_t *) stopItem->data)->nType = RDA_STOP_MSG;
					stopItem->blockNumber = numBlock++;
					stopItem->next = NULL;	/* restart from opState=0 here */
					latestItem = stopItem;
			
					/* Add stop item to data packet list */
					if (firstDataItem == NULL) {
						firstDataItem = stopItem;
					} else {
						rda_buffer_item_t *last = firstDataItem;
						while (last->next != NULL) {
							last = last->next;
						}
						last->next = stopItem;
					}
					/* ... and also add it to clients that are currently waiting */
					for (i=0;i<numClients;i++) {
						if (clients[i].item == NULL) {
							clients[i].item = stopItem;
							stopItem->refCount++;
						}
					}
					/* switch to waiting-for-stop operation mode */
					opState = 2;
				} else {
					if (SC->verbosity > 4) {
						printf("Sample count in FieldTrip buffer decreased, will re-read header.\n");
					}
				}
			}
			
			if (opState != 1) {
				/* if not running, just take note of the updated quantities */
				lastNum = curNum;
			} else {
				/* If the type is right (for 16 bit servers) AND we've got a new block,
					then read this block and start streaming it out
				*/
				if (typeOk && newBlock && opState==1) {
					/* There's new data to stream out */
					rda_buffer_item_t *item;
				
					item = rda_aux_get_samples_and_markers(SC->ft_buffer, &lastNum, &curNum, numBlock, SC->use16bit);
					if (item != NULL) {
						if (firstDataItem == NULL) {
							firstDataItem = item;
						} else {
							rda_buffer_item_t *last = firstDataItem;
							while (last->next != NULL) last = last->next;
							last->next = item;
						}
						for (i=0;i<numClients;i++) {
							if (clients[i].item == NULL) {
								if (SC->verbosity>6) {
									printf("Adding new job for client %i\n", clients[i].sock);
								}
								clients[i].item = item;
								item->refCount++;
							}
						}
						lastNum = curNum;
						latestItem = item;
						/* modify startItem's block number (for new clients),
						   then increase block number
						*/
						startItem->blockNumber = numBlock++;
					} else {
						fprintf(stderr, "Could not get data from FT buffer or allocate memory\n");
					}
				}
			}
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
		while (firstDataItem != NULL && firstDataItem->refCount == 0) { 
			rda_buffer_item_t *next = firstDataItem->next;
			free(firstDataItem->data);
			free(firstDataItem);
			firstDataItem = next;
		}
		
		if (firstDataItem == NULL) {
			latestItem = NULL;
			if (opState == 2) {
				opState = 0;
			}
		}
		
		/* Update the startItem's (=header) next pointer */
		if (startItem != NULL) startItem->next = latestItem;
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
	if (startItem) {
		free(startItem->data);
		free(startItem);
	}
	/* ... and firstDataItem (including following list elements) */
	while (firstDataItem!=NULL) {
		rda_buffer_item_t *next = firstDataItem->next;
		free(firstDataItem->data);
		free(firstDataItem);
		firstDataItem = next;
	}
	SC->is_running = 0;
	return NULL;
}
 
/* see header file for documentation */
rda_server_ctrl_t *rda_start_server(int ft_buffer, int use16bit, int port, int blocksize, int *errval) {
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
	SC->blocksize = (blocksize < 0) ? 0 : blocksize;

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
