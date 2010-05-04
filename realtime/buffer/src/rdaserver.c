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

#define RDA_FIXED_NAME_LEN   7
#define RDA_MAX_NUM_CLIENTS  32			/* On Windows, this number must be smaller than 64 */

const UINT8_T _rda_guid[16]={
	0x8E,0x45,0x58,0x43,0x96,0xC9,0x86,0x4C,0xAF,0x4A,0x98,0xBB,0xF6,0xC9,0x14,0x50
};

/** Retrieve the header from a FieldTrip buffer and put the definition part into hdr 
	@param ft_buffer  	FieldTrip connection
	@param hdr			Receives header defintion (must not be NULL!)
	@return 	0 on success, -1 on error
*/
int rda_aux_get_hdr(int ft_buffer, headerdef_t *hdr) {
	int r;
	message_t req, *resp = NULL;
	messagedef_t msg_def;
	
	
	req.def = &msg_def;
	req.buf = NULL;
	msg_def.version = VERSION;
	msg_def.command = GET_HDR;
	msg_def.bufsize = 0;
	
	r = clientrequest(ft_buffer, &req, &resp);
	if (r<0) return r;
	
	if (resp == NULL) return -1;
	if (resp->def == NULL) {
		if (resp->buf != NULL) free(resp->buf);
		free(resp);
		return -1;
	}
	if (resp->buf == NULL) {
		free(resp->def);
		free(resp);
		return -1;
	}
	if (resp->def->command == GET_OK) {
		memcpy(hdr, resp->buf, sizeof(headerdef_t));
		r = 0;
	} else {
		r = -1;
	}
	free(resp->buf);
	free(resp->def);
	free(resp);
	return r;
}

/** Prepares "start" RDA packet with fake channel names (for now)
	@param hdr	Fieldtrip header definition
	@param item	Receives the packet (usually "startItem")
	@return 0 on success, -1 if out of memory 
*/
int rda_aux_prep_start(const headerdef_t *hdr, rda_buffer_item_t *item) {
	rda_msg_start_t *R;
	size_t bytesTotal = sizeof(rda_msg_start_t) + hdr->nchans*(sizeof(double) + RDA_FIXED_NAME_LEN*sizeof(char));
	char *str;
	double *dRes;
	int i;
		
	if (bytesTotal > item->sizeAlloc) {
		if (item->data != NULL) free(item->data);
		item->data = malloc(bytesTotal);
		if (item->data == NULL) {
			item->sizeAlloc = 0;
			return -1;
		}
		item->sizeAlloc = bytesTotal;
	}
	item->size = bytesTotal;
	
	R = (rda_msg_start_t *) item->data;
	memcpy(R->hdr.guid, _rda_guid, sizeof(_rda_guid));
	R->hdr.nSize = bytesTotal;
	R->hdr.nType = RDA_START_MSG;
	R->nChannels = hdr->nchans;
	R->dSamplingInterval = 1.0/hdr->fsample;
	
	/* R+1 points to first byte after header info */
	dRes = (double *) ((void *)(R+1)); 
	/* Fill with resolution = 1.0 -- we have nothing better */
	for (i=0;i<hdr->nchans;i++) dRes[i]=1.0;
	
	/* Let 'str' point to first byte after the resolution values */
	str = (char *) ((void *)(dRes + hdr->nchans));
	for (i=0;i<hdr->nchans;i++) {
		snprintf(str, RDA_FIXED_NAME_LEN, "%06d", i+1);
		str+=RDA_FIXED_NAME_LEN;
	}
	/* done */
	return 0;
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


/** Get data from FT buffer and convert it into a single precision RDA packet
	@return  0 on success, -1 on connection errors, -2 on out of memory
*/
int rda_aux_get_float_data(int ft_buffer, unsigned int start, unsigned int end, rda_buffer_item_t *item) {
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
		rda_msg_data_t *R;
		datadef_t *ddef = (datadef_t *) resp->buf;
		unsigned int numTotal = ddef->nchans * ddef->nsamples;
		size_t bytesTotal = numTotal*sizeof(float) + sizeof(rda_msg_data_t);
		
		if (bytesTotal > item->sizeAlloc) {
			if (item->data != NULL) free(item->data);
			item->data = malloc(bytesTotal);
			if (item->data == NULL) {
				item->sizeAlloc = 0;
				r = -2;
				goto cleanup;
			}
			item->sizeAlloc = bytesTotal;
		}
		item->size = bytesTotal;
		/* Ok, we've got the data and the space to put it into
		   but first we set up the header */
		R = (rda_msg_data_t *) item->data;
		memcpy(R->hdr.guid, _rda_guid, sizeof(_rda_guid));
		R->hdr.nSize = bytesTotal;
		R->hdr.nType = RDA_FLOAT_MSG;
		R->nBlock = 0; /* will be filled later */
		R->nPoints = ddef->nsamples;
		R->nMarkers = 0; /* maybe later */
		
		/* convert the samples, note that (R+1) points to the first byte after the RDA data header
		   and (ddef+1) points to the first byte after the FieldTrip data header */
		rda_aux_convert_to_float(numTotal, (void *)(R+1), ddef->data_type, (void *)(ddef+1));
		
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
int rda_aux_get_int16_data(int ft_buffer, unsigned int start, unsigned int end, rda_buffer_item_t *item) {
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
		rda_msg_data_t *R;
		datadef_t *ddef = (datadef_t *) resp->buf;
		unsigned int numTotal = ddef->nchans * ddef->nsamples;
		size_t bytesTotal = numTotal*sizeof(INT16_T) + sizeof(rda_msg_data_t);
		
		if ((ddef->data_type != DATATYPE_INT16) && (ddef->data_type != DATATYPE_UINT16)) {
			r = -1;
			goto cleanup;
		}
		
		if (bytesTotal > item->sizeAlloc) {
			if (item->data != NULL) free(item->data);
			item->data = malloc(bytesTotal);
			if (item->data == NULL) {
				item->sizeAlloc = 0;
				r = -2;
				goto cleanup;
			}
			item->sizeAlloc = bytesTotal;
		}
		item->size = bytesTotal;
		/* Ok, we've got the data and the space to put it into
		   but first we set up the header */
		R = (rda_msg_data_t *) item->data;
		memcpy(R->hdr.guid, _rda_guid, sizeof(_rda_guid));
		R->hdr.nSize = bytesTotal;
		R->hdr.nType = RDA_INT_MSG;
		R->nBlock = 0; /* will be filled later */
		R->nPoints = ddef->nsamples;
		R->nMarkers = 0; /* maybe later */
			
		/* copy the samples, note that (R+1) points to the first byte after the RDA data header
		   and (ddef+1) points to the first byte after the FieldTrip data header */
		memcpy((void *)(R+1), (void *)(ddef+1), numTotal*sizeof(INT16_T));
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
	rda_buffer_item_t *startItem    = NULL; 
	rda_buffer_item_t *firstDataItem = NULL;
	
	headerdef_t ftHdr;
	int i,typeOk;
	
	unsigned int lastNumSamples = 0;
	unsigned int lastNumEvents  = 0;
	
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
			/* If we can't get the header, go back to the start of the loop */
			if (rda_aux_get_hdr(SC->ft_buffer, &ftHdr)) continue;
			/* If this is a 16-bit server, we can only handle 16-bit FieldTrip buffers */
			if (SC->use16bit) {
				typeOk = (ftHdr.data_type == DATATYPE_INT16) || (ftHdr.data_type == DATATYPE_UINT16);
				if (!typeOk) continue;	/* TODO: think about a possibly better way of handling this */
			}
			/* Is this the first time we read a header (and succesfully convert it to an RDA packet)? */
			if (startItem == NULL) {
				/* Yes, prepare startItem */
				startItem = (rda_buffer_item_t *) malloc(sizeof(rda_buffer_item_t));
				if (startItem == NULL) { 
					fprintf(stderr,"rdaserver_thread: Out of memory\n");
					break;
				}
				startItem->sizeAlloc = 0;
				startItem->data = NULL;
				startItem->next = NULL;
				/* startItem does not need a refCount: it's always kept */
				
				if (rda_aux_prep_start(&ftHdr, startItem)) {
					fprintf(stderr,"rdaserver_thread: Out of memory\n");
					break;
				}
				/* Now the start item should have proper size/sizeAlloc/data fields */
				
				if (SC->verbosity > 4) {
					printf("Picked up FT header for first time!\n"); 
				}
				
				/* Add the start item to every client */
				for (i=0;i<numClients;i++) {
					clients[i].item = startItem;
					clients[i].written = 0;
				}
				lastNumSamples = ftHdr.nsamples;
				lastNumEvents  = ftHdr.nevents;
			} else if (typeOk && ftHdr.nsamples > lastNumSamples) {
				/* No, but there's new data to stream out */
				rda_buffer_item_t *item, *last;
				int newItem = 0;
				
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
					printf("Trying to get data [%i ; %i]\n", lastNumSamples, ftHdr.nsamples); 
				}
				if (SC->use16bit) {
					i = rda_aux_get_int16_data(SC->ft_buffer, lastNumSamples, ftHdr.nsamples-1, item);
				} else {
					i = rda_aux_get_float_data(SC->ft_buffer, lastNumSamples, ftHdr.nsamples-1, item);
				}
				if (i!=0) {
					fprintf(stderr, "Error in aux_get_float_data!\n");
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
				lastNumSamples = ftHdr.nsamples;
				lastNumEvents  = ftHdr.nevents;
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
			startItem->next = firstDataItem;
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
	if (startItem) {
		if (startItem->data) free(startItem->data);
		free(startItem);
	}
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
	
#ifdef PLATFORM_WIN32
	WSADATA wsa;
 	if(WSAStartup(MAKEWORD(1, 1), &wsa))
	{
		fprintf(stderr, "tcpserver: cannot start sockets\n");
		goto cleanup;
	}
#endif

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
