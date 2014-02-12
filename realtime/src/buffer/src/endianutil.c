/* 
 * Copyright (C) 2010, Stefan Klanke
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 */

#include "buffer.h"

/* TODO: see if these can be optimized using compiler intrinsics etc. */
void ft_swap16(unsigned int numel, void *data) {
	unsigned int n;
	char *d = (char *) data;
	for (n=0;n<numel;n++) {
		char t = d[0];
		d[0] = d[1];
		d[1] = t;
		d+=2;
	}
}

void ft_swap32(unsigned int numel, void *data) {
	unsigned int n;
	char *d = (char *) data;
	for (n=0;n<numel;n++) {
		char t0 = d[0];
		char t1 = d[1];
		d[0] = d[3];
		d[1] = d[2];
		d[2] = t1;
		d[3] = t0;
		d+=4;
	}
}
		
void ft_swap64(unsigned int numel, void *data) {
	unsigned int n;
	char *d = (char *) data;
	for (n=0;n<numel;n++) {
		char t0 = d[0];
		char t1 = d[1];
		char t2 = d[2];
		char t3 = d[3];
		d[0] = d[7];
		d[1] = d[6];
		d[2] = d[5];
		d[3] = d[4];
		d[4] = t3;
		d[5] = t2;
		d[6] = t1;
		d[7] = t0;
		d+=8;
	}
}

void ft_swap_data(UINT32_T numel, UINT32_T datatype, void *data) {
	switch(datatype) {
		case DATATYPE_CHAR:
		case DATATYPE_UINT8:
		case DATATYPE_INT8:
			return;
		case DATATYPE_UINT16:
		case DATATYPE_INT16:
			ft_swap16(numel, data);
			return;
		case DATATYPE_UINT32:
		case DATATYPE_INT32:
		case DATATYPE_FLOAT32:
			ft_swap32(numel, data);
			return;
		case DATATYPE_UINT64:
		case DATATYPE_INT64:
		case DATATYPE_FLOAT64:
			ft_swap64(numel, data);
			return;
	}	
}

int ft_swap_chunks_to_native(UINT32_T size, UINT32_T nchans, void *buf) {
	UINT32_T offset = 0;
	while (offset + sizeof(ft_chunkdef_t) <= size) {
		ft_chunk_t *chunk = (ft_chunk_t *) ((char *) buf + offset);
		
		ft_swap32(2, &(chunk->def));
		
		offset += sizeof(ft_chunkdef_t) + chunk->def.size;
		/* chunk definition fault (=too big) ? */
		if (offset > size) return -1;
		
		switch(chunk->def.type) {
			case FT_CHUNK_RESOLUTIONS:
				if (chunk->def.size >= nchans*sizeof(FLOAT64_T)) {
					ft_swap64(nchans, chunk->data);
				}
				break;
			/* Add other cases here as needed */
		}		
		
		offset += sizeof(ft_chunkdef_t) + chunk->def.size;
	}
	return 0;
}

/* returns 0 on success, -1 on error */
int ft_swap_events_to_native(UINT32_T size, void *buf) {
	UINT32_T offset = 0;

	while (offset + sizeof(eventdef_t) <= size) {
		unsigned int wst, wsv;
		
		eventdef_t *edef = (eventdef_t *) ((char *) buf + offset);
		ft_swap32(8, edef); /* all fields are 32-bit */
		
		/* Increase offset to beginning of next event */
		offset += sizeof(eventdef_t) + edef->bufsize;
		if (offset > size) return -1;	/* this event is too big for "buf" */
		
		wst = wordsize_from_type(edef->type_type);
		wsv = wordsize_from_type(edef->value_type);
		
		/* check if type and value fit into this event's local buffer */
		if (wst*edef->type_numel + wsv*edef->value_numel > edef->bufsize) return -1;
		
		ft_swap_data(edef->type_numel, edef->type_type, (char *) buf + offset);
		ft_swap_data(edef->value_numel, edef->value_type, (char *) buf + offset + wst*edef->type_numel);
	}
	return 0;
}


/* returns 0 on success, -1 on error */
int ft_swap_buf_to_native(UINT16_T command, UINT32_T bufsize, void *buf) {
	datadef_t *ddef;
	
	switch(command) {
		case GET_HDR:
			/* This should not have a buf attached */
			return 0;
		case GET_DAT:
			/* buf contains a datsel_t = 2x UINT32_T */
			if (bufsize == 8) ft_swap32(2, buf);
			return 0;
		case GET_EVT:
			/* buf contains a datsel_t = 2x UINT32_T */
			if (bufsize == 8) ft_swap32(2, buf);
			return 0;
		case WAIT_DAT:
			/* buf contains a waitdef_t = 3x UINT32_T */
			ft_swap32(3, buf);
			return 0;
		case PUT_DAT:
			/* buf contains a datadef_t and after that the data */
			ddef = (datadef_t *) buf;
			ft_swap32(4, ddef);	/* this is for datadef_t */
			ft_swap_data(ddef->nchans*ddef->nsamples, ddef->data_type, (char *)ddef + sizeof(datadef_t)); /* ddef+1 points to first data byte */
			return 0;
		case PUT_HDR:
			/* buf contains a headerdef_t and optionally chunks */
			ft_swap32(6, buf);	/* all fields are 32-bit values */
			return ft_swap_chunks_to_native(bufsize - sizeof(headerdef_t), ((headerdef_t *) buf)->nchans, (char *) buf + sizeof(headerdef_t));
		case PUT_EVT:
			/* buf contains multiple eventdef_t and buf's */
			return ft_swap_events_to_native(bufsize, buf);
	}
	return -1;
}



int ft_swap_chunks_from_native(UINT32_T size, UINT32_T nchans, void *buf) {
	UINT32_T offset = 0;
	while (offset + sizeof(ft_chunkdef_t) <= size) {
		ft_chunk_t *chunk = (ft_chunk_t *) ((char *) buf + offset);
		offset += sizeof(ft_chunkdef_t) + chunk->def.size;
		
		/* chunk definition fault (=too big) ? */
		if (offset > size) return -1;
		
		switch(chunk->def.type) {
			case FT_CHUNK_RESOLUTIONS:
				if (chunk->def.size >= nchans*sizeof(FLOAT64_T)) {
					ft_swap64(nchans, chunk->data);
				}
				break;
			/* Add other cases here as needed */
		}
		
		ft_swap32(2, &(chunk->def));
	}
	return 0;
}

int ft_swap_events_from_native(UINT32_T size, void *buf) {
	UINT32_T offset = 0;
	
	while (offset + sizeof(eventdef_t) <= size) {
		unsigned int wst;
		
		eventdef_t *edef = (eventdef_t *) ((char *) buf + offset);
		offset += sizeof(eventdef_t) + edef->bufsize;
		
		wst = wordsize_from_type(edef->type_type);
		
		ft_swap_data(edef->type_numel, edef->type_type, (char *) buf + offset);
		ft_swap_data(edef->value_numel, edef->value_type, (char *) buf + offset + wst*edef->type_numel);
		ft_swap32(8, edef); /* all fields are 32-bit */
	}
	return 0;
}


int ft_swap_from_native(UINT16_T orgCommand, message_t *msg) {
	datadef_t *ddef;
	UINT32_T nchans;
	UINT32_T bufsize = msg->def->bufsize;
	
	ft_swap16(1, &msg->def->version);
	ft_swap16(1, &msg->def->command);
	ft_swap32(1, &msg->def->bufsize);
	
	if (bufsize == 0) return 0;
	
	switch(orgCommand) {
		case GET_HDR:
			nchans = ((headerdef_t *) msg->buf)->nchans;
			ft_swap32(6, msg->buf);	/* all fields are 32-bit values */
			return ft_swap_chunks_from_native(bufsize - sizeof(headerdef_t), nchans, (char *) msg->buf + sizeof(headerdef_t));
		case GET_DAT:
			ddef = (datadef_t *) msg->buf;
			ft_swap_data(ddef->nchans*ddef->nsamples, ddef->data_type, (char *)ddef + sizeof(datadef_t)); /* ddef+1 points to first data byte */
			ft_swap32(4, ddef); /* all fields are 32-bit */
			return 0;
		case GET_EVT:
			return ft_swap_events_from_native(bufsize, msg->buf);
		case WAIT_DAT:
			ft_swap32(2, msg->buf);	/* nsamples + nevents = 32bit */
			return 0;
	}
	return -1;
}
