/*
 * Copyright (C) 2010, Stefan Klanke
 * Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
 * Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 */

#ifndef __FtBuffer_h
#define __FtBuffer_h

#include <buffer.h>
#include <SimpleStorage.h>

/** Simple wrapper class for FieldTrip buffer request. Not complete yet.
*/

class FtBufferRequest {
	public:
	
	FtBufferRequest() : m_buf() {
		m_def.version = VERSION;
		m_msg.def = &m_def;
	}
	
	void prepGetHeader() {
		m_def.command = GET_HDR;
		m_def.bufsize = 0;
		m_msg.buf = NULL;		
	}
	
	bool prepPutHeader(UINT32_T numChannels, UINT32_T dataType, float fsample) {
		// This is for safety: If a user ignores this function returning false,
		// and just sends this request to the buffer server, we need to make sure
		// no harm is done - so we sent a small invalid packet
		m_def.command = GET_ERR; 
		m_def.bufsize = 0;
		headerdef_t *hd;
		
		if (!m_buf.resize(sizeof(headerdef_t))) return false;
		m_msg.buf = m_buf.data();
		hd = (headerdef_t *) m_buf.data();
		hd->nchans    = numChannels;
		hd->nsamples  = 0;
		hd->nevents   = 0;
		hd->data_type = dataType;
		hd->fsample   = fsample;
		hd->bufsize   = 0;
		m_def.command = PUT_HDR;
		m_def.bufsize = sizeof(headerdef_t);
		return true;
	}	
	
	bool prepPutHeaderAddChunk(UINT32_T chunkType, UINT32_T chunkSize, const void *data) {
		if (m_def.command != PUT_HDR) return false;
	
		UINT32_T oldSize = m_buf.size();
		UINT32_T newSize = oldSize + sizeof(ft_chunkdef_t) + chunkSize;
		
		if (!m_buf.resize(newSize)) return false;
		
		m_def.bufsize = newSize;
		m_msg.buf = m_buf.data();
		
		ft_chunk_t *chunk = (ft_chunk_t *) ((char *) m_msg.buf + oldSize);
		chunk->def.type = chunkType;
		chunk->def.size = chunkSize;
		memcpy(chunk->data, data, chunkSize);
		// increment headerdef->bufsize by total size of chunk
		headerdef_t *hd = (headerdef_t *) m_msg.buf;
		hd->bufsize += chunkSize + sizeof(ft_chunkdef_t);
		return true;
	}		

	bool prepPutData(UINT32_T numChannels, UINT32_T numSamples, UINT32_T dataType, const void *data) {
		// This is for safety: If a user ignores this function returning false,
		// and just sends this request to the buffer server, we need to make sure
		// no harm is done - so we sent a small invalid packet
		m_def.command = GET_ERR; 
		m_def.bufsize = 0;
		
		UINT32_T  wordSize = wordsize_from_type(dataType);
		if (wordSize == 0) return false;
		
		UINT32_T  dataSize = wordSize * numSamples * numChannels;
		UINT32_T  totalSize = dataSize + sizeof(datadef_t);
		datadef_t *dd;
		
		if (!m_buf.resize(totalSize)) return false;
		m_msg.buf = m_buf.data();
		dd = (datadef_t *) m_buf.data();
		dd->nchans = numChannels;
		dd->nsamples = numSamples;
		dd->data_type = dataType;
		dd->bufsize = dataSize;
		// dd+1 points to the next byte after the data_def
		memcpy(dd+1, data, dataSize);
		m_def.command = PUT_DAT;
		m_def.bufsize = totalSize;
		return true;
	}
	
	bool prepPutEvent(INT32_T sample, INT32_T offset, INT32_T duration, const char *type=NULL, const char *value=NULL) {
		// This is for safety: If a user ignores this function returning false,
		// and just sends this request to the buffer server, we need to make sure
		// no harm is done - so we sent a small invalid packet
		m_def.command = GET_ERR; 
		m_def.bufsize = 0;
		
		int len_type  = (type == NULL) ? 0 : strlen(type);
		int len_value = (value == NULL) ? 0 : strlen(value);
		
		UINT32_T  totalSize = sizeof(eventdef_t) + len_type + len_value;
		eventdef_t *ed;
		
		if (!m_buf.resize(totalSize)) return false;
		m_msg.buf = m_buf.data();
		
		ed = (eventdef_t *) m_buf.data();
		ed->type_type = DATATYPE_CHAR;
		ed->type_numel = len_type;
		ed->value_type = DATATYPE_CHAR;
		ed->value_numel = len_value;
		ed->sample = sample;
		ed->offset = offset;
		ed->duration = duration;
		ed->bufsize  = ed->type_numel + ed->value_numel;
		
		if (len_type>0) {
			char *dest = (char *) m_buf.data() + sizeof(eventdef_t);
			memcpy(dest, type, len_type);
			if (len_value>0) {
				memcpy(dest + len_type, value, len_value);
			}
		}
		m_def.command = PUT_EVT;
		m_def.bufsize = totalSize;
		return true;
	}
		
	void prepGetData(UINT32_T begsample, UINT32_T endsample) {
		m_def.command = GET_DAT;
		m_msg.buf = &m_extras.ds;
		m_def.bufsize = sizeof(datasel_t);
		m_extras.ds.begsample = begsample;
		m_extras.ds.endsample = endsample;
	}
		
	void prepWaitData(UINT32_T threshold, UINT32_T milliseconds) {
		m_def.command = WAIT_DAT;
		m_msg.buf = &m_extras.wd;
		m_def.bufsize = sizeof(waitdef_t);
		m_extras.wd.threshold = threshold;
		m_extras.wd.milliseconds = milliseconds;
	}

	const message_t *out() const {
		return &m_msg;
	}
		
	protected:
	
	SimpleStorage m_buf;
	message_t m_msg;
	messagedef_t m_def;
	// we'll always use only one of these at a time
	// so let's define a union to save some (stack) space
	union {	
		waitdef_t wd;
		datasel_t ds;
	} m_extras;
};


/** Simple wrapper class for FieldTrip buffer responses. Not complete yet.
*/
class FtBufferResponse {
	public:
	
	FtBufferResponse() {
		m_response = NULL;
	}
		
	bool checkGetHeader(headerdef_t &hdr, SimpleStorage *bufStore = NULL) const {
		if (m_response == NULL) return false;
		if (m_response->def == NULL) return false;
		if (m_response->def->version != VERSION) return false;
		if (m_response->def->command != GET_OK) return false;
		if (m_response->def->bufsize < sizeof(headerdef_t)) return false;
		if (m_response->buf == NULL) return false;
		
		memcpy(&hdr, m_response->buf, sizeof(headerdef_t));
		if (bufStore != NULL) {
			unsigned int len = m_response->def->bufsize - sizeof(headerdef_t);
			char *src = (char *) m_response->buf + sizeof(headerdef_t);
			if (!bufStore->resize(len)) return false;
			memcpy(bufStore->data(), src, len);
		}		
		return true;
	}
	
	bool checkGetData(datadef_t &datadef, SimpleStorage *datStore = NULL) const {
		if (m_response == NULL) return false;
		if (m_response->def == NULL) return false;
		if (m_response->def->version != VERSION) return false;
		if (m_response->def->command != GET_OK) return false;
		if (m_response->def->bufsize < sizeof(datadef_t)) return false;
		if (m_response->buf == NULL) return false;
		
		memcpy(&datadef, m_response->buf, sizeof(datadef_t));
		if (datStore != NULL) {
			unsigned int len = m_response->def->bufsize - sizeof(datadef_t);
			char *src = (char *) m_response->buf + sizeof(datadef_t);
			if (!datStore->resize(len)) return false;
			memcpy(datStore->data(), src, len);
		}		
		return true;
	}
	
	bool checkWait(unsigned int &nSamples) const {
		if (m_response == NULL) return false;
		if (m_response->def == NULL) return false;
		if (m_response->def->version != VERSION) return false;
		if (m_response->def->command != WAIT_OK) return false;
		if (m_response->def->bufsize != sizeof(UINT32_T)) return false;
		if (m_response->buf == NULL) return false;
		nSamples = *((UINT32_T *)m_response->buf);
		return true;
	}	
	
	message_t **in() {
		if (m_response != NULL) clearResponse();
		return &m_response;
	}
	
	bool checkPut() {
		if (m_response == NULL) return false;
		if (m_response->def == NULL) return false;
		if (m_response->def->version != VERSION) return false;
		return m_response->def->command == PUT_OK;
	}
	
	bool checkFlush() {
		if (m_response == NULL) return false;
		if (m_response->def == NULL) return false;
		if (m_response->def->version != VERSION) return false;
		return m_response->def->command == FLUSH_OK;
	}	
	
	~FtBufferResponse() {
		if (m_response != NULL) clearResponse();
	}

	void clearResponse() {
		if (m_response->buf != NULL) free(m_response->buf);
		if (m_response->def != NULL) free(m_response->def);
		free(m_response);
		m_response = NULL;
	}
	
	message_t *m_response;
};

class FtChunkIterator {
	public:
	
	FtChunkIterator(SimpleStorage& buf) : store(buf) {
		pos = 0;
	}
	
	ft_chunk_t *getNext() {
		if (pos + sizeof(ft_chunkdef_t) > store.size()) return NULL;
		ft_chunk_t *chunk = (ft_chunk_t * ) ((char *) store.data() + pos);
		pos += sizeof(ft_chunkdef_t) + chunk->def.size;
		return chunk;
	}
	
	protected:
	SimpleStorage &store;
	int pos;
};

#endif
