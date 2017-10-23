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

struct FtDataType {
	template<typename T> static UINT32_T getType(T dummy) {
		return DATATYPE_UNKNOWN;
	}

	// "specialised versions" of the above template function, returning the
	// right GDF_Type for the given input argument
	static UINT32_T getType(INT8_T dummy)   { return DATATYPE_INT8; }
	static UINT32_T getType(UINT8_T dummy)  { return DATATYPE_UINT8; }
	static UINT32_T getType(INT16_T dummy)  { return DATATYPE_INT16; }
	static UINT32_T getType(UINT16_T dummy) { return DATATYPE_UINT16; }
	static UINT32_T getType(INT32_T dummy)  { return DATATYPE_INT32; }
	static UINT32_T getType(UINT32_T dummy) { return DATATYPE_UINT32; }
	static UINT32_T getType(INT64_T dummy)  { return DATATYPE_INT64; }
	static UINT32_T getType(UINT64_T dummy) { return DATATYPE_UINT64; }
	static UINT32_T getType(float dummy)    { return DATATYPE_FLOAT32; }
	static UINT32_T getType(double dummy)   { return DATATYPE_FLOAT64; }
};


/** Simple wrapper class for FieldTrip buffer requests and responses.
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

		unsigned int wordSize = wordsize_from_type(dataType);
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

	bool prepPutEvent(INT32_T sample, INT32_T offset, INT32_T duration, const char *type=NULL, INT32_T value = 0) {
		// This is for safety: If a user ignores this function returning false,
		// and just sends this request to the buffer server, we need to make sure
		// no harm is done - so we sent a small invalid packet
		m_def.command = GET_ERR;
		m_def.bufsize = 0;

		int len_type  = (type == NULL) ? 0 : strlen(type);

		UINT32_T  totalSize = sizeof(eventdef_t) + len_type + sizeof(INT32_T);
		eventdef_t *ed;

		if (!m_buf.resize(totalSize)) return false;
		m_msg.buf = m_buf.data();

		ed = (eventdef_t *) m_buf.data();
		ed->type_type = DATATYPE_CHAR;
		ed->type_numel = len_type;
		ed->value_type = DATATYPE_INT32;
		ed->value_numel = 1;
		ed->sample = sample;
		ed->offset = offset;
		ed->duration = duration;
		ed->bufsize  = ed->type_numel + sizeof(DATATYPE_INT32);

		char *dest = (char *) m_buf.data() + sizeof(eventdef_t);

		if (len_type>0) {
			memcpy(dest, type, len_type);
			dest+=len_type;
		}
		memcpy(dest, &value, sizeof(INT32_T));

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

	void prepGetEvents(UINT32_T begevent, UINT32_T endevent) {
		m_def.command = GET_EVT;
		m_msg.buf = &m_extras.es;
		m_def.bufsize = sizeof(eventsel_t);
		m_extras.es.begevent = begevent;
		m_extras.es.endevent = endevent;
	}

	void prepWaitData(UINT32_T nSamples, UINT32_T nEvents, UINT32_T milliseconds) {
		m_def.command = WAIT_DAT;
		m_msg.buf = &m_extras.wd;
		m_def.bufsize = sizeof(waitdef_t);
		m_extras.wd.threshold.nsamples = nSamples;
		m_extras.wd.threshold.nevents = nEvents;
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
		eventsel_t es;
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

	int checkGetEvents(SimpleStorage *evtStore = NULL) const {
		if (m_response == NULL) return false;
		if (m_response->def == NULL) return false;
		if (m_response->def->version != VERSION) return false;
		if (m_response->def->command != GET_OK) return false;

		int check = check_event_array(m_response->def->bufsize, m_response->buf);
		if (check < 0) return check;

		if (evtStore != NULL && evtStore->resize(m_response->def->bufsize)) {
			memcpy(evtStore->data(), m_response->buf, m_response->def->bufsize);
		}
		return check;
	}

	bool checkWait(unsigned int &nSamples, unsigned int &nEvents) const {
		if (m_response == NULL) return false;
		if (m_response->def == NULL) return false;
		if (m_response->def->version != VERSION) return false;
		if (m_response->def->command != WAIT_OK) return false;
		if (m_response->def->bufsize != sizeof(samples_events_t)) return false;
		if (m_response->buf == NULL) return false;
		samples_events_t *nse = (samples_events_t *) m_response->buf;
		nSamples = nse->nsamples;
		nEvents  = nse->nevents;
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


/*
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
*/


/** Small class for generating PUT_EVT requests.
	In comparison to the prepPutEvent method in FtBufferRequest,
	this class allows multiple events in one request.
	The class is designed to be used in realtime acquisition loops,
	and	will keep and recycle allocated memory until explicitly disposed().
*/
class FtEventList {
	public:

	/* Creates an empty event list */
	FtEventList() {
		numEvs=sizeAlloc=0;
		buf=NULL;
		request.def = &reqdef;
		reqdef.command = PUT_EVT;
		reqdef.version = VERSION;
		reqdef.bufsize = 0;
	}

	~FtEventList() {
		if (buf!=NULL) free(buf);
	}

	/* Clears all events from the list, but does not free the memory */
	void clear() {
		reqdef.bufsize = 0;
		numEvs = 0;
	}

	/* Clears all events from the list and frees used memory */
	void dispose() {
		if (buf!=NULL) free(buf);
		buf = 0;
		reqdef.bufsize = 0;
		numEvs = 0;
	}

	/* Add an event of general type to the list, will silently ignore invalid types */
	void add(int sample, UINT32_T type_type, UINT32_T type_numel, const void *type, UINT32_T value_type, UINT32_T value_numel, const void *value) {
		eventdef_t *ne;
		unsigned int typeSize  = type_numel*wordsize_from_type(type_type);
		unsigned int valueSize = value_numel*wordsize_from_type(value_type);
		unsigned int newSize   = reqdef.bufsize + sizeof(eventdef_t) + typeSize + valueSize;

		if (typeSize == 0 || valueSize == 0) return;

		if (sizeAlloc < newSize) {
			char *newBuf;
			if (sizeAlloc == 0) {
				newBuf = (char *) malloc(newSize);
			} else {
				newBuf = (char *) realloc(buf, newSize);
			}
			if (newBuf == NULL) {
				fprintf(stderr, "Warning: out of memory in re-allocating event list.\n");
				return;
			}
			sizeAlloc = newSize;
			buf = newBuf;
		}

		ne = (eventdef_t *) (buf + reqdef.bufsize);
		ne->type_type   = type_type;
		ne->type_numel  = type_numel;
		ne->value_type  = value_type;
		ne->value_numel = value_numel;
		ne->sample = sample;
		ne->offset = 0;
		ne->duration = 0;
		ne->bufsize = typeSize + valueSize;
		memcpy(buf + reqdef.bufsize + sizeof(eventdef_t), type, typeSize);
		memcpy(buf + reqdef.bufsize + sizeof(eventdef_t) + typeSize, value, valueSize);
		reqdef.bufsize = newSize;
		numEvs++;
	}

	/* Add an event where both type and value are C strings */
	void add(int sample, const char *type, const char *value) {
		add(sample, DATATYPE_CHAR, strlen(type), type, DATATYPE_CHAR, strlen(value), value);
	}
	/* Add an event with a string type, and signed 16-bit value */
	void add(int sample, const char *type, INT16_T value) {
		add(sample, DATATYPE_CHAR, strlen(type), type, DATATYPE_INT16, 1, &value);
	}
	/* Add an event with a string type, and signed 32-bit value */
	void add(int sample, const char *type, INT32_T value) {
		add(sample, DATATYPE_CHAR, strlen(type), type, DATATYPE_INT32, 1, &value);
	}
	/* Add an event with a string type, and double precision value */
	void add(int sample, const char *type, double value) {
		add(sample, DATATYPE_CHAR, strlen(type), type, DATATYPE_FLOAT64, 1, &value);
	}
	/* Add an event with a string type, and single precision value */
	void add(int sample, const char *type, float value) {
		add(sample, DATATYPE_CHAR, strlen(type), type, DATATYPE_FLOAT32, 1, &value);
	}

	/* Returns the request corresponding to PUT_EVT of the current list */
	const message_t *asRequest() {
		request.buf = (reqdef.bufsize == 0) ? NULL : buf;
		return &request;
	}

	/* Returns the number of events currently in the list */
	int count() const {
		return numEvs;
	}

	/* Transforms sample indices of events by first adding "offset", then dividing by divisor */
	void transform(int offset, int divisor) {
		unsigned int pos = 0;
		for (unsigned int i=0;i<numEvs;i++) {
			eventdef_t *ne = (eventdef_t *) (buf + pos);

			ne->sample = (ne->sample + offset) / divisor;
			pos += sizeof(eventdef_t) + ne->bufsize;
		}
	}

	protected:

	messagedef_t reqdef;
	message_t request;
	char *buf;
	unsigned int numEvs, sizeAlloc;
};

/** Small class for generating PUT_DAT requests.
	The class is designed to be used in realtime acquisition loops,
	and	will keep and recycle allocated memory until explicitly disposed().
*/
class FtSampleBlock {
	public:
	/* Create a sample block class for the given data type */
	FtSampleBlock(UINT32_T datatype) {
		sizeAlloc=0;
		request.def = &reqdef;
		reqdef.command = PUT_DAT;
		reqdef.version = VERSION;
		reqdef.bufsize = 0;
		this->datatype = datatype;
		wordsize = wordsize_from_type(datatype);
		ddef = NULL;
	}

	~FtSampleBlock() {
		if (ddef!=NULL) free(ddef);
	}

	/* force free'ing memory */
	void dispose() {
		if (ddef!=NULL) free(ddef);
		ddef = 0;
		reqdef.bufsize = 0;
	}

	/* Tries to allocate memory for a block for numChannels x numSamples for the
		data type specified during construction, or recycles already existing
		memory. Returns a pointer to the first byte of the sample block, or NULL
		in case memory could not be allocated.
	*/
	void *getMatrix(int numChannels, int numSamples) {
		unsigned int newSize = sizeof(datadef_t) + numChannels*numSamples*wordsize;

		if (newSize > sizeAlloc) {
			if (ddef!=NULL) free(ddef);
			ddef = (datadef_t *) malloc(newSize);
			if (ddef == NULL) {
				reqdef.bufsize = 0;
				sizeAlloc = 0;
				return NULL;
			}
			sizeAlloc = newSize;
		}
		ddef->nchans    = numChannels;
		ddef->nsamples  = numSamples;
		ddef->data_type = datatype;
		ddef->bufsize   = numChannels * numSamples * wordsize;
		reqdef.bufsize  = ddef->bufsize + sizeof(datadef_t);
		return (void *) (ddef+1); // first byte after datadef_t
	}

	const message_t *asRequest() {
		request.buf = (reqdef.bufsize == 0) ? NULL : ddef;
		return &request;
	}

	protected:

	UINT32_T datatype;
	unsigned int wordsize;
	datadef_t *ddef;
	unsigned int sizeAlloc;
	messagedef_t reqdef;
	message_t request;
};


class FtConnection {
	public:

	FtConnection(int retry=0) {
		this->retry  = (retry < 0) ? 0 : retry;
		this->sock   = -1;
		this->type   = -1;
		++numConnections;
	}

	~FtConnection();

	void setRetry(int retry) 	{ this->retry = (retry < 0) ? 0 : retry; }
	int getSocket() const 		{ return sock; }
	int getType() const 		{ return type; }
	bool isOpen() const 		{ return (sock > -1); }

	bool connect(const char *address);
	bool connectDirect() {
		disconnect();
		sock = type = 0;
		return true;
	}

	bool connectTcp(const char *hostname, int port);
	bool connectUnix(const char *pathname);
	void disconnect() {
		if (sock > 0) closesocket(sock);
		sock = -1;
		type = -1;
	}

	protected:

	#ifdef WIN32
	static WSADATA wsa;
	#endif
	static int numConnections;
	int sock, type, retry;
};

#endif
