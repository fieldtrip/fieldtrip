#include <BioSemiClient.h>
#include <stdio.h>
#include <signal.h>
#include <FtBuffer.h>

volatile bool keepRunning = true;

class EventChain {
	public:
	
	EventChain(const char *name=NULL) {
		numEvs=sizeAlloc=0;
		buf=NULL;
		request.def = &reqdef;
		reqdef.command = PUT_EVT;
		reqdef.version = VERSION;
		reqdef.bufsize = 0;
		
		if (name==NULL) {
			strcpy(this->name, "TRIGGER");
			lenName = 7;
		} else {
			lenName = strlen(name);
			if (lenName > sizeof(this->name)) lenName = sizeof(this->name);
			memcpy(this->name, name, lenName);
		}
	}
	
	~EventChain() {
		if (buf!=NULL) free(buf);
	}
	
	void clear() {
		reqdef.bufsize = 0;
		numEvs = 0;
	}
	
	void dispose() {
		if (buf!=NULL) free(buf);
		buf = 0;
		reqdef.bufsize = 0;
		numEvs = 0;
	}
	
	void addTrigger(int sample, INT32_T value) {
		eventdef_t *ne;
		unsigned int newSize = reqdef.bufsize + sizeof(eventdef_t) + lenName + sizeof(INT32_T); // 7 for "TRIGGER" plus 4 for the INT32 value
		if (sizeAlloc < newSize) {
			char *newBuf;
			if (sizeAlloc == 0) {
				newBuf = (char *) malloc(newSize);
			} else {
				newBuf = (char *) realloc(buf, newSize);
			}
			if (newBuf == NULL) {
				fprintf(stderr, "Warning: out of memory in re-allocating event chain.\n");
				return;
			}
			sizeAlloc = newSize;
			buf = newBuf;
		}
		
		ne = (eventdef_t *) (buf + reqdef.bufsize);
		ne->type_type  = DATATYPE_CHAR;
		ne->type_numel = lenName;
		ne->value_type = DATATYPE_INT32;
		ne->value_numel = 1;
		ne->sample = sample;
		ne->offset = 0;
		ne->duration = 0;
		ne->bufsize = lenName + sizeof(INT32_T);
		memcpy(buf + reqdef.bufsize + sizeof(eventdef_t), name, lenName);
		memcpy(buf + reqdef.bufsize + sizeof(eventdef_t) + lenName, &value, sizeof(int));
		reqdef.bufsize = newSize;
		numEvs++;
	}
	
	const message_t *asRequest() { 
		request.buf = (reqdef.bufsize == 0) ? NULL : buf;
		return &request; 
	}	
	
	int count() const {
		return numEvs;
	}
	
	protected:
	
	messagedef_t reqdef;
	message_t request;
	char *buf;
	char name[32];
	unsigned int lenName;
	unsigned int numEvs, sizeAlloc;
};


class SampleBlock {
	public:
	SampleBlock() {
		sizeAlloc=0;
		request.def = &reqdef;
		reqdef.command = PUT_DAT;
		reqdef.version = VERSION;
		reqdef.bufsize = 0;
		ddef = NULL;
	}
	
	~SampleBlock() {
		if (ddef!=NULL) free(ddef);
	}
	
	void dispose() {
		if (ddef!=NULL) free(ddef);
		ddef = 0;
		reqdef.bufsize = 0;
	}
	
	int *getMatrix(int numChannels, int numSamples) {
		unsigned int newSize = sizeof(datadef_t) + numChannels*numSamples*sizeof(INT32_T);
		
		if (newSize > sizeAlloc) {
			if (ddef!=NULL) free(ddef);
			ddef = (datadef_t *) malloc(newSize);
			if (ddef == NULL) {
				reqdef.bufsize = 0;
				return NULL;
			}
		}
		ddef->nchans    = numChannels;
		ddef->nsamples  = numSamples;
		ddef->data_type = DATATYPE_INT32;
		ddef->bufsize   = numChannels * numSamples * sizeof(INT32_T);
		return (int *) (ddef+1); // first bytes after datadef_t
	}
	
	const message_t *asRequest() { 
		request.buf = (reqdef.bufsize == 0) ? NULL : ddef;
		return &request; 
	}		
	
	protected:
	
	datadef_t *ddef;
	unsigned int sizeAlloc;
	messagedef_t reqdef;
	message_t request;
};

void abortHandler(int sig) {
	printf("Ctrl-C pressed -- stopping...\n");
	keepRunning = false;
}

int main(int argc, char *argv[]) {
	int numEEG, numAIB, numTrigger;
	int port, ftSocket;
	int sampleCounter = 0;
	char hostname[256];
	BioSemiClient BS;
	FtBufferResponse resp;
	EventChain  eventChain("TRIGGER");
	SampleBlock sampleBlock;
	bool triggerState[16];
	
	for (int i=0;i<16;i++) triggerState[i] = false;
	
	if (argc<4) {
		printf("Usage: biosemi2ft <#EEG-channels> <#AIB-channels> <#trigger-channels> [hostname=localhost [port=1972]]\n");
		return 0;
	}
	
	numEEG = atoi(argv[1]);
	if (numEEG < 0) {
		fprintf(stderr, "Number of EEG channels (first parameter) must be >=0\n");
		return 1;
	}
	
	numAIB = atoi(argv[2]);
	if (numAIB < 0) {
		fprintf(stderr, "Number of AIB channels (second parameter) must be >=0\n");
		return 1;
	}
	
	numTrigger = atoi(argv[3]);
	if (numTrigger < 0 || numTrigger>16) {
		fprintf(stderr, "Number of trigger channels (third parameter) must be between 0 and 16\n");
		return 1;
	}
	
	if (argc>4) {
		strncpy(hostname, argv[4], sizeof(hostname));
	} else {
		strcpy(hostname, "localhost");
	}
	
	if (argc>5) {
		port = atoi(argv[5]);
	} else {
		port = 1972;
	}
	
	if (numEEG + numAIB + numTrigger > 0) {
		ftSocket = open_connection(hostname, port);
		if (ftSocket < 0) {
			fprintf(stderr, "Could not connect to FieldTrip buffer on %s:%i\n", hostname, port);
			return 1;
		}
	} else {
		// a -1 here means do not write anything, only for inspection
		ftSocket = -1;
	}
	
	if (!BS.openDevice()) return 1;
	
	if (BS.getNumChannels() < numEEG) numEEG = BS.getNumChannels();
	if (BS.getNumChanAIB() < numAIB)  numAIB = BS.getNumChanAIB();
	
	printf("Speed mode.........: %i\n", BS.getSpeedMode());
	printf("Number of channels.: %i, selected first %i\n", BS.getNumChannels(), numEEG);
	printf("Number of AIB chns.: %i, selected first %i\n", BS.getNumChanAIB(), numAIB);
	printf("Sampling frequency.: %i Hz\n", BS.getSamplingFreq());
	
	if (ftSocket!=-1) {
		FtBufferRequest req;

		req.prepPutHeader(numEEG + numAIB, DATATYPE_INT32, BS.getSamplingFreq());
		
		int err = clientrequest(ftSocket, req.out(), resp.in());
		if (err || !resp.checkPut()) {
			fprintf(stderr, "Could not write header to FieldTrip buffer\n.");
			BS.closeDevice();
			close_connection(ftSocket);
		}
	}
		
	/* register CTRL-C handler */
	signal(SIGINT, abortHandler);
	printf("Starting to listen - press CTRL-C to quit\n");
	
	while (keepRunning) {
		BioSemiBlock block;
		bool newBlock;
		int numChan = numEEG + numAIB;
		
		newBlock = BS.checkNewBlock(block);
		
		if (!newBlock) {
			BS.msleep(1);
			continue;
		}
		
		printf("Ptr = %6i,  samples = %4i,   synced = %4i\n", block.startPtr, block.numSamples, block.numInSync);
		
		if (block.numSamples != block.numInSync) continue; // replace by break later
		
		if (numChan > 0) {
			int *dest = sampleBlock.getMatrix(numEEG + numAIB, block.numSamples);
			if (dest==NULL) {
				printf("Out of memory\n");
				break;
			}
			for (int j=0;j<block.numSamples;j++) {
				int ptr = block.startPtr + 8 + j*block.stride;
				for (int i=0;i<numEEG;i++) {
					dest[i + j*numChan] = BS.getValue(ptr + i*4);
				}
				for (int i=0;i<numAIB;i++) {
					dest[i+numEEG+j*numChan] = BS.getValue(ptr + (i+BS.getNumChannels())*4);
				}
			}
			int err = clientrequest(ftSocket, sampleBlock.asRequest(), resp.in());
			if (err || !resp.checkPut()) {
				fprintf(stderr, "Could not write samples to FieldTrip buffer\n.");
			}
			
		}
		if (numTrigger > 0) {
			eventChain.clear();
			for (int j=0;j<block.numSamples;j++) {
				int value = BS.getValue(block.startPtr + 4 + j*block.stride);
				int bv = 0x100; // trig channel 1
				for (int i=0;i<numTrigger;i++) {
					if (value & bv) {
						if (!triggerState[i]) {
							eventChain.addTrigger(sampleCounter+j, i+1);
							triggerState[i] = true;
							printf("-!- Trigger at sample %i, channel %i\n", sampleCounter+j, i+1);
						}
					} else {
						triggerState[i] = false;
					}
				}
			}
			
			if (eventChain.count() > 0) {
				int err = clientrequest(ftSocket, sampleBlock.asRequest(), resp.in());
				if (err || !resp.checkPut()) {
					fprintf(stderr, "Could not write events to FieldTrip buffer\n.");
				}
			}
		}
		
		sampleCounter += block.numSamples;
	}
	
	BS.closeDevice();
	if (ftSocket != -1) close_connection(ftSocket);
	
	return 0;
}
