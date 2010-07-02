/*
 * Copyright (C) 2010, Stefan Klanke
 * Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
 * Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 */
#include <stdio.h>
#include <math.h>
#include <signal.h>
#include <buffer.h>
#include <FtBuffer.h>
#include <MatlabNetComClient.h>
#include <pthread.h>


#define MAX_OBJ  5000
#define MAX_LEN  100
#define RB_SIZE  50

int keepRunning = 1;
char hostname[256];
int port, numChannels, blockSize;
int ftSocket = -1;
int maxCscSamples, recordBufferSize, maxEventStringLength;
HANDLE pipeWrite, pipeRead;

char cheetahObjects[MAX_OBJ][MAX_LEN];
char cheetahTypes[MAX_OBJ][MAX_LEN];
int cscIndex[MAX_OBJ];
int evcIndex[MAX_OBJ];
int numCSC = 0, numEvC = 0;
int fSamp = -1;

char *channelNamesForChunk = NULL;
int channelNamesSize;
short *retBufSamples;
uint64_t *retBufTimeStamps;
int *retBufChannelNumbers;
int *retBufValidSamples;
int *retBufSamplingFreq;
int *retBufEventIDs;
int *retBufTTL;
char **retBufEventString;


typedef struct {
	int status; 	// 0=empty, 1=being filled, 2=being sent
	int chansFilled;
	uint64_t timeStamp;
	short *samples; /* numChannels x maxCscSamples */
} MultiChannelBlock;

typedef struct {
	int status;		// 0=empty, 1=being sent
	uint64_t timeStamp;
	char *name;
	int eventID;
	int ttl;
} EventBlock;

int *rbStart; /* start of ringbuffer per channel */
int *rbNumEl; /* number of elements in RB per channel */
MultiChannelBlock mcbRingBuf[RB_SIZE];
int evStart = 0, evNumEl = 0;
EventBlock evRingBuf[RB_SIZE];

int getObjectsAndTypes() {
	char *objPtr[MAX_OBJ], *typPtr[MAX_OBJ];
	int n;
	
	for (n=0;n<MAX_OBJ;n++) {
		cheetahObjects[n][0] = 0;
		cheetahTypes[n][0] = 0;
		objPtr[n] = cheetahObjects[n];
		typPtr[n] = cheetahTypes[n];
	}
	
	numCSC = 0;
	
	n = GetCheetahObjectsAndTypes(objPtr, typPtr);
	if (n==0) return -1;
	
	channelNamesSize = 0;
	
	for (n=0;n<MAX_OBJ;n++) {
		if (cheetahObjects[n][0] == 0) break;
		printf("%3i: %-20s  %s\n", n, cheetahObjects[n], cheetahTypes[n]);
		
		if (!strcmp(cheetahTypes[n], "CscAcqEnt")) {
			cscIndex[numCSC++] = n;
			channelNamesSize += strlen(cheetahObjects[n]) + 1; // trailing 0
			continue;
		}
		
		if (!strcmp(cheetahTypes[n], "EventAcqEnt")) {
			evcIndex[numEvC++] = n;
			continue;
		}
	}
	
	if (numCSC>0) {
		int offset = 0;
		channelNamesForChunk = new char[channelNamesSize];
		for (int i=0;i<numCSC;i++) {
			const char *name = cheetahObjects[cscIndex[i]];
			int len = strlen(name);
			memcpy(channelNamesForChunk + offset, name, len+1); // includes \0
			offset += len+1;
		}
	}

	return n;
}


void *dataToFieldTripThread(void *arg) {
	int slot;
	DWORD numRead;
	FtBufferRequest req;
	FtBufferResponse resp;	
	
	printf("Thread is running...\n");
	
	while (1) {
		ReadFile(pipeRead, &slot, sizeof(slot), &numRead, NULL);
		if (numRead == 0) continue;
		if (numRead != sizeof(slot)) break;
		
		if (slot < 0 || slot >= RB_SIZE) break;
		
		printf("Sending out slot %i...\n", slot);
		
		req.prepPutData(numCSC, maxCscSamples, DATATYPE_INT16, mcbRingBuf[slot].samples);
		int r = tcprequest(ftSocket, req.out(), resp.in());
		if (r < 0 || !resp.checkPut()) {
			printf("Could not write samples to FieldTrip buffer - aborting\n");
			break;
		}						
		mcbRingBuf[slot].status = 0;
	}
	
	printf("Leaving thread...\n");
	
	return NULL;
}


void abortHandler(int sig) {
	printf("Ctrl-C pressed -- stopping...\n");
	keepRunning = 0;
}

int main(int argc, char *argv[]) {
	FtBufferRequest req;
	FtBufferResponse resp;
	int numRecordsReturned, numRecordsDropped;
	int r;
	pthread_t sendThread;
	
	if (argc<2) {
		strcpy(hostname, "localhost");
	} else {
		strncpy(hostname, argv[1], 256);
		hostname[255] = 0;
	}
	if (argc<3) {
		port = 1972;
	} else {
		port = atoi(argv[2]);
	}
	
	r = CreatePipe(&pipeRead, &pipeWrite, NULL, 0);
	if (!r) {
		printf("Could not create internal pipe\n");
		exit(1);
	}
	
	printf("Trying to connect to FieldTrip buffer on %s:%i\n", hostname, port);
	ftSocket = open_connection(hostname, port);
	if (ftSocket < 0) {
		printf("Connection failed - aborting\n");
		exit(1);
	}
	
	printf("Trying to connect to Neuralynx server on localhost\n");
	r = ConnectToServer("localhost");
	
	if (!r) {
		printf("Connection failed\n");
		return 1;
	}
	printf("Connected!\n");
	
	r = SetApplicationName("Neuralynx-to-FieldTrip");
	
	r = getObjectsAndTypes();
	if (r==-1) {
		printf("Could not get Cheetah objects and types - aborting\n");
		exit(1);
	}
	
	recordBufferSize = GetRecordBufferSize();
	maxCscSamples = GetMaxCSCSamples();
	maxEventStringLength = GetMaxEventStringLength();
	
	printf("Number of CSC streams......: %i\n", numCSC);
	printf("Number of event streams....: %i\n", numEvC);
	printf("Samples per block..........: %i\n", maxCscSamples);
	printf("Record buffer size.........: %i\n", recordBufferSize);
	printf("Max. length event strings..: %i\n", maxEventStringLength);
	printf("Server PC name.............: %s\n", GetServerPCName());
	printf("Server IP address..........: %s\n", GetServerIPAddress());
	printf("Server application.........: %s\n", GetServerApplicationName());
	
	if (numCSC==0) {
		printf("No continuous channels - exiting\n");
		exit(1);
	}
	
	printf("Allocationg temporary buffers...\n");
	
	retBufSamples         = new short[recordBufferSize * maxCscSamples];
	retBufTimeStamps      = new uint64_t[recordBufferSize];
	retBufChannelNumbers  = new int[recordBufferSize];
	retBufSamplingFreq    = new int[recordBufferSize];	
	retBufValidSamples    = new int[recordBufferSize];
	retBufEventIDs        = new int[recordBufferSize];
	retBufTTL             = new int[recordBufferSize];
	retBufEventString     = new char*[recordBufferSize];
	for (int j=0;j<recordBufferSize;j++) {
		retBufEventString[j] = new char[maxEventStringLength+1];
	}
	rbStart = new int[numCSC];
	rbNumEl = new int[numCSC];
	
	printf("Preparing internal ring buffer\n");
	
	for (int j=0;j<RB_SIZE;j++) {
		mcbRingBuf[j].status = 0; // empty
		mcbRingBuf[j].timeStamp = 0;
		mcbRingBuf[j].chansFilled = 0;
		mcbRingBuf[j].samples = new short[numCSC * maxCscSamples];
		
		evRingBuf[j].status = 0;
		evRingBuf[j].timeStamp = 0;
		evRingBuf[j].name = new char[maxEventStringLength+1];
	}
	
	printf("Opening CSC channels...\n");
	for (int j=0;j<numCSC;j++) {
		rbStart[j] = rbNumEl[j] = 0;
		r = OpenStream(cheetahObjects[cscIndex[j]]);
		if (!r) {
			printf("Could not open CSC channel %s\n", cheetahObjects[cscIndex[j]]);
			keepRunning = false;
			break;
		}
	}

	printf("Opening event channels...\n");
	for (int j=0;j<numEvC;j++) {
		r = OpenStream(cheetahObjects[evcIndex[j]]);
		if (!r) {
			printf("Could not open event channel %s\n", cheetahObjects[evcIndex[j]]);
			keepRunning = false;
			break;
		}
	}
	
	/* register CTRL-C handler */
	signal(SIGINT, abortHandler);
	
	if (pthread_create(&sendThread, NULL, dataToFieldTripThread, NULL)) {
		fprintf(stderr, "Could not spawn conversion thread\n");
		return 1;
	}
	
	printf("Starting acquisition loop - press CTRL-C to stop.\n");	
	
	while (keepRunning) {
		int doneSomething = 0;
		
		// loop over continuous channels
		for (int nCSC = 0;nCSC < numCSC; nCSC++) {
			char *stream = cheetahObjects[cscIndex[nCSC]];
			
			r = GetNewCSCData(stream, &retBufTimeStamps, &retBufChannelNumbers, &retBufSamplingFreq, &retBufValidSamples, &retBufSamples, &numRecordsReturned, &numRecordsDropped);
		
			if (r==0) {
				printf("Error in GetNewCSCData from stream '%s' - aborting\n", stream);
				keepRunning = false;
				break;
			}
	
			if (numRecordsDropped > 0) {
				printf("Records dropped in stream '%s' - cancelling operation.\n", stream);
				keepRunning = false;
				break;
			}
			
			if (numRecordsReturned == 0) continue;
			
			doneSomething = 1;
	
			for (int i=0;i<numRecordsReturned;i++) {
				// TODO: check if using this channel number is right and preferable to nCSC
				int chan = retBufChannelNumbers[i];
				int freq = retBufSamplingFreq[i];
				int vals = retBufValidSamples[i]; 
				uint64_t ts = retBufTimeStamps[i];

				// printf("#%i/%i  chan = %3i  Time = %8u  fSamp = %6i  valid = %4i\n", i+1, numRecordsReturned, chan, (uint32_t) ts, freq, vals);
				if (fSamp == -1) {
					fSamp = freq;
					// write out header to Fieldtrip buffer here
					req.prepPutHeader(numCSC, DATATYPE_INT16, (float) fSamp);
					if (channelNamesForChunk != NULL) {
						req.prepPutHeaderAddChunk(FT_CHUNK_CHANNEL_NAMES, channelNamesSize, channelNamesForChunk);
					}
					int r = tcprequest(ftSocket, req.out(), resp.in());
					if (r < 0 || !resp.checkPut()) {
						printf("Could not write header to FieldTrip buffer - aborting\n");
						keepRunning = false;
						break;
					}						
				} else if (fSamp != freq) {
					printf("Inconsistent sampling frequency in stream '%s' - aborting!\n", stream);
					keepRunning = false;
					break;
				}
				if (vals != maxCscSamples) {
					printf("Got less than %i samples in record from '%s' - aborting\n", maxCscSamples, stream);
					keepRunning = false;
					break;
				}
				
				// ok, everything looks fine, add to ring buffer
				// this is the index where we need to put it, wrap around if necessary
				int index = rbStart[chan] + rbNumEl[chan];
				if (index>=RB_SIZE) index-=RB_SIZE; 
				
				if (mcbRingBuf[index].status == 2) {
					printf("Streaming thread cannot keep up with load - aborting\n");
					keepRunning = false;
					break;
				}
				if (mcbRingBuf[index].status == 0) {
					mcbRingBuf[index].timeStamp = ts;
					mcbRingBuf[index].status = 1;
				} else if (mcbRingBuf[index].timeStamp != ts) {
					printf("Inconsistent time stamps in stream '%s' records and internal ring buffer - aborting\n", stream);
					keepRunning = false;
					break;
				}
				// fine, write samples at right channel
				// have 'dest' point to proper row in ring buffer slot
				short *dest = mcbRingBuf[index].samples + chan;
				// have 'src' point to first sample in current record
				short *src  = retBufSamples + i*maxCscSamples;
				for (int n=0;n<maxCscSamples;n++) {
					dest[n*numCSC] = src[n];
				}
				rbNumEl[chan]++;
				// see if all channels are filled in this slot, and if so, write out
				if (++mcbRingBuf[index].chansFilled == numCSC) {
					DWORD numWritten = 0;
					
					printf("Slot %i filled!\n", index);
					mcbRingBuf[index].chansFilled = 0;
					mcbRingBuf[index].status = 2;
					
					// write index to pipe so the other thread picks it up
					int ok = WriteFile(pipeWrite, &index, sizeof(index), &numWritten, NULL);
					if (!ok || numWritten != sizeof(index)) {
						printf("Error when writing to pipe\n");
						keepRunning = false;
						break;
					}
					
					for (int j=0;j<numCSC;j++) {
						if (++rbStart[j] == RB_SIZE) rbStart[j] = 0;
						rbNumEl[j]--;
					}
				}
			}
		}  // end of loop over cont. channels
		
		for (int nEvC = 0; nEvC < numEvC; nEvC++) {
			char *stream = cheetahObjects[evcIndex[nEvC]];
			
			r = GetNewEventData(stream, &retBufTimeStamps, &retBufEventIDs, &retBufTTL, retBufEventString, &numRecordsReturned, &numRecordsDropped);
			
		
			if (r==0) {
				printf("Error in GetNewEventData from stream '%s' - aborting\n", stream);
				keepRunning = false;
				break;
			}
	
			if (numRecordsDropped > 0) {
				printf("Records dropped in stream '%s' - cancelling operation.\n", stream);
				keepRunning = false;
				break;
			}
			
			if (numRecordsReturned == 0) continue;
			
			doneSomething = 1;

			for (int i=0;i<numRecordsReturned;i++) {
				printf("Received event %s\n", retBufEventString[i]);
				printf("TTL: %i,  Timestamp: %i,  ID: %i\n", retBufTTL[i], (int) retBufTimeStamps[i], retBufEventIDs[i]);
			}
			
			// TODO: write out events (ringbuffer first)

		}

		// if we're not to busy, let other processes use the machine
		if (!doneSomething) Sleep(1);
	} 
	
	DisconnectFromServer();

	if (ftSocket != -1) close_connection(ftSocket);
		
	delete[] retBufSamples;
	delete[] retBufTimeStamps;
	delete[] retBufChannelNumbers;
	delete[] retBufSamplingFreq;
	delete[] retBufValidSamples;
	delete[] retBufTTL;
	delete[] retBufEventIDs;
	delete[] rbStart;
	delete[] rbNumEl;
	for (int j=0;j<RB_SIZE;j++) {
		delete[] mcbRingBuf[j].samples;
		delete[] evRingBuf[j].name;
	}
	if (channelNamesForChunk) {
		delete[] channelNamesForChunk;
	}
	
	// this is for closing the background thread
	//ftSocket = -1;
	//WriteFile(pipeWrite, &ftSocket, sizeof(ftSocket), NULL, NULL);
		
	return 0;
}
