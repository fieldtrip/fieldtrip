/* 
 * Copyright (C) 2010, Stefan Klanke
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 * This program is used for setting up a shared memory segment where CTFs Acq can
 * write data to. This program picks up the data, copies it into a short internal
 * ringbuffer, and subsequently frees up the corresponding slot in the shared memory
 * segment. Notably this should also work in those cases where Acq writes more samples
 * than would normally fit into each slot, that is, this program also cleans up the 
 * next slot in those cases, such that Acq finds a perfectly valid slot for the next
 * data packet.
 *
 * Using a socket pair, a separate thread is notified of the incoming data packets,
 * which are then being pulled out of the internal ringbuffer, analysed for triggers,
 * and then written to a local or remote FieldTrip buffer (samples + events).
 * 
 * The program is stopped by pressing CTRL-C after which a clean shutdown is performed.
 */
#include <stdio.h>
#include "ctf.h"
#include <buffer.h>
#include <sys/shm.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <unistd.h>
#include <signal.h>
#include <pthread.h>

#define CTF_TRIGGER_TYPE   11

#define OVERALLOC   1000    /* to overcome Acq bug */
#define INT_RB_SIZE   10	/* Length of internal ring buffer of (overallocated) packets */
#define MAX_CHANNEL  512
#define MAX_TRIGGER   16

#ifdef FAKE_MEG
/* modify the definitions for testing on different machines with fake_meg */
#undef  ACQ_MSGQ_SIZE   
#define ACQ_MSGQ_SIZE 10
#undef  ACQ_MSGQ_SHMKEY 
#define ACQ_MSGQ_SHMKEY    0x08150842 
#endif

static char usage[] = 
"Usage: acq2ft [hostname port]\n" \
"Just calling 'acq2ft' without parameters starts a local FieldTrip buffer on port 1972.\n" \
"Using '-' for the hostname starts a local buffer on the given port.\n";

/*  This is used for the internal ringbuffer to overcome the bug in Acq.
    Same structure but bigger data field...
*/
typedef struct {
  ACQ_MessageType message_type;
  int messageId;
  int sampleNumber;
  int numSamples;
  int numChannels;
  int data[28160 + OVERALLOC];
} ACQ_OverAllocType;

/* This is used internally to store events in the same format that is being sent to the buffer.
   The memory pointed to by 'evs' is not free'd between handling different slots, but instead 
   reused. 'sizeAlloc' keeps track of the amount of reserved memory, whereas 'size' contains
   the size in bytes of the actual event content. */
typedef struct {
	void *evs;
	int num, size, sizeAlloc;
} EventChain;

volatile int keepRunning = 1;
int ftSocket = -1;
int mySockets[2];	/* a socket pair for communication between the two threads */

ACQ_OverAllocType intPackets[INT_RB_SIZE];


int sensType[MAX_CHANNEL];
int triggerChannel[MAX_TRIGGER];
char triggerChannelName[MAX_TRIGGER][32];
int triggerChannelNameLen[MAX_TRIGGER];
int lastValue[MAX_TRIGGER];
int numTriggerChannels;


/* prototypes for helper functions defined below */
ACQ_MessagePacketType *createSharedMem();
void closeSharedMem(ACQ_MessagePacketType *packet);
void initSharedMem(ACQ_MessagePacketType *packet);
ACQ_MessageType waitPacket(volatile ACQ_MessageType *msgtyp);
void abortHandler(int sig);
void *dataToFieldTripThread(void *arg);
headerdef_t *handleRes4(const char *dsname, UINT32_T *size);
void addTriggerEvent(EventChain *EC, int trigChan, UINT32_T sample, int value);


int main(int argc, char **argv) {
	ACQ_MessagePacketType *packet;
	int currentPacket, intSlot;
	int lastId;
	int numChannels, numSamples, sampleNumber;
	host_t host;
	pthread_t tcpserverThread, convertThread;
	
	/* Parse command line arguments */
	if (argc == 1) {
		/* Put defaults in case no arguments given */
		strcpy(host.name, "-");
		host.port = 1972;
	} else if (argc == 3) {
		strncpy(host.name, argv[1], sizeof(host.name));
		host.name[sizeof(host.name)-1]=0;
		
		host.port = atoi(argv[2]);
		if (host.port <= 0) {
			fputs("Port number must be positive\n\n",stderr);
			fputs(usage, stderr);
			return 1;
		}
	} else {
		fputs(usage, stderr);
		return 1;
	}
	
	/* Create socket pair */
	if (socketpair(AF_UNIX, SOCK_STREAM, 0, mySockets)) {
		perror("Could not create socket pair for internal communication: ");
		return 1;
	}
	
	/*  Internal ringbuffer uses same mechanism as Acq: free slots are
		marked with ACQ_MSGQ_INVALID. To avoid polling in the converter thread,
		the main thread sends the index of the latest slot over the socket pair.
		If the converter thread does not manage to clear the ringbuffer quickly
		enough, the program will report an error. */
	for (intSlot=0;intSlot<INT_RB_SIZE;intSlot++) {
		intPackets[intSlot].message_type = ACQ_MSGQ_INVALID;
	}
	
	/* Spawn conversion thread */
	if (pthread_create(&convertThread, NULL, dataToFieldTripThread, NULL)) {
		fprintf(stderr, "Could not spawn conversion thread\n");
		return 1;
	}
	
	/* Spawn tcpserver or connect to remote buffer */
	if (strcmp(host.name, "-") == 0) {
		int rc;
		
		rc = pthread_create(&tcpserverThread, NULL, tcpserver, &host);
		if (rc) {
			fprintf(stderr, "Could not spawn tcpserver thread (%d)\n", rc);
			return 1;
		}
		ftSocket = 0; /* dma */
	} else {
		ftSocket = open_connection(host.name, host.port);
		
		if (ftSocket < 0) {
			fprintf(stderr, "Could not connect to remote buffer\n");
			return 1;
		}
	}		
	
	printf("Trying to set up shared memory...\n");
	packet = createSharedMem();
	if (packet == NULL) {
		fprintf(stderr, "Could not set up shared memory - exiting\n");
		exit(1);
	}
	initSharedMem(packet);
	
	/* register CTRL-C handler */
	signal(SIGINT, abortHandler);
	
	/* Both ringbuffers start at 0 */
	currentPacket = intSlot = 0;
	
	printf("Entering main loop. Press CTRL-C to stop operation.\n\n");
	while (keepRunning) {
		int size;
		char *dataset;
		int code = waitPacket(&packet[currentPacket].message_type);
		
		switch(code) {
			case ACQ_MSGQ_SETUP_COLLECTION:
				lastId = packet[currentPacket].messageId;
				dataset = (char *) packet[currentPacket].data;
				printf("Setup | ID=%i | %.255s\n", lastId, dataset);
				
				size = 5*sizeof(int) + strlen(dataset) + 1;
				if (size >= sizeof(ACQ_MessagePacketType)) {
					fprintf(stderr, "Dataset name not zero-terminated -- ignoring\n");
					continue;
				}
				
				if (intPackets[intSlot].message_type != ACQ_MSGQ_INVALID) {
					fprintf(stderr, "Internal converter thread does not keep up with the load.\n");
				} else {
					memcpy(&intPackets[intSlot], &packet[currentPacket], size);
					write(mySockets[0], &intSlot, sizeof(int));
					if (++intSlot == INT_RB_SIZE) intSlot=0;
				}
				
				packet[currentPacket].message_type = ACQ_MSGQ_INVALID;
				if (++currentPacket == ACQ_MSGQ_SIZE) currentPacket=0;
				break;
			case ACQ_MSGQ_DATA:
				lastId = packet[currentPacket].messageId;
				numChannels  = packet[currentPacket].numChannels;
				numSamples   = packet[currentPacket].numSamples;
				sampleNumber = packet[currentPacket].sampleNumber;
				printf("Data | %3i channels x %3i samples | nr = %6i | ID=%4i | slot=%3i\n", numChannels, numSamples, sampleNumber, lastId, currentPacket);
								
				size = 5*sizeof(int) + numSamples*numChannels*sizeof(int);
				
				if (size > sizeof(ACQ_OverAllocType)) {
					fprintf(stderr, "Acq wrote far too much data -- cannot handle this\n");
				} else {
					if (intPackets[intSlot].message_type != ACQ_MSGQ_INVALID) {
						fprintf(stderr, "Internal converter thread does not keep up with the load.\n");
					} else {
						memcpy(&intPackets[intSlot], &packet[currentPacket], size);
						write(mySockets[0], &intSlot, sizeof(int));
						if (++intSlot == INT_RB_SIZE) intSlot=0;
					}
				}				
				
				if (numSamples * numChannels > 28160) {
					fprintf(stderr, "Warning: Acq wrote too much data into this block!\n");
            
					/* clear this packet */
					packet[currentPacket].message_type = ACQ_MSGQ_INVALID;
					/* next packet last in ringbuffer? */ 
					if (++currentPacket == ACQ_MSGQ_SIZE) {
						currentPacket = 0;
					} else {
						/* make sure the next packet is marked as free */
						packet[currentPacket].message_type = ACQ_MSGQ_INVALID;
					}
				} else {
					packet[currentPacket].message_type = ACQ_MSGQ_INVALID;
					if (++currentPacket == ACQ_MSGQ_SIZE) currentPacket=0;
				}
				break;
			case ACQ_MSGQ_CLOSE_CONNECTION:
				printf("Stop | ID=%i \n", lastId);
				
				/* Does this mean anything for the buffer as well ? */
				
				packet[currentPacket].message_type = ACQ_MSGQ_INVALID;
				currentPacket = 0;
				break;
			case ACQ_MSGQ_INVALID:
				printf("Waiting (on %i)\n",currentPacket);
				break;
			default:
				fprintf(stderr, "Warning: Unexpected ID in packet %i\n", currentPacket);
				break;
		}
	}
	printf("Closing sockets / stopping tcpserver...\n");
	close(mySockets[0]); /* the other one will be closed in the thread */
	if (ftSocket > 0) {
		close_connection(ftSocket);
	} else {
		pthread_cancel(tcpserverThread);
		pthread_detach(tcpserverThread);
	}
	printf("Closing shared memory...\n");
	closeSharedMem(packet);
	printf("Joining conversion thread...\n");
	pthread_join(convertThread, NULL);
	printf("Done.\n");
	return 0;
}


ACQ_MessagePacketType *createSharedMem() {
	void *ptr;
	int shmid;
	size_t siz = sizeof(ACQ_MessagePacketType)*ACQ_MSGQ_SIZE;
   
    siz += OVERALLOC*sizeof(int); /* to overcome Acq bug */
		
	shmid = shmget(ACQ_MSGQ_SHMKEY, siz, 0666|IPC_CREAT);
	if (shmid == -1) {
		perror("shmget");
		return NULL;
	}

	/* attach to the segment to get a pointer to it */
	ptr = shmat(shmid, NULL, 0);
	if (ptr == (void *) -1) {
		perror("shmat");
		return NULL;
	}
	return (ACQ_MessagePacketType *) ptr;
}

void closeSharedMem(ACQ_MessagePacketType *packet) {
	shmdt(packet);
}

void initSharedMem(ACQ_MessagePacketType *packet) {
	int i;
	for (i=0;i<ACQ_MSGQ_SIZE;i++) {
		packet[i].message_type =  ACQ_MSGQ_INVALID;
	}	
}

/* wait up to 1 second for new packet */
ACQ_MessageType waitPacket(volatile ACQ_MessageType *msgtyp) {
	int i;
	volatile ACQ_MessageType t;
	for (i=0;i<1000 && keepRunning;i++) {
		t = *msgtyp;
		if (t != ACQ_MSGQ_INVALID) break;
		usleep(1000);
	}
	return t;
}

/** Sets 'keepRunning' to 0 in case the uses presses CTRL-C - this makes
    the main thread leave its processing loop, after which the program
	exits smoothly.
*/
void abortHandler(int sig) {
	printf("Ctrl-C pressed -- stopping acq2ft...\n");
	keepRunning = 0;
}

/** Background thread for grabbing setup and data packets from the internal, overallocated ringbuffer.
    This thread stops automatically if the socket pair is closed in the main thread.
*/
void *dataToFieldTripThread(void *arg) {
	ACQ_OverAllocType *pack;
	EventChain EC;
	int slot, res, i,j;
	int numChannels = 0, numSamples;
	int warningGiven = 0;
	messagedef_t reqdef;
	message_t request, *response;
	
	EC.evs = NULL;
	EC.sizeAlloc = 0;
	
	request.def = &reqdef;
	
	while (1) {
		res = read(mySockets[1], &slot, sizeof(int));
		
		if (res == 0) break; /* socket pair was closed - exit */
		
		if (res != sizeof(int)) {
			fprintf(stderr, "Error when reading from socket pair\n");
			break;
		}
		
		if (slot<0 || slot>=INT_RB_SIZE) {
			fprintf(stderr, "Got errorneous slot number from socket pair\n");
			break;
		}
		
		pack = &intPackets[slot];
		
		if (pack->message_type == ACQ_MSGQ_SETUP_COLLECTION) {
			UINT32_T size;
			headerdef_t *hdef;	/* contains header information + chunks !!! */
			
			hdef = handleRes4((const char *) pack->data, &size);
			/* clear internal ringbuffer slot */
			pack->message_type = ACQ_MSGQ_INVALID;
			
			if (hdef == NULL) continue;	/* problem while picking up header -- ignore this packet */
			
			/* prepare PUT_HDR request here, but only send it along with first data block */
			reqdef.version = VERSION;
			reqdef.command = PUT_HDR;
			reqdef.bufsize = size;
			request.buf = hdef;
				
			res = clientrequest(ftSocket, &request, &response);
				
			if (res < 0) {
				fprintf(stderr, "Error in FieldTrip connection\n");
			} else if (response) {
				if (response->def->command != PUT_OK) {
					fprintf(stderr, "Error in PUT_HDR\n");
				} else {
					/* printf("FT: Transmitted header\n"); */
					
					/* set numChannels variable to the value we picked up
					   this also enables transmitting data in following packets
					*/
					numChannels = hdef->nchans;
					warningGiven = 0;
				}
				cleanup_message((void **) &response);
			}
			free(hdef);
			
		} else if (pack->message_type == ACQ_MSGQ_DATA) {
			datadef_t *ddef; 
			int sampleNumber;
			
			if (numChannels == 0) {
				fprintf(stderr, "No header written yet -- ignoring data packet\n");
				pack->message_type = ACQ_MSGQ_INVALID;
				continue;
			}
			
			if (!warningGiven && pack->numChannels != numChannels) {
				printf("\nWARNING: Channel count in first data packet does not equal header information from .res4 file (%i channels)\n\n", numChannels);
				warningGiven = 1;
			}
			
			sampleNumber = pack->sampleNumber;
			numSamples   = pack->numSamples;
			
			/* Put the FT datadef at the location of the current ACQ packet definition.
			   This just fits, no memcpy'ing of the samples again... */
			ddef = (datadef_t *) &pack->messageId; 
			
			ddef->nsamples = numSamples;
			ddef->nchans = numChannels;
			ddef->data_type = DATATYPE_INT32;
			ddef->bufsize = sizeof(int) * numSamples * numChannels;
			
			reqdef.version = VERSION;
			reqdef.command = PUT_DAT;
			reqdef.bufsize = ddef->bufsize + sizeof(datadef_t);
			request.buf = ddef; /* data is still behind that */
			
			res = clientrequest(ftSocket, &request, &response);
			
			if (res < 0) {
				fprintf(stderr, "Error in FieldTrip connection (writing data)\n");
			} else if (response) {
				if (response->def->command != PUT_OK) {
					fprintf(stderr, "Error in PUT_DAT\n");
				} else {
					/* printf("FT: Transmitted samples\n"); */
				}
				cleanup_message((void **) &response);
			}
			
			/* look at trigger channels and add events to chain, clear this first */
			EC.size = EC.num = 0;
			for (j=0;j<numSamples;j++) {
				int *sj = pack->data + j*numChannels;
				for (i=0;i<numTriggerChannels;i++) {
					int sji = sj[triggerChannel[i]];
					if (sji != lastValue[i] && sji > 0) addTriggerEvent(&EC, i, sampleNumber + j, sji);
					lastValue[i] = sji;
				}
			}
			if (EC.size > 0) {
				reqdef.version = VERSION;
				reqdef.command = PUT_EVT;
				reqdef.bufsize = EC.size;
				request.buf = EC.evs;
				
				res = clientrequest(ftSocket, &request, &response);
			
				if (res < 0) {
					fprintf(stderr, "Error in FieldTrip connection (writing events)\n");
				} else if (response) {
					if (response->def->command != PUT_OK) {
						fprintf(stderr, "Error in PUT_EVT\n");
					} else {
						printf("Wrote %i events (%i bytes)\n", EC.num, reqdef.bufsize);
						/* printf("FT: Transmitted samples\n"); */
					}
				}
				cleanup_message((void **) &response);
			}
		} else {
			fprintf(stderr,"Converter thread: Packet contains neither SETUP nor DATA (%i)...\n", pack->message_type);
		}
		pack->message_type = ACQ_MSGQ_INVALID;
	}
	printf("Leaving converter thread...\n");
	close(mySockets[1]);
	if (EC.sizeAlloc > 0) free(EC.evs);
	return NULL;
}

short getInt16(char *src) {
	union {
		short s;
		char b[2];
	} u;
	u.b[0] = src[1];
	u.b[1] = src[0];
	return u.s;
}

int getInt32(char *src) {
	union {
		int  i;
		char b[4];
	} u;
	u.b[0] = src[3];
	u.b[1] = src[2];
	u.b[2] = src[1];
	u.b[3] = src[0];
	return u.i;
}

/** Parses .res4 file corresponding to dsname and returns a headerdef and
	chunks suitable for sending off to FieldTrip. Also returns size of that
	message payload in 'size'. Returns NULL on error.
	
	Returned pointer needs to be disposed of using free() later.
*/
headerdef_t *handleRes4(const char *dsname, UINT32_T *size) {
	headerdef_t *hdef;
	char *res4name;
	ft_chunk_t *chunk;
	int i, len, pdot, pstart, rlen, offset, offsetNames, numFilters;
	FILE *f;
	union {
		double d;
		char b[8];
	} fsamp;
				
	/* printf("Picked up dataset name : %s\n", dsname); */
	len = strlen(dsname);
			
	pdot = len;
	while (--pdot>0) {
		if (dsname[pdot] == '.') break;
	}
	if (pdot == 0) {
		fprintf(stderr, "No . found in dataset name - don't know where to pick up .res4 file\n");
		return NULL;
	}
	pstart = pdot;
	while (--pstart>0) {
		if (dsname[pstart] == '/') {
			/* slash at pstart -> increase because filename comes 1 character later */
			pstart++;
			break;
		}
	}
	/* no slash found? then pstart = 0, which is fine - treat as relative path */
			
	/* compose .res4 file name from dsname/dsname[pslash+1:pdot].res4 */
	rlen = len + 1 + pdot - pstart + 5 + 1;
			
	res4name = (char *) malloc(rlen);
	if (res4name == NULL) {
		fprintf(stderr, "Out of memory -- could not compose filename\n");
		return NULL;
	}
			
	/* compose .res4 file name from dsname/dsname[pstart:pdot-1].res4 */
	memcpy(res4name, dsname, len);
	res4name[len] = '/';
	memcpy(res4name + len + 1, dsname + pstart, pdot - pstart);
	rlen = len + 1 + pdot - pstart;
	res4name[rlen++] = '.';
	res4name[rlen++] = 'r';
	res4name[rlen++] = 'e';
	res4name[rlen++] = 's';
	res4name[rlen++] = '4';
	res4name[rlen++] = 0;
			
	/* printf("Trying to open %s\n", res4name); */
	f = fopen(res4name, "rb");
	if (f==NULL) {
		fprintf(stderr, "File %s could not be opened\n", res4name);
		free(res4name);
		return NULL;
	}
			
	fseek(f, 0, SEEK_END);
	len = ftell(f);
	fseek(f, 0, SEEK_SET);
			
	printf("\nCTF RES4 file %s contains %i bytes.\n", res4name, len);
	free(res4name); /* not needed anymore */

	/* get space for headerdef, 1x chunkdef, size of res4-file */
	rlen = sizeof(headerdef_t) + sizeof(ft_chunkdef_t) + len;
	hdef = (headerdef_t *) malloc(rlen);
	
	if (hdef == NULL) {
		fprintf(stderr, "Out of memory - can not allocate space for reading .res4 file\n");
		fclose(f);
		return NULL;
	}
					
	chunk = (ft_chunk_t *) (hdef+1); /* chunk starts directly after hdef */
	if (fread(chunk->data, 1, len, f) != len) {
		fprintf(stderr, "Could not read complete .res4 file\n");
		fclose(f);
		free(hdef);
		return NULL;
	} 
	fclose(f);			
	
	hdef->nsamples = 0;
	hdef->nevents = 0;
	/* .res4 file is big-endian, but we assume this machine to be little-endian */
	hdef->nchans  = getInt16(chunk->data + 1292);
	
	printf("Number of channels: %i\n", hdef->nchans);
	
	for (i=0;i<8;i++) {
		fsamp.b[7-i] = chunk->data[1296+i];
	}
	hdef->fsample = (float) fsamp.d;
	
	if (hdef->fsample < 1.0) {
		printf("\nWARNING: suspicious sampling frequency (%f Hz) picked from .res4 -- setting to 1200 Hz\n\n", hdef->fsample);
		hdef->fsample = 1200.0;
	} else {
		printf("Sampling frequency: %f Hz\n", hdef->fsample);
	}
	
	hdef->data_type = DATATYPE_INT32;
	hdef->bufsize = sizeof(ft_chunkdef_t) + len;
	chunk->def.type = FT_CHUNK_CTF_RES4;
	chunk->def.size = len;
	if (size != NULL) *size = rlen;

	/* Ok, the header + chunk is ready, now parse some additional information
		used internally
	*/
	
	/* set rlen to "run description length", see read_ctf_res4, line 92f */
	rlen = getInt32(chunk->data + 1836);
	
	/* offset points at first byte after run_desc */
	offset = rlen + 1844;
	
	/* "number of filters" */
	numFilters = getInt16(chunk->data + offset);
	offset += 2;
	
	/* printf("numFilters = %i\n", numFilters); */
	
	for (i=0;i<numFilters;i++) {
		len = getInt16(chunk->data + offset+16);
		offset += 18 + len;
		/* printf("Filter %i: offset = %i,  len = %i\n", i, offset, len); */
	}
	
	/* next we've got 32 bytes per channel (name) */
	offsetNames = offset;
	offset += 32 * hdef->nchans; 
	
	numTriggerChannels = 0;
	
	for (i=0;i<hdef->nchans;i++) {
		unsigned char ct = chunk->data[offset + 1 + 1328*i]; 
		if (i<MAX_CHANNEL) {
			sensType[i] = ct;
		}
		if (ct == CTF_TRIGGER_TYPE) {
			if (numTriggerChannels < MAX_TRIGGER) {
				int j;
				char *tn = triggerChannelName[numTriggerChannels];
				
				memcpy(tn, chunk->data + offsetNames + 32*i, 32);
				for (j=0;j<32;j++) {
					if (tn[j] <= 32 || tn[j] == '-') break;
				}
				triggerChannelNameLen[numTriggerChannels] = j;
			
				triggerChannel[numTriggerChannels] = i;
				printf("Trigger channel @ %i: %.*s\n", i, j, tn);
				
				++numTriggerChannels;
			}
		}	
	}
	
	for (i=0;i<numTriggerChannels;i++) {
		lastValue[i] = 0;
	}
	
	return hdef;
}

/** Add one event to the chain, increases 'size' and 'num' fields of EC if successful.
	If memory can not be allocated, the event is dropped and an error message is printed.
*/
void addTriggerEvent(EventChain *EC, int trigChan, UINT32_T sample, int value) {
	eventdef_t *ne;
	char *ntype;
	int addSize = sizeof(eventdef_t) + triggerChannelNameLen[trigChan] + sizeof(int);
	
	if (EC->sizeAlloc == 0) {
		EC->evs = malloc(addSize);
		if (EC->evs == NULL) {
			fprintf(stderr, "Cannot add trigger event - out of memory...\n");
			return;
		}
		EC->sizeAlloc = addSize;
	} else if (EC->size + addSize > EC->sizeAlloc) {
		/* try to get more space */
		void *nevs =  realloc(EC->evs, EC->size + addSize);
		if (nevs == NULL) {
			fprintf(stderr, "Cannot add trigger event - out of memory...\n");
			return;
		}
		EC->evs = nevs;
		EC->sizeAlloc = EC->size + addSize;
	}
	
	ntype = (char *) EC->evs + EC->size;
	
	EC->size+=addSize;
	EC->num++;
	
	ne = (eventdef_t *) ntype;
	ntype += sizeof(eventdef_t);
	
	ne->type_type   = DATATYPE_CHAR;
	ne->type_numel  = triggerChannelNameLen[trigChan];
	ne->value_type  = DATATYPE_INT32;
	ne->value_numel = 1;
	ne->sample = sample;
	ne->offset = 0;
	ne->duration = 0;
	ne->bufsize = triggerChannelNameLen[trigChan] + sizeof(int);
	memcpy(ntype, triggerChannelName[trigChan], triggerChannelNameLen[trigChan]);
	memcpy(ntype + triggerChannelNameLen[trigChan], &value, sizeof(int));
}
