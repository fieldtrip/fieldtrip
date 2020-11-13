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
#include <sys/time.h>
#include <unistd.h>
#include <signal.h>
#include <pthread.h>
#include <math.h>

#define CTF_TRIGGER_TYPE   11

#define OVERALLOC   1000    /* to overcome Acq bug */
#define INT_RB_SIZE   10    /* Length of internal ring buffer of (overallocated) packets */
#define MAX_CHANNEL  512
#define MAX_TRIGGER   16
#define MAX_OUT        4

#ifdef FAKE_MEG
/* modify the definitions for testing on different machines with fake_meg */
#undef  ACQ_MSGQ_SIZE   
#define ACQ_MSGQ_SIZE 10
#undef  ACQ_MSGQ_SHMKEY 
#define ACQ_MSGQ_SHMKEY    0x08150842 
#endif

static char usage[] = 
"Usage: acq2ftx hostname:port:flags:decimation:channels\n";

/*  This is used for the internal ringbuffer to overcome the bug in Acq.
    Same structure but bigger data field...
    */
typedef struct {
	ACQ_MessageType message_type;
	int messageId;
	int sampleNumber;
	int numSamples;
	int numChannels;	/* Note: this is not the proper number of channels - bug in Acq */
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

typedef struct {
	int ftSocket;
	int writeEvents;
	int downSample;
	int writeRes4;
	int applyGains;
	int numChannelsFound, numChannelsGiven;
	char *channelList;
	char *channelNames[MAX_CHANNEL];
	int channels[MAX_CHANNEL];
	int skipSamples;
	int headerOk;
	float fSample;
} OutputConfig;



volatile int keepRunning = 1;
int ownBuffer = 0;
int mySockets[2];	/* a socket pair for communication between the two threads */

ACQ_OverAllocType intPackets[INT_RB_SIZE];
OutputConfig outConf[MAX_OUT];
int numOutputs = 0;

unsigned char channelSensType[MAX_CHANNEL];
double channelGainV[MAX_CHANNEL];
char channelNames[MAX_CHANNEL][32];
int channelNameLen[MAX_CHANNEL];

int numTriggerChannels;
int triggerChannel[MAX_TRIGGER];
char *triggerChannelName[MAX_TRIGGER];	/* points into channelNames */
int triggerChannelNameLen[MAX_TRIGGER];
int lastValue[MAX_TRIGGER];

struct {
	datadef_t ddef;
	union {
		int data[28160 + OVERALLOC];
		float fData[28160 + OVERALLOC];
	};
} putDatBuf;

/* prototypes for helper functions defined below */
ACQ_MessagePacketType *createSharedMem();
void closeSharedMem(ACQ_MessagePacketType *packet);
void initSharedMem(ACQ_MessagePacketType *packet);
ACQ_MessageType waitPacket(volatile ACQ_MessageType *msgtyp, double *time);
void abortHandler(int sig);
void alarmHandler(int sig);
void *dataToFieldTripThread(void *arg);
ft_chunk_t *handleRes4(const char *dsname, int *numChannels, float *fSample);
void addTriggerEvent(EventChain *EC, int trigChan, UINT32_T sample, int value);
void matchChannels(int numChannels);
void writeHeader(OutputConfig *oc, int numChannels, float fSample, const ft_chunk_t *res4chunk);
void writeSamples(OutputConfig *oc, const ACQ_OverAllocType *pack);
void writeEvents(OutputConfig *oc, const EventChain *EC, EventChain *ECd);

int main(int argc, char **argv) {
	ACQ_MessagePacketType *packet;
	int currentPacket, intSlot;
	int lastId;
	int numChannels, numSamples, sampleNumber;
	pthread_t tcpserverThread, convertThread;
	int i, rc;
	double sumT, sumT2;
	int numT;
	double tNow, tLast;
	struct timeval timerInterval = {0, 10000}; /*  10ms */
	struct timeval timerValue    = {0, 10000}; /*  10ms */
	struct itimerval timerOpts;
	sigset_t   signal_mask;

	timerOpts.it_value    = timerValue;
	timerOpts.it_interval = timerInterval;	

	/* clear all output definitions */
	memset(outConf, sizeof(outConf), 0);

	if (argc==1) {
		fputs(usage, stderr);
		return 1;
	}

	/* Block the ALARM signal here in order to make sure that any client thread we start (converterThread + tcpserver)
	   will not receive it. Child thread inherit the signal mask of the parent. We will unblock the signal later
	   once the child threads have started.
	   */
	sigemptyset(&signal_mask);
	sigaddset(&signal_mask, SIGALRM);
	rc = pthread_sigmask(SIG_BLOCK, &signal_mask, NULL);
	if (rc!=0) {
		fprintf(stderr, "Cannot change signal mask to block ALARM\n");
		exit(1);
	}

	for (i=1;i<argc;i++) {
		int numFields;
		host_t host;
		int decimate;
		char flags[11], *flag;
		char channelList[30001];

		numFields = sscanf(argv[i], "%255[^:]:%i:%10[^:]:%i:%30000s", host.name, &host.port, flags, &decimate, channelList);
		if (numFields < 5) {
			fprintf(stderr, "Invalid format of output definition (argument) %i -- ignoring\n", i);
			continue;
		}

		if (numOutputs == MAX_OUT) {
			fprintf(stderr, "Cannot define more than %i output buffers - ignoring\n", MAX_OUT);
			continue;
		}

		if (host.port <= 0) {
			fprintf(stderr, "Port number must be positive -- ignoring output definition %i.", i);
			continue;
		}

		if (strlen(host.name)==0) {
			fprintf(stderr, "Hostname cannot be empty -- ignoring output definition %i.", i);
			continue;
		}

		if (decimate<1) {
			fprintf(stderr, "Decimation factor must be >= 1 -- ignoring output definition %i.", i);
			continue;
		}

		/* Spawn tcpserver or connect to remote buffer */
		if (strcmp(host.name, "-") == 0) {
			int rc;

			if (ownBuffer) {
				fprintf(stderr, "Cannot spawn more than one buffer locally -- ignoring output definition %i.", i);
				continue;
			}

			rc = pthread_create(&tcpserverThread, NULL, tcpserver, &host);
			if (rc) {
				fprintf(stderr, "Could not spawn tcpserver thread (%d)\n", rc);
				return 1;
			}
			outConf[numOutputs].ftSocket = 0; /* dma */
			ownBuffer = 1;
		} else {
			int ftSocket = open_connection(host.name, host.port);

			if (ftSocket < 0) {
				fprintf(stderr, "Could not connect to remote buffer on %s:%i\n", host.name, host.port);
				return 1;
			}
			outConf[numOutputs].ftSocket = ftSocket; /* tcp */
		}		

		for (flag = flags; *flag != 0; flag++) {
			switch(*flag) {
				case 'r':
				case 'R':
					outConf[numOutputs].writeRes4 = 1;
					break;
				case 'e':
				case 'E':
					outConf[numOutputs].writeEvents = 1;
					break;
				case 'g':
				case 'G':
					outConf[numOutputs].applyGains = 1;
					break;
				case '-':
					/* silently ignore */
					break;
				default:
					fprintf(stderr, "Warning: Ignoring unknown flag '%c' in output definition %i.\n", *flag, i);
					break;
			}
		}
		outConf[numOutputs].downSample = decimate;
		if (channelList[0]=='*') {
			outConf[numOutputs].channelList = NULL;
			outConf[numOutputs].numChannelsGiven = -1; /* means all channels!!! */
			outConf[numOutputs].numChannelsFound = -1; /* means all channels!!! */
		} else if (channelList[0]==0) {
			outConf[numOutputs].channelList = NULL;
			outConf[numOutputs].numChannelsGiven = 0; /* no channels, possibly events */
			outConf[numOutputs].numChannelsFound = 0;
		} else {
			int n, nChans = 0;
			char *list = strdup(channelList);

			if (list == NULL) {
				fprintf(stderr, "Cannot duplicate channel list -- out of memory\n");
				return 1;
			}

			outConf[numOutputs].channelList = list;

			while (*list) {
				char *cPos = strchr(list, ',');

				outConf[numOutputs].channelNames[nChans++] = list;
				if (cPos == NULL) break; /* end of string */

				/* replace comma by \0 */
				*cPos = 0;
				/* set list pointer to position after that */
				list = cPos+1;
				if (nChans == MAX_CHANNEL) {
					fprintf(stderr, "Cannot handle more than %i channels - skipping rest\n", MAX_CHANNEL);
					break;
				}
			}

			outConf[numOutputs].numChannelsGiven = nChans;
			printf("Selected %i channels:\n", nChans);
			for (n=0;n<nChans-1;n++) {
				printf("'%s', ",outConf[numOutputs].channelNames[n]);
			}
			if (nChans>0) {
				printf("'%s'\n", outConf[numOutputs].channelNames[nChans-1]);
			}
		}
		numOutputs++;
	}

	if (numOutputs == 0) {
		fprintf(stderr, "Warning: No output buffers defined - will only monitor acquisition...\n");
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



	printf("Trying to set up shared memory...\n");
	packet = createSharedMem();
	if (packet == NULL) {
		fprintf(stderr, "Could not set up shared memory - exiting\n");
		exit(1);
	}
	initSharedMem(packet);

	/* Unblock the ALARM signal in this (the main) thread */
	sigemptyset(&signal_mask);
	sigaddset(&signal_mask, SIGALRM);
	rc = pthread_sigmask(SIG_UNBLOCK, &signal_mask, NULL);
	if (rc!=0) {
		fprintf(stderr, "Cannot change signal mask to unblock ALARM in main thread\n");
		exit(1);
	}

	/* register CTRL-C handler */
	signal(SIGINT,  abortHandler);
	/* setup ALARM timer and handler */
	setitimer(ITIMER_REAL, &timerOpts, NULL);
	signal(SIGALRM, alarmHandler);

	/* Both ringbuffers start at 0 */
	currentPacket = intSlot = 0;

	printf("Entering main loop. Press CTRL-C to stop operation.\n\n");
	while (keepRunning) {
		int size;
		char *dataset;
		int code = waitPacket(&packet[currentPacket].message_type, &tNow);

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
				sumT = 0.0;
				sumT2 = 0.0;
				numT = 0;
				break;
			case ACQ_MSGQ_DATA:
				lastId = packet[currentPacket].messageId;
				numChannels  = packet[currentPacket].numChannels;
				numSamples   = packet[currentPacket].numSamples;
				sampleNumber = packet[currentPacket].sampleNumber;

				if (numT > 0) {
					double deltaTM; 
					double deltaT = tNow - tLast;
					sumT += deltaT;
					deltaTM = (deltaT - sumT/numT); 
					sumT2 += deltaTM*deltaTM;
					printf("Data | %3i samples | ID=%4i | slot=%3i | dT=%5.1f  mT=%5.1f  sT=%5.1f\n", 
							numSamples, lastId, currentPacket, deltaT, sumT/numT, sqrt(sumT2/numT));
				} else {
					printf("Data | %3i samples | ID=%4i | slot=%3i\n", 
							numSamples, lastId, currentPacket);
				}
				numT++;
				tLast = tNow;

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
				sumT = 0.0;
				sumT2 = 0.0;
				numT = 0;

				/* Does this mean anything for the buffer as well ? */

				packet[currentPacket].message_type = ACQ_MSGQ_INVALID;
				currentPacket = 0;
				break;
			case ACQ_MSGQ_INVALID:
				printf("Waiting (on %i)\n", currentPacket);
				break;
			default:
				fprintf(stderr, "Warning: Unexpected ID in packet %i (%i)\n", currentPacket, code);
				break;
		}
	}
	printf("Closing sockets / stopping tcpserver...\n");
	close(mySockets[0]); /* the other one will be closed in the thread */
	for (i=0;i<numOutputs;i++) {
		if (outConf[i].ftSocket > 0) {
			close_connection(outConf[i].ftSocket);
		} else {
			pthread_cancel(tcpserverThread);
			pthread_detach(tcpserverThread);
		}
		if (outConf[i].channelList != NULL) {
			free(outConf[i].channelList);
		}
	}
	printf("Closing shared memory...\n");
	closeSharedMem(packet);
	printf("Joining conversion thread...\n");
	pthread_join(convertThread, NULL);
	printf("Done.\n");
	return 0;
}

/** Create the shared memory segment that Acq can write data to.
  In order to overcome a bug in Acq, we overallocate the segment
  by OVERALLOC*sizeof(int), so Acq doesn't segfault when it tries
  to write too much data into the last shared memory slot.
  */
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

/** Clear all slots by setting the message_type to INVALID
*/
void initSharedMem(ACQ_MessagePacketType *packet) {
	int i;
	size_t siz = sizeof(ACQ_MessagePacketType)*ACQ_MSGQ_SIZE;
	siz += OVERALLOC*sizeof(int); /* to overcome Acq bug */

	bzero(packet, siz);
	for (i=0;i<ACQ_MSGQ_SIZE;i++) {
		packet[i].message_type =  ACQ_MSGQ_INVALID;
	}	
}

/** Wait up to 1 second = 100 loops for new packet from Acq.
*/
ACQ_MessageType waitPacket(volatile ACQ_MessageType *msgtyp, double *time) {
	volatile ACQ_MessageType t;
	struct timeval tv;
	int iter = 0;

	while (keepRunning) {
		t = *msgtyp;
		if (t != ACQ_MSGQ_INVALID) break;
		if (++iter >= 100) break;
		sleep(1); /* will actually only sleep max. 10ms due to ALARM signal */
	}
	gettimeofday(&tv, NULL);
	*time = tv.tv_sec * 1000 + tv.tv_usec * 0.001;
	return t;
}

/** Sets 'keepRunning' to 0 in case the uses presses CTRL-C - this makes
  the main thread leave its processing loop, after which the program
  exits smoothly.
  */
void abortHandler(int sig) {
	static int cancelCount = 0;

	printf("Ctrl-C pressed -- stopping acq2ft...\n");
	keepRunning = 0;

	/* In case something goes wrong, exit the hard way */
	if (++cancelCount > 4) exit(1);
}

/** Does not do anything by itself, but having a 100 Hz alarm timer means that
  the sleep() in waitPacket is interrupted after 10ms maximum. This is a
  workaround to cope with the low timer resolution on Odin's 2.4 Linux kernel.
  */
void alarmHandler(int sig) {
}

/** Background thread for grabbing setup and data packets from the internal, overallocated ringbuffer.
  This thread stops automatically if the socket pair is closed in the main thread.
  */
void *dataToFieldTripThread(void *arg) {
	ACQ_OverAllocType *pack;
	EventChain EC, ECd;
	int slot, res, i,j;
	int numChannels = 0, numSamples;
	int warningGiven = 0;

	/*  Initialize event chains. The second one (ECd) is for duplicating
	    the first (EC) such that sample indices can be modified for down-sampling.
	    */
	EC.evs = NULL;
	EC.sizeAlloc = 0;
	ECd.evs = NULL;
	ECd.sizeAlloc = 0;

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
			int i, nChans;
			float fSample;
			ft_chunk_t *chunk;

			chunk = handleRes4((const char *) pack->data, &nChans, &fSample);
			/* clear internal ringbuffer slot */
			pack->message_type = ACQ_MSGQ_INVALID;

			if (chunk == NULL) {
				/* problem while picking up header -- ignore this packet */
				numChannels = 0;
				fprintf(stderr, "SETUP packet received, but could not read .res4 file!\n");
			} else {
				numChannels = nChans;
				matchChannels(numChannels);
				for (i=0;i<numOutputs;i++) {
					writeHeader(&outConf[i], nChans, fSample, chunk);
				}
				warningGiven = 0;
				free(chunk);
			}
		} else if (pack->message_type == ACQ_MSGQ_DATA) {
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
			pack->numChannels = numChannels;

			sampleNumber = pack->sampleNumber;
			numSamples   = pack->numSamples;

			for (i=0;i<numOutputs;i++) {
				if (outConf[i].numChannelsFound != 0) {
					writeSamples(&outConf[i], pack);
				}
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
				for (i=0;i<numOutputs;i++) {
					if (outConf[i].writeEvents) {
						writeEvents(&outConf[i], &EC, &ECd);
					}
				}
			}
		} else {
			fprintf(stderr,"Converter thread: Packet contains neither SETUP nor DATA (%i)...\n", pack->message_type);
		}
		pack->message_type = ACQ_MSGQ_INVALID;
	}
	printf("Leaving converter thread...\n");
	close(mySockets[1]);
	if (EC.sizeAlloc > 0) free(EC.evs);
	if (ECd.sizeAlloc > 0) free(ECd.evs);
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
ft_chunk_t *handleRes4(const char *dsname, int *numChannels, float *fSample) {
	char *res4name;
	ft_chunk_t *chunk;
	int i, len, pdot, pstart, rlen, offset, numFilters;
	FILE *f;
	union {
		double d;
		char b[8];
	} fsamp;
	int nchans;

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
	rlen = sizeof(ft_chunkdef_t) + len;
	chunk = (ft_chunk_t *) malloc(rlen);

	if (chunk == NULL) {
		fprintf(stderr, "Out of memory - can not allocate space for reading .res4 file\n");
		fclose(f);
		return NULL;
	}

	if (fread(chunk->data, 1, len, f) != len) {
		fprintf(stderr, "Could not read complete .res4 file\n");
		fclose(f);
		free(chunk);
		return NULL;
	} 
	fclose(f);			


	/* .res4 file is big-endian, but we assume this machine to be little-endian */
	nchans = getInt16(chunk->data + 1292);
	printf("Number of channels: %i\n", nchans);

	if (numChannels != NULL) *numChannels = nchans;

	for (i=0;i<8;i++) {
		fsamp.b[7-i] = chunk->data[1296+i];
	}
	if (fsamp.d < 1.0) {
		printf("\nWARNING: suspicious sampling frequency (%f Hz) picked from .res4 -- setting to 1200 Hz\n\n", fsamp.d);
		fsamp.d = 1200.0;
	} else {
		printf("Sampling frequency: %f Hz\n", fsamp.d);
	}

	if (fSample != NULL) *fSample = (float) fsamp.d;

	chunk->def.type = FT_CHUNK_CTF_RES4;
	chunk->def.size = len;

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
		len = getInt16(chunk->data +offset + 16);
		offset += 18 + 8*len;
		/* printf("Filter %i: offset = %i,  len = %i\n", i, offset, len); */
	}

	/* next we've got 32 bytes per channel (name) */
	for (i=0;i<nchans;i++) {
		int j;
		char *name_i = chunk->data + offset;

		for (j=0;j<31;j++) {
			if (name_i[j] <= 32 || name_i[j] == '-') break;
		}
		name_i[j] = 0; /* terminate the name with a \0 */
		if (i<MAX_CHANNEL) {
			memcpy(channelNames[i], name_i, j+1); /* copy the string */
			channelNameLen[i] = j;
		}
		offset += 32;
	}

	numTriggerChannels = 0;

	for (i=0;i<nchans && i<MAX_CHANNEL; i++) {
		union {
			char bytes[8];
			double value;
		} sensGain, qGain, ioGain;
		unsigned char ct = chunk->data[offset + 1 + 1328*i]; 

		if (i<MAX_CHANNEL) {
			int j;
			/* turn around big-endian gain values */
			for (j=0;j<8;j++) {
				sensGain.bytes[7-j] = chunk->data[offset + 8+j + 1328*i];
				qGain.bytes[7-j]    = chunk->data[offset + 16+j + 1328*i];
				ioGain.bytes[7-j]   = chunk->data[offset + 24+j + 1328*i];
			}

			channelSensType[i] = ct;
			channelGainV[i] = ioGain.value / (qGain.value * sensGain.value);
			printf("Ch %i: %s  Gain: %g\n", i+1, channelNames[i], channelGainV[i]);
		}
		if (ct == CTF_TRIGGER_TYPE) {
			if (numTriggerChannels < MAX_TRIGGER) {
				triggerChannel[numTriggerChannels] = i;
				triggerChannelNameLen[numTriggerChannels] = strlen(channelNames[i]);
				triggerChannelName[numTriggerChannels] = channelNames[i];
				printf("Trigger channel @ %i: %s\n", i, channelNames[i]);
				++numTriggerChannels;
			}
		}	
	}

	for (i=0;i<numTriggerChannels;i++) {
		lastValue[i] = 0;
	}

	return chunk;
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
	ne->sample      = sample;
	ne->offset      = 0;
	ne->duration    = 0;
	ne->bufsize     = triggerChannelNameLen[trigChan] + sizeof(int);
	memcpy(ntype, triggerChannelName[trigChan], triggerChannelNameLen[trigChan]);
	memcpy(ntype + triggerChannelNameLen[trigChan], &value, sizeof(int));
}

/** Iterates through all output definitions and sets up a table of indices for selected channels.
  Gets called after channel labels are available (handleRes4).
  */
void matchChannels(int numChannels) {
	int i,j,k;

	for (i=0;i<numOutputs;i++) {
		OutputConfig *oc = &outConf[i];
		if (oc->numChannelsGiven <= 0) continue;
		/*
		   printf("%i. output definition:\n", i+1);
		   */
		oc->numChannelsFound = 0;
		for (j=0;j<oc->numChannelsGiven;j++) {
			for (k=0;k<numChannels;k++) {
				if (!strcmp(channelNames[k], oc->channelNames[j])) {
					oc->channels[oc->numChannelsFound++] = k;
					/* 
					   printf("Channel %s has gain %f\n", channelNames[k], channelGainV[k]); 
					   */
				}
			}
		}
	}
}

/** Writes a FieldTrip header to the buffer pointed to by the given output configuration.
  Will transmit the basic header, the channel names chunk, and the res4-file chunk if
  specified by the 'R' flag.
  */
void writeHeader(OutputConfig *oc, int numChannels, float fSample, const ft_chunk_t *res4chunk) {
	int j,res;
	int lenNames = 0, totalLen;
	char *buffer, *bptr;
	int nChans;

	messagedef_t reqdef;
	message_t request, *response;
	headerdef_t *hdef;

	request.def = &reqdef;
	request.buf = NULL;

	/* first, determine the size of the channel name chunk */
	if (oc->numChannelsFound == -1) {
		nChans = numChannels;
		for (j=0;j<nChans;j++) {
			lenNames += channelNameLen[j] + 1;
		}
	} else {
		nChans = oc->numChannelsFound;
		for (j=0;j<nChans;j++) {
			int ch = oc->channels[j];
			lenNames += channelNameLen[ch] + 1;
		}
	}

	totalLen = sizeof(headerdef_t) + sizeof(ft_chunkdef_t) + lenNames;
	if (oc->writeRes4) totalLen += sizeof(ft_chunkdef_t) + res4chunk->def.size;

	buffer = (char *) malloc(totalLen);
	hdef = (headerdef_t *) buffer;

	hdef->nchans  = nChans;
	hdef->fsample = fSample/oc->downSample;
	hdef->data_type = oc->applyGains ? DATATYPE_FLOAT32 : DATATYPE_INT32;
	hdef->nsamples = 0;
	hdef->nevents  = 0;
	hdef->bufsize  = totalLen - sizeof(headerdef_t);

	/* now, fill channel names chunk */
	bptr = buffer + sizeof(headerdef_t);
	((ft_chunkdef_t *) bptr)->type = FT_CHUNK_CHANNEL_NAMES;
	((ft_chunkdef_t *) bptr)->size = lenNames;

	bptr += sizeof(ft_chunkdef_t);
	if (oc->numChannelsFound == -1) {
		for (j=0;j<numChannels;j++) {
			/* copy including trailing zero */
			memcpy(bptr, channelNames[j], channelNameLen[j] + 1);
			bptr+=channelNameLen[j] + 1;
		}
	} else {
		for (j=0;j<nChans;j++) {
			int ch = oc->channels[j];
			memcpy(bptr, channelNames[ch], channelNameLen[ch] + 1);
			bptr+=channelNameLen[ch] + 1;
		}
	}

	/* add RES4 chunk if required */
	if (oc->writeRes4) {
		memcpy(bptr, res4chunk, sizeof(ft_chunkdef_t) + res4chunk->def.size);
	}

	reqdef.version = VERSION;
	reqdef.command = PUT_HDR;
	reqdef.bufsize = totalLen;
	request.buf = buffer;

	oc->headerOk = 0;
	oc->skipSamples = 0;

	res = clientrequest(oc->ftSocket, &request, &response);

	if (res < 0) {
		fprintf(stderr, "Error in FieldTrip connection\n");
	} else if (response) {
		if (response->def->command != PUT_OK) {
			fprintf(stderr, "Error in PUT_HDR\n");
		} else {
			/* printf("FT: Transmitted header\n"); */
			oc->headerOk = 1;
		}
		cleanup_message((void **) &response);
	}
	free(buffer);
}

/** Writes samples to the buffer pointed to by the given output configuration.
  Will downsample (without filtering!) if specified by the decimation value.
  */
void writeSamples(OutputConfig *oc, const ACQ_OverAllocType *pack) {
	int nsamp, nchans;
	messagedef_t reqdef;
	message_t request, *response;
	int i,j, res;

	request.def = &reqdef;
	request.buf = &putDatBuf;

	nsamp = 0;

	if (oc->applyGains) {
		/* write samples to floating point data buffer pointed to by 'dest' */
		float *dest = putDatBuf.fData;
		putDatBuf.ddef.data_type = DATATYPE_FLOAT32;

		if (oc->numChannelsFound == -1) {
			/* transmit all channels, but with gains applied */
			nchans = pack->numChannels;
			for (j=oc->skipSamples; j<pack->numSamples; j+=oc->downSample) {
				const int *source_j = pack->data + j*pack->numChannels;
				for (i=0;i<nchans;i++) {
					*dest++ = channelGainV[i]*source_j[i];
				}
				++nsamp;
			}
		} else {
			/* only transmit selected channels with gains applied */
			nchans = oc->numChannelsFound;
			for (j=oc->skipSamples; j<pack->numSamples; j+=oc->downSample) {
				const int *source_j = pack->data + j*pack->numChannels;
				for (i=0;i<nchans;i++) {
					int ch = oc->channels[i];
					*dest++ = channelGainV[ch]*source_j[ch];
				}
				++nsamp;
			}
		}

	} else {
		/* write samples to integer data buffer pointed to by 'dest' */
		int *dest = putDatBuf.data;
		putDatBuf.ddef.data_type = DATATYPE_INT32;

		if (oc->numChannelsFound == -1) {
			/* transmit all channels without gains */
			nchans = pack->numChannels;
			for (j=oc->skipSamples; j<pack->numSamples; j+=oc->downSample) {
				const int *source_j = pack->data + j*pack->numChannels;
				for (i=0;i<nchans;i++) {
					*dest++ = source_j[i];
				}
				++nsamp;
			}
		} else {
			/* only transmit selected channels without gains */
			nchans = oc->numChannelsFound;
			for (j=oc->skipSamples; j<pack->numSamples; j+=oc->downSample) {
				const int *source_j = pack->data + j*pack->numChannels;
				for (i=0;i<nchans;i++) {
					*dest++ = source_j[oc->channels[i]];
				}
				++nsamp;
			}
		}
	}
	oc->skipSamples = j - pack->numSamples;
	putDatBuf.ddef.nchans   = nchans;
	putDatBuf.ddef.nsamples = nsamp;
	putDatBuf.ddef.bufsize  = 4 * nsamp * nchans; /* both int32 + float32 are 4 bytes wide */

	reqdef.version = VERSION;
	reqdef.command = PUT_DAT;
	reqdef.bufsize = putDatBuf.ddef.bufsize + sizeof(datadef_t);
	request.buf    = &putDatBuf.ddef; 

	res = clientrequest(oc->ftSocket, &request, &response);

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
}

/** Writes events to the buffer pointed to by the given output configuration.
  Will correct sample indices according to the specified decimation value.
  The logic here is that first the EventChain (list of events) will be
  duplicated into ECd (which should be pre-allocated most of the time).
  */
void writeEvents(OutputConfig *oc, const EventChain *EC, EventChain *ECd) {
	messagedef_t reqdef;
	message_t request, *response;
	int res;

	request.def = &reqdef;
	request.buf = &putDatBuf;

	reqdef.version = VERSION;
	reqdef.command = PUT_EVT;
	reqdef.bufsize = EC->size;

	if (oc->downSample == 1) {
		/* no need to translate */
		request.buf = EC->evs;
	} else {
		int i, offset = 0;

		/* first see if ECd is big enough */
		if (ECd->sizeAlloc < EC->size) {
			if (ECd->evs != NULL) {
				free(ECd->evs);
			}
			ECd->evs = malloc(EC->size);
			if (ECd->evs == NULL) {
				ECd->sizeAlloc = 0;
				fprintf(stderr, "Cannot duplicate events for downsampling -- out of memory\n");
				return;
			}
			ECd->sizeAlloc = EC->size;
		}
		memcpy(ECd->evs, EC->evs, EC->size);

		for (i=0;i<EC->num;i++) {
			eventdef_t *evdef = (eventdef_t *) ((char *) ECd->evs + offset);
			evdef->sample /= oc->downSample;
			offset += evdef->bufsize + sizeof(eventdef_t);
		}
		request.buf = ECd->evs;
	}

	res = clientrequest(oc->ftSocket, &request, &response);

	if (res < 0) {
		fprintf(stderr, "Error in FieldTrip connection (writing events)\n");
	} else if (response) {
		if (response->def->command != PUT_OK) {
			fprintf(stderr, "Error in PUT_EVT\n");
		} else {
			printf("Wrote %i events (%i bytes)\n", EC->num, reqdef.bufsize);
			/* printf("FT: Transmitted samples\n"); */
		}
	}
	cleanup_message((void **) &response);
}

