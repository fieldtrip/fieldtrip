#include <stdio.h>
#include "AcqBuffer.h"
#include <buffer.h>
#include <sys/shm.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <unistd.h>
#include <signal.h>
#include <pthread.h>


#define OVERALLOC  10000   /* to overcome Acq bug */

/*	 uncomment this for testing on different machines with fake_meg
*/
#undef  ACQ_MSGQ_SIZE   
#define ACQ_MSGQ_SIZE 10

#undef  ACQ_MSGQ_SHMKEY 
#define ACQ_MSGQ_SHMKEY    0x08150816
/*
*/

static char usage[] = 
"Usage: acq2ft [hostname port]\n" \
"Just calling 'acq2ft' without parameters starts a local FieldTrip buffer on port 1972.\n" \
"Using '-' for the hostname starts a local buffer on the given port.\n";

/* this is only used internally to retrieve data from the socket pair */
typedef struct {
  ACQ_MessageType message_type;
  int messageId;
  int sampleNumber;
  int numSamples;
  int numChannels;
  int data[28160 + OVERALLOC];
} ACQ_OverAllocType;

volatile int keepRunning = 1;
int ftSocket = -1;
int mySockets[2];	/* keeps a socket pair */

/* prototypes for helper functions defined below */
ACQ_MessagePacketType *createSharedMem();
void closeSharedMem(ACQ_MessagePacketType *packet);
void initSharedMem(ACQ_MessagePacketType *packet);
int waitPacket(volatile ACQ_MessagePacketType *pi);
void abortHandler(int sig);
void *dataToFieldTripThread(void *arg);

int main(int argc, char **argv) {
	ACQ_MessagePacketType *packet;
	int currentPacket = 0;
	int lastId, optval;
	int numChannels, numSamples, sampleNumber;
	char dataset[1024] = {0};
	host_t host;
	pthread_t tcpserverThread, convertThread;
	
	/* Parse command line arguments */
	if (argc == 1) {
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
   
   /* set large send buffer for mySockets[0] */
   optval = 4*sizeof(ACQ_MessagePacketType);
   if (setsockopt(mySockets[0], SOL_SOCKET, SO_SNDBUF, (const char*)&optval, sizeof(optval)) < 0) {
      perror("setsockopt/SO_SNDBUF");
	}
   if (setsockopt(mySockets[1], SOL_SOCKET, SO_RCVBUF, (const char*)&optval, sizeof(optval)) < 0) {
      perror("setsockopt/SO_RCVBUF");
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
	
	printf("Entering main loop. Press CTRL-C to stop operation.\n\n");
	while (keepRunning) {
		int size;
		int code = waitPacket(&packet[currentPacket]);
		switch(code) {
			case ACQ_MSGQ_SETUP_COLLECTION:
				lastId = packet[currentPacket].messageId;
				numChannels = packet[currentPacket].numChannels;
				strncpy(dataset, (char *) packet[currentPacket].data, 1024);
				dataset[1023] = 0;
				printf("Setup | ID=%i | %i channels | %s\n", lastId, numChannels, dataset);
				
				size = 5*sizeof(int) + strlen(dataset) + 1;
				bufwrite(mySockets[0], &size, sizeof(size));
				bufwrite(mySockets[0], &packet[currentPacket], size);
				
				packet[currentPacket].message_type = ACQ_MSGQ_INVALID;
				
				/* TODO: grab the data here and put it into FieldTrip */
				
				if (++currentPacket == ACQ_MSGQ_SIZE) currentPacket=0;
				break;
			case ACQ_MSGQ_DATA:
				lastId = packet[currentPacket].messageId;
				numChannels  = packet[currentPacket].numChannels;
				numSamples   = packet[currentPacket].numSamples;
				sampleNumber = packet[currentPacket].sampleNumber;
				printf("Data | %3i channels x %3i samples | nr = %6i | ID=%i | slot=%i\n", numChannels, numSamples, sampleNumber, lastId, currentPacket);
				
				size = 5*sizeof(int) + numSamples*numChannels*sizeof(int);
				bufwrite(mySockets[0], &size, sizeof(size));
				bufwrite(mySockets[0], &packet[currentPacket], size);
				
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
int waitPacket(volatile ACQ_MessagePacketType *pi) {
	int i,t;
	for (i=0;i<1000 && keepRunning;i++) {
		t = pi->message_type;
		if (t != ACQ_MSGQ_INVALID) break;
		usleep(1000);
	}
	return t;
}

void abortHandler(int sig) {
	printf("Ctrl-C pressed -- stopping acq2ft...\n");
	keepRunning = 0;
}


void *dataToFieldTripThread(void *arg) {
	ACQ_OverAllocType packet;
	int size,sr;
	int numChannels = 0, numSamples;
	messagedef_t reqdef;
	message_t request, *response;
	
	request.def = &reqdef;
	
	while (1) {
		sr = bufread(mySockets[1], &size, sizeof(size));
		if (sr <= 0) break;
		
		if (size > sizeof(packet)) {
			fprintf(stderr, "Packet too large to fit\n");
			exit(2);
		}
		
		sr = bufread(mySockets[1], &packet, size);
		/* read on a socket pait should always give desired size */
		if (sr != size) {
			fprintf(stderr, "Unexpected read error (%i vs. %i)\n", size, sr);
			exit(2);
		}
		
		if (packet.message_type == ACQ_MSGQ_SETUP_COLLECTION) {
			char *dsname = (char *) packet.data;
			char *res4name, *aux;
			headerdef_t *hdef;
			ft_chunk_t *chunk;
			int len, pdot, pstart, rlen;
			FILE *f;
			
			/* printf("Picked up dataset name : %s\n", dsname); */
			len = strlen(dsname);
			
			pdot = len;
			while (--pdot>0) {
				if (dsname[pdot] == '.') break;
			}
			if (pdot == 0) {
				fprintf(stderr, "No . found in dataset name - don't know where to pick up .res4 file\n");
				continue;
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
				continue;
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
				continue;
			}
			
			fseek(f, 0, SEEK_END);
			len = ftell(f);
			fseek(f, 0, SEEK_SET);
			
			printf("CTF RES4 file %s contains %i bytes.\n", res4name, len);
			free(res4name); /* not needed anymore */
			
			rlen = len + sizeof(headerdef_t) + sizeof(ft_chunkdef_t);
			aux = (char *) malloc(rlen);
			if (aux == NULL) {
				fprintf(stderr, "Out of memory - can not allocate space for reading .res4 file\n");
				fclose(f);
				continue;
			}
			
			hdef = (headerdef_t *) aux;
			chunk = (ft_chunk_t *) (aux + sizeof(headerdef_t));
			
			if (fread(chunk->data, 1, len, f) != len) {
				fprintf(stderr, "Could not read complete .res4 file\n");
				fclose(f);
				free(aux);
				continue;
			}
			fclose(f);
			
			hdef->nchans = packet.numChannels;
			hdef->nsamples = 0;
			hdef->nevents = 0;
			hdef->fsample = 1200.0; /* TODO: check this */
			hdef->data_type = DATATYPE_INT32;
			hdef->bufsize = sizeof(ft_chunkdef_t) + len;
			chunk->def.type = FT_CHUNK_CTF_RES4;
			chunk->def.size = len;
			
			reqdef.version = VERSION;
			reqdef.command = PUT_HDR;
			reqdef.bufsize = sizeof(headerdef_t) + rlen;
			request.buf = aux;
			
			sr = clientrequest(ftSocket, &request, &response);
			
			free(aux);
			
			if (sr < 0) {
				fprintf(stderr, "Error in FieldTrip connection\n");
			} else if (response) {
				if (response->def->command != PUT_OK) {
					fprintf(stderr, "Error in PUT_HDR\n");
				} else {
					/* printf("FT: Transmitted header\n"); */
					numChannels = packet.numChannels;
				}
				cleanup_message((void **) &response);
			}
		} else if (packet.message_type == ACQ_MSGQ_DATA) {
			datadef_t *ddef = (datadef_t *) &packet.messageId; /* This just fits ! */
			
			if (numChannels == 0) {
				fprintf(stderr, "No header written yet -- ignoring data packet\n");
				continue;
			}
			if (numChannels != packet.numChannels) {
				fprintf(stderr, "Number of channels in data packet (%i) does not match setup (%i)\n", packet.numChannels, numChannels);
				continue;
			}
			
			numChannels = packet.numChannels;
			numSamples = packet.numSamples;
			
			ddef->nsamples = numSamples;
			ddef->nchans = numChannels;
			ddef->data_type = DATATYPE_INT32;
			ddef->bufsize = sizeof(int) * numSamples * numChannels;
			
			reqdef.version = VERSION;
			reqdef.command = PUT_DAT;
			reqdef.bufsize = ddef->bufsize + sizeof(datadef_t);
			request.buf = ddef; /* data is still behind that */
			
			sr = clientrequest(ftSocket, &request, &response);
			
			if (sr < 0) {
				fprintf(stderr, "Error in FieldTrip connection\n");
			} else if (response) {
				if (response->def->command != PUT_OK) {
					fprintf(stderr, "Error in PUT_DAT\n");
				} else {
					/* printf("FT: Transmitted samples\n"); */
				}
				cleanup_message((void **) &response);
			}
		} else {
			fprintf(stderr,"Converter thread: Packet contains neither SETUP nor DATA (%i)...\n", packet.message_type);
		}
	}
	printf("Leaving converter thread...\n");
	close(mySockets[1]);
	return NULL;
}
