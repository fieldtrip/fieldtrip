/* 
 * Copyright (C) 2010, Stefan Klanke
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 * This is just a small test program for writing (too big) packets to a
 * shared memory segment similar to how CTF's Acq operates.
 */
#include <stdio.h>
#include "AcqBuffer.h"
#include <sys/shm.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/ipc.h>
#include <unistd.h>
#include <signal.h>
#include <stdlib.h>
#include <string.h>

#undef  ACQ_MSGQ_SIZE   
#define ACQ_MSGQ_SIZE 10

#undef  ACQ_MSGQ_SHMKEY 
#define ACQ_MSGQ_SHMKEY    0x08150842

/* prototypes for helper functions defined below */
ACQ_MessagePacketType *createSharedMem();
void closeSharedMem(ACQ_MessagePacketType *packet);
void initSharedMem(ACQ_MessagePacketType *packet);
volatile int keepRunning = 1;

void alarmHandler(int signal) {
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


int main(int argc, char **argv) {
	ACQ_MessagePacketType *packet;
	int currentPacket = 0;
	int numChannels = 356;
	int blockSize = 70;  /* > 28160/356 is too much for regular operation! */
	int sampleNumber = 0;
	int ID = 0;
	int sampleFreq = 1200;
	struct timeval timerInterval = {0,50000}; /*  50ms */
	struct timeval timerValue    = {0,50000};  /* 200ms */
	struct itimerval timerOpts;
	
	timerOpts.it_value = timerValue;
	timerOpts.it_interval = timerInterval;
	
	if (argc>=2) {
		sampleFreq = atoi(argv[1]);
		if (sampleFreq>0) {
			timerOpts.it_interval.tv_usec = timerOpts.it_value.tv_usec = 1000000*blockSize/sampleFreq;
		}
	}
		
	printf("Trying to set up shared memory...\n");
	packet = createSharedMem();
	if (packet == NULL) {
		fprintf(stderr, "Could not set up shared memory - exiting\n");
		exit(1);
	}
	initSharedMem(packet);
	
	packet[0].messageId    = 0;
	packet[0].sampleNumber = 0;
	packet[0].numSamples   = 0;
	packet[0].numChannels  = numChannels;
	strcpy((char *)(packet[0].data), "test.ds");
	packet[0].message_type = ACQ_MSGQ_SETUP_COLLECTION;

	currentPacket = 1;
	
	printf("Wrote setup packet at 0\n");
	
	
	usleep(200000);  /* 0.2 sec */
	setitimer(ITIMER_REAL, &timerOpts, NULL);
	signal(SIGALRM, alarmHandler);
	/* register CTRL-C handler */
	signal(SIGINT,  abortHandler);

	
	while (keepRunning) {
		++ID;
      
		if (packet[currentPacket].message_type != ACQ_MSGQ_INVALID) {
			printf("Packet #%i not free...\n", currentPacket);
			break;
		}
	
		packet[currentPacket].messageId = ID;
		packet[currentPacket].sampleNumber = sampleNumber;
		packet[currentPacket].numSamples = blockSize;
		packet[currentPacket].numChannels = numChannels;
		memset(packet[currentPacket].data, 0, numChannels * blockSize * sizeof(int));
		//packet[currentPacket].data[0] = i+1;
		//packet[currentPacket].data[10*numChannels] = 10;
		packet[currentPacket].numChannels = numChannels;
		packet[currentPacket].message_type = ACQ_MSGQ_DATA;
		
		printf("Wrote data packet %i at %i\n", ID, currentPacket);
		
		sampleNumber += blockSize;
		
		if (++currentPacket == ACQ_MSGQ_SIZE) currentPacket = 0;
		
		sleep(1);
	}
	
	packet[currentPacket].message_type = ACQ_MSGQ_CLOSE_CONNECTION;
	usleep(100000);

	closeSharedMem(packet);
	return 0;
}


ACQ_MessagePacketType *createSharedMem() {
	void *ptr;
	int shmid;
	size_t siz = sizeof(ACQ_MessagePacketType)*ACQ_MSGQ_SIZE;
		
	shmid = shmget(ACQ_MSGQ_SHMKEY, siz, 0666);
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

