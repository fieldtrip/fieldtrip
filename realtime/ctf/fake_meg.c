#include <stdio.h>
#include "AcqBuffer.h"
#include <sys/shm.h>
#include <sys/types.h>
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

  
int main(int argc, char **argv) {
	ACQ_MessagePacketType *packet;
	int currentPacket = 0;
	int numChannels = 356;
	int blockSize = 60;  /* this is too much for regular operation! */
	int numPackets = 42;
	int sampleNumber = 0;
	int i;
		
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
	strcpy((char *)(packet[0].data), "/home/mrphys/stekla/meg/buffertest_1200hz_20100503_01.ds");
	packet[0].message_type = ACQ_MSGQ_SETUP_COLLECTION;

	currentPacket = 1;
	
	printf("Wrote setup packet at 0\n");
	
	usleep(200000); /* 0.2 sec */
	
	for (i=0;i<numPackets;i++) {
      
		if (packet[currentPacket].message_type != ACQ_MSGQ_INVALID) {
			printf("Packet #%i not free...\n", currentPacket);
			break;
		}
	
		packet[currentPacket].messageId = i+1;
		packet[currentPacket].sampleNumber = sampleNumber;
		packet[currentPacket].numSamples = blockSize;
		packet[currentPacket].numChannels = numChannels;
		memset(packet[currentPacket].data, 0, numChannels * blockSize * sizeof(int));
		packet[currentPacket].data[0] = i+1;
		packet[currentPacket].data[10*numChannels] = 10;
		packet[currentPacket].numChannels = numChannels;
		packet[currentPacket].message_type = ACQ_MSGQ_DATA;
		
		printf("Wrote data packet %i at %i\n", i, currentPacket);
		
		sampleNumber += blockSize;
		
		if (++currentPacket == ACQ_MSGQ_SIZE) currentPacket = 0;
		
		usleep(50000);
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

