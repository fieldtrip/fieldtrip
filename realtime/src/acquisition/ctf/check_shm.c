/* 
 * Copyright (C) 2010, Stefan Klanke
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 * Small program to inspect the CTF shared memory segment, and to delete it if desired.
 */
#include <stdio.h>
#include "AcqBuffer.h"
#include <sys/shm.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <unistd.h>

#define OVERALLOC   0    /* to overcome Acq bug */

/*	 uncomment this for testing on different machines with fake_meg
#undef  ACQ_MSGQ_SIZE   
#define ACQ_MSGQ_SIZE 10
#undef  ACQ_MSGQ_SHMKEY 
#define ACQ_MSGQ_SHMKEY    0x08150842 
*/

int main(int argc, char **argv) {
	void *ptr;
	int shmid, r;
	size_t siz = sizeof(ACQ_MessagePacketType)*ACQ_MSGQ_SIZE;
	struct shmid_ds DS;
   
    siz += OVERALLOC*sizeof(int); /* to overcome Acq bug */
	
	
	siz = 10;
	
	printf("Trying to access shared memory, asking for %i bytes\n", siz);
		
	shmid = shmget(ACQ_MSGQ_SHMKEY, 10, 0666|IPC_CREAT);
	if (shmid == -1) {
		perror("shmget");
		return 1;
	}
	
	r = shmctl(shmid, IPC_STAT, &DS);
	if (r==0) {
		printf("SHM info\n");
		printf("size (bytes)     = %i\n", DS.shm_segsz);
		printf("current attaches = %i\n", DS.shm_nattch);
		printf("PID of creator   = %i\n", DS.shm_cpid);
		printf("PID of last opr. = %i\n", DS.shm_lpid);
		printf("UID of owner     = %i\n", DS.shm_perm.uid);
		printf("UID of creator   = %i\n", DS.shm_perm.cuid);
	}

	printf("Trying to attach...\n");

	/* attach to the segment to get a pointer to it */
	ptr = shmat(shmid, NULL, 0);
	if (ptr == (void *) -1) {
		perror("shmat");
		return 1;
	}
	
	r = shmctl(shmid, IPC_STAT, &DS);
	if (r==0) {
		printf("Updated SHM info\n");
		printf("current attaches = %i\n", DS.shm_nattch);
		printf("PID of creator   = %i\n", DS.shm_cpid);
		printf("PID of last opr. = %i\n", DS.shm_lpid);
	}
	
	printf("Detaching...\n");
	shmdt(ptr);
	
	printf("Do you want to destroy the SHM (y/n) ?\n");
	r = getchar();
	if (r=='y') {
		r = shmctl(shmid, IPC_RMID, &DS);
		if (r) {
			perror("shmctl");
		} else {
			printf("Ok!\n");
		}
	}
	
	return 0;
}
