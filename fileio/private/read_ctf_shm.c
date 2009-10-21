/*
 * READ_CTF_SHM reads metainformation or selected blocks of data from
 * shared memory. This function can be used for real-time processing of
 * data while it is being acquired.
 *
 * Use as
 *   [msgType msgId sampleNumber numSamples numChannels] = read_ctf_shm;
 * or
 *   [data] = read_ctf_shm(msgNumber);
 *   [data] = read_ctf_shm(msgNumber, numValues);
 *
 * Copyright (C) 2007, Robert Oostenveld
 *
 * $Log: read_ctf_shm.c,v $
 * Revision 1.1  2009/01/14 09:43:37  roboos
 * moved source code for mex file from fileio/mex to file/private
 * compiling the source code from within Matlab now ensures that the mex file will be located immediately at the right position
 *
 * Revision 1.3  2007/08/01 09:36:16  roboos
 * changed some whitespace
 *
 * Revision 1.2  2007/07/27 12:35:24  roboos
 * implemented reading a selected number of samples (default is all)
 *
 * Revision 1.1  2007/07/24 11:17:07  roboos
 * first implementation, tested and works on odin
 *
 */

#include <string.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include "mex.h"
#include "matrix.h"

#define ACQ_MSGQ_SIZE      600
#define ACQ_MSGQ_SHMKEY    0x39457f73
#define ACQ_MSGQ_SHMPROJID 12345
#define ACQ_MSGQ_SHMPATH   "/opt/ctf/bin/Acq"

typedef enum
{
  ACQ_MSGQ_SETUP_COLLECTION,
  ACQ_MSGQ_DATA,
  ACQ_MSGQ_CLOSE_CONNECTION,
  ACQ_MSGQ_INVALID = INT_MAX
} ACQ_MessageType;

typedef struct
{
  ACQ_MessageType message_type;
  int messageId;
  int sampleNumber;
  int numSamples;
  int numChannels;
  int data[28160];
} ACQ_MessagePacketType;

void
mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
  key_t key;
  ACQ_MessagePacketType *packet;
  int *msgType, *msgId, *sampleNumber, *numSamples, *numChannels;
  int shmid, shmsize, i;
  int numValues = 28160;

  shmsize = sizeof(ACQ_MessagePacketType) * ACQ_MSGQ_SIZE;

  /* make the key
  if ((key = ftok(ACQ_MSGQ_SHMPATH, ACQ_MSGQ_SHMPROJID)) == -1) {
    perror("ftok");
    exit(1);
    }
   */

  /* use the pre-defined key */
  key = ACQ_MSGQ_SHMKEY;

  /* connect to (and possibly create) the segment */
  if ((shmid = shmget(key, shmsize, 0644 | IPC_CREAT )) == -1)
    mexErrMsgTxt("shmget");

  /* attach to the segment to get a pointer to it */
  packet = shmat(shmid, (void *)0, 0);
  if ((char *)packet == (char *)(-1))
    mexErrMsgTxt("shmat");

  if (nrhs==0) {
    /* read the meta information from all packets */
    plhs[0] = mxCreateNumericMatrix(1, ACQ_MSGQ_SIZE, mxINT32_CLASS, mxREAL);
    plhs[1] = mxCreateNumericMatrix(1, ACQ_MSGQ_SIZE, mxINT32_CLASS, mxREAL);
    plhs[2] = mxCreateNumericMatrix(1, ACQ_MSGQ_SIZE, mxINT32_CLASS, mxREAL);
    plhs[3] = mxCreateNumericMatrix(1, ACQ_MSGQ_SIZE, mxINT32_CLASS, mxREAL);
    plhs[4] = mxCreateNumericMatrix(1, ACQ_MSGQ_SIZE, mxINT32_CLASS, mxREAL);

    msgType      = mxGetData(plhs[0]);
    msgId        = mxGetData(plhs[1]);
    sampleNumber = mxGetData(plhs[2]);
    numSamples   = mxGetData(plhs[3]);
    numChannels  = mxGetData(plhs[4]);

    for (i=0; i<ACQ_MSGQ_SIZE; i++) {
      msgType[i]      = (int)(packet[i].message_type);
      msgId[i]        = packet[i].messageId;
      sampleNumber[i] = packet[i].sampleNumber;
      numSamples[i]   = packet[i].numSamples;
      numChannels[i]  = packet[i].numChannels;
    }
  }
  else {
    if (nrhs==2) {
      numValues = (int)mxGetScalar(prhs[1]);
      numValues = ( numValues>28160 ? 28160 : numValues );  /* check boundary */
      numValues = ( numValues<0     ? 0     : numValues );  /* check boundary */
    }

    /* read the data from the selected packet */
    /* one-offset in Matlab, zero-offset in C */
    i = mxGetScalar(prhs[0]) - 1;
    if (i<0)
      mexErrMsgTxt("Cannot read before the first packet");
    if (i>=ACQ_MSGQ_SIZE)
      mexErrMsgTxt("Cannot read after the last packet");
    plhs[0] = mxCreateNumericMatrix(1, numValues, mxINT32_CLASS, mxREAL);
    memcpy(mxGetData(plhs[0]), packet[i].data, numValues*sizeof(int));
  }

  /* detach from the segment */
  if (shmdt(packet) == -1)
    mexErrMsgTxt("shmdt");

} /* end of mexFunction */

