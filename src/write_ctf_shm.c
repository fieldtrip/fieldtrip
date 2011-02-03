/*
 * WRITE_CTF_SHM writes a packet to shared memory. This function can be
 * used for real-time processing of data while it is being acquired.
 *
 * Use as
 *   write_ctf_shm(msgNumber, msgType, msgId, sampleNumber, numSamples, numChannels, data);
 *
 * Copyright (C) 2007, Robert Oostenveld
 *
 * $Log: write_ctf_shm.c,v $
 * Revision 1.1  2009/01/14 09:43:37  roboos
 * moved source code for mex file from fileio/mex to file/private
 * compiling the source code from within Matlab now ensures that the mex file will be located immediately at the right position
 *
 * Revision 1.2  2007/08/01 09:36:16  roboos
 * changed some whitespace
 *
 * Revision 1.1  2007/07/30 13:36:04  roboos
 * new implementation
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
  int msgNumber, msgType, msgId, sampleNumber, numSamples, numChannels, *data;
  int shmid, shmsize;
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

  if (nrhs<7)
    mexErrMsgTxt("Not enough input arguments");

  msgNumber    = (int)mxGetScalar(prhs[0])-1; /* one offset in Matlab, zero offset in C */
  /*
  msgType      = (int)mxGetScalar(prhs[1]);
  msgId        = (int)mxGetScalar(prhs[2]);
  sampleNumber = (int)mxGetScalar(prhs[3]);
  numSamples   = (int)mxGetScalar(prhs[4]);
  numChannels  = (int)mxGetScalar(prhs[5]);
  */

  if (mxGetClassID(prhs[6]) != mxINT32_CLASS)
    mexErrMsgTxt("Invalid type of data, should be int32");

  if (msgNumber<0)
    mexErrMsgTxt("Cannot write before the first packet");
  if (msgNumber>=ACQ_MSGQ_SIZE)
    mexErrMsgTxt("Cannot write after the last packet");

  numValues = mxGetNumberOfElements(prhs[6]);
  numValues = ( numValues>28160 ? 28160 : numValues );  /* check boundary */
  numValues = ( numValues<0     ? 0     : numValues );  /* check boundary */

   /* write the meta-information to the packet */
  packet[msgNumber].message_type = (int)mxGetScalar(prhs[1]);
  packet[msgNumber].messageId    = (int)mxGetScalar(prhs[2]);
  packet[msgNumber].sampleNumber = (int)mxGetScalar(prhs[3]);
  packet[msgNumber].numSamples   = (int)mxGetScalar(prhs[4]);
  packet[msgNumber].numChannels  = (int)mxGetScalar(prhs[5]);

   /* write the data to the packet */
  memcpy(packet[msgNumber].data, mxGetData(prhs[6]), numValues*sizeof(int));

  /* detach from the segment */
  if (shmdt(packet) == -1)
    mexErrMsgTxt("shmdt");

} /* end of mexFunction */

