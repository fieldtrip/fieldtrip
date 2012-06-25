/*
 * AcqBuffer.c
 *
 * Copyright (C) 2007, Robert Oostenveld
 *
 * $Id$
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <ncurses.h>
#include <unistd.h>  /* only for usleep */
#include "AcqBuffer.h"

/*#define ACQ_MSGQ_SIZE 20*/
#define SHOWPACKET 20
#define MINFREE 30                  /* number of records to keep free, 30 seems to be enough up to 12kHz */
#define PERMISSION 0666             /* for the shared memory */
#define PAUSE 32000                 /* in micro seconds */
#define DATASIZE 28160              /* number of data samples in each packet */
#define TRIGSIZE 3000               /* number of triggers that is maintained */
#define MAXNUMTRIGCHAN 10           /* maximum number of trigger channels that can be maintained */

#define limitrange(n, low, high) (n<low ? low : (n>high ? high : n))
#define wrapnumpacket(n)   (n<0 ? n+ACQ_MSGQ_SIZE : (n>=ACQ_MSGQ_SIZE ? 0 : n))
#define wrapnumtrigger(n)  (n<0 ? n+(3*TRIGSIZE)  : (n>=(3*TRIGSIZE) ? 0 : n))
#define display_print      if (showDisplay) printw

ACQ_MessagePacketType *packet;

void write_setup(int);
void write_data(int, int);
void clear_packet(int);
void clear_all(void);
void clear_all_except_setup(void);
void print_packet(int);
void display_start(void);
void display_stop(void);
void display_refresh(void);

int numRealChan = 0, numTrigChan = 0, numTrigger = 0, *trigChan, *trigPointer, lastValue[MAXNUMTRIGCHAN], lastPacket = -1;
int showDisplay = 1;

int main(int argc, char *argv[]) {
  key_t shmkey;
  int shmid, shmsize;
  int i, j, k, t = 0, c = 0;
  int pause = PAUSE, again = 1;
  int setupIndex = -1;
  int minSampleIndex = -1, minMessageIndex = -1;
  int maxSampleIndex = -1, maxMessageIndex = -1;
  int minSampleValue = INT_MAX, minMessageValue = INT_MAX;
  int maxSampleValue = INT_MIN, maxMessageValue = INT_MIN;
  int maxCountCleared = 0, countCleared = 0, countLoop = 0;
  int countInvalid = 0, countData = 0, countSetup = 0, countCancel = 0, countUnknown = 0;

  for (i=0; i<MAXNUMTRIGCHAN; i++)
    lastValue[i] = 0;

  /* make the shared memory key
   * if ((shmkey = ftok(ACQ_MSGQ_SHMPATH, ACQ_MSGQ_SHMPROJID)) == -1) {
   *   perror("ftok");
   *   exit(1);
   * }
   */

  /* use the pre-defined shared memory key */
  shmkey = ACQ_MSGQ_SHMKEY;

  /* determine the size of the shared memory buffer */
  shmsize = sizeof(ACQ_MessagePacketType) * ACQ_MSGQ_SIZE;

  /* connect to (and possibly create) the segment */
  if ((shmid = shmget(shmkey, shmsize, PERMISSION | IPC_CREAT )) == -1) {
    perror("shmget");
    exit(1);
  }

  /* attach to the segment to get a pointer to it */
  packet = shmat(shmid, (void *)0, 0);
  if ((char *)packet == (char *)(-1)) {
    perror("shmat");
    exit(1);
  }

  display_start();
  display_print("initializing\n");
  /* clear_all_except_setup();
   clear_all();*/

  display_refresh();

  while (again) {
    clear();

    setupIndex = -1;
    minSampleIndex = -1, minMessageIndex = -1;
    maxSampleIndex = -1, maxMessageIndex = -1;
    minSampleValue = INT_MAX, minMessageValue = INT_MAX;
    maxSampleValue = INT_MIN, maxMessageValue = INT_MIN;
    countInvalid = 0, countData = 0, countSetup = 0, countCancel = 0, countUnknown = 0;

    /*******************************************************************************
     * collect information about all packets in the buffer
     *******************************************************************************/

    for (i=0; i<ACQ_MSGQ_SIZE; i++) {
      if (i<SHOWPACKET)
        print_packet(i);

      switch (packet[i].message_type) {
        case ACQ_MSGQ_SETUP_COLLECTION:
          if (countSetup>0) {
            /* multiple packets with a setup are not allowed, only keep the first */
            display_print("clearing superfluous setup packet at %d (setup @ %d)\n", i, setupIndex);
            clear_packet(i);
            countCleared++;
            break;
          }
          countSetup++;
          setupIndex = i;
          /* update the specifications that relate to trigger maintenance */
          numRealChan = packet[i].data[DATASIZE-TRIGSIZE*3-MAXNUMTRIGCHAN-2];  /* update the value */
          numTrigChan = packet[i].data[DATASIZE-TRIGSIZE*3-MAXNUMTRIGCHAN-1];  /* update the value */
          trigChan    = packet[i].data+DATASIZE-TRIGSIZE*3-MAXNUMTRIGCHAN;     /* update the pointer */
          trigPointer = packet[i].data+DATASIZE-TRIGSIZE*3;                    /* update the pointer */
          if (numRealChan<0)
            numRealChan = 0;
          if (numTrigChan<0)
            numTrigChan = 0;
          if (numTrigChan>MAXNUMTRIGCHAN) {
            display_print("cannot maintain more than %d trigger channels in real time\n", MAXNUMTRIGCHAN);
            numTrigChan = MAXNUMTRIGCHAN;
          }
          break;

        case ACQ_MSGQ_DATA:
          countData++;
          if (packet[i].sampleNumber<minSampleValue) {
            minSampleIndex = i;
            minSampleValue = packet[i].sampleNumber;
          }
          if (packet[i].sampleNumber>maxSampleValue) {
            maxSampleIndex = i;
            maxSampleValue = packet[i].sampleNumber;
          }
          if (packet[i].messageId<minMessageValue) {
            minMessageIndex = i;
            minMessageValue = packet[i].messageId;
          }
          if (packet[i].messageId>maxMessageValue) {
            maxMessageIndex = i;
            maxMessageValue = packet[i].messageId;
          }
          if ((packet[i].sampleNumber/ packet[i].numSamples)>lastPacket) {
            /* detect the flanks in the trigger channels */
            lastPacket = packet[i].sampleNumber/ packet[i].numSamples;
            for (j=0; j<numTrigChan; j++) {
              if (trigChan[j]>-1 && trigChan[j]<numRealChan)
                for (k=0; k<packet[i].numSamples; k++) {
                  int sample = k*numRealChan + trigChan[j];
                  /* detect changes in the value of the trigger channel */
                  if (packet[i].data[sample]!=lastValue[j]) {
                    lastValue[j] = packet[i].data[sample];
                    if (lastValue[j]) {
                      display_print("trigger detected\n");
                      /* only store it if it is an upgoing flank */
                      trigPointer[numTrigger+0] = trigChan[j];              /* channel number */
                      trigPointer[numTrigger+1] = packet[i].sampleNumber+k; /* sample number  */
                      trigPointer[numTrigger+2] = lastValue[j];             /* sample value   */
                      numTrigger = wrapnumtrigger(numTrigger+3);
                    }
                  }
                }
            }
          }
          break;

        case ACQ_MSGQ_CLOSE_CONNECTION:
          countCancel++;
          break;

        case ACQ_MSGQ_INVALID:
          countInvalid++;
          break;

        default:
          countUnknown++;
          clear_packet(i);
          countCleared++;
          break;

      } /* end switch */
    } /* end for */

    /*******************************************************************************
     * print information about all packets in the buffer
     *******************************************************************************/

    display_print("\n");
    display_print("buffer size  = %d\n", ACQ_MSGQ_SIZE);
    display_print("shm size     = %d\n", shmsize);
    display_print("shm key      = %#x\n", shmkey);
    display_print("pause        = %d\n", pause);
    if (setupIndex>=0) {
      display_print("dataset      = %s @ %d\n", (char *)packet[setupIndex].data, setupIndex);
    }
    else {
      display_print("dataset      = <unknown>\n");
    }
    display_print("\n");

    display_print("countSetup   = %d\n", countSetup);
    display_print("countData    = %d\n", countData);
    display_print("countCancel  = %d\n", countCancel);
    display_print("countInvalid = %d\n", countInvalid);
    display_print("countUnknown = %d\n", countUnknown);
    display_print("countCleared = %d\n", countCleared);
    display_print("\n");

    /* this might look like a weird location to reinitialize */
    maxCountCleared = (countCleared>maxCountCleared ? countCleared : maxCountCleared);
    countCleared = 0;

    display_print("min(sampleNumber) = %d @ %d\n", minSampleValue, minSampleIndex);
    display_print("max(sampleNumber) = %d @ %d\n", maxSampleValue, maxSampleIndex);
    /* display_print("min(messageId)    = %d @ %d\n", minMessageValue, minMessageIndex);
    display_print("max(messageId)    = %d @ %d\n", maxMessageValue, maxMessageIndex);*/
    display_print("max(countCleared) = %d\n", maxCountCleared);
    display_print("\n");

    if (maxSampleIndex>=0) {
      display_print("current trial = %d @ %d\n", maxSampleValue/packet[maxSampleIndex].numSamples, maxSampleIndex);
    }
    else {
      display_print("current trial = <unknown>\n");
    }
    display_print("\n");

    display_print("numRealChan  = %d\n", numRealChan);
    display_print("numTrigChan  = %d\n", numTrigChan);
    display_print("numTrigger   = %d\n", numTrigger/3);
    display_print("lastPacket   = %d\n", lastPacket);
    display_print("trigChan     = ");
    for (i=0; i<numTrigChan; i++)
      display_print("%d\t", trigChan[i]);
    display_print("\n");
    display_print("lastValue    = ");
    for (i=0; i<numTrigChan; i++)
      display_print("%d\t", lastValue[i]);
    display_print("\n");
    display_print("\n");

    display_refresh();

    /*******************************************************************************
     * do the desired maintenance on the packet buffer (only if a key was pressed)
     *******************************************************************************/

    k = getch();
    switch (k) {
    case 's':
      write_setup(0);
      break;

    case 'p':
      pause = PAUSE;
      break;

    case 'f':
      while (countInvalid && packet[t].message_type!=ACQ_MSGQ_INVALID)
        t = wrapnumpacket(t+1);
      while (packet[t].message_type==ACQ_MSGQ_INVALID) {
        write_data(t, c);
        t = wrapnumpacket(t+1);
        c++;
      }
      break;

    case 'd':
      while (countInvalid && packet[t].message_type!=ACQ_MSGQ_INVALID)
        t = wrapnumpacket(t+1);
      write_data(t, c);
      t = wrapnumpacket(t+1);
      c++;
      break;

    case 'c':
      maxCountCleared = 0;
      clear_all_except_setup();
      break;

    case 'h':
      showDisplay = (!showDisplay);
      break;

    case 'x':
      maxCountCleared = 0;
      clear_all();
      break;

    case 'q':
      again = 0;
      break;
    }
 
    if (k!=ERR) {
      display_print("key pressed = %d\n", k);
      display_refresh();
      continue;
    }

    /*******************************************************************************
     * do the regular maintenance on the packet buffer
     *******************************************************************************/

    if (countCancel) {
      display_print("initializing all packets\n");
      clear_all();
      countInvalid = ACQ_MSGQ_SIZE;
      countData    = 0;
      countSetup   = 0;
      countCancel  = 0;
    }

    /* if ((MINFREE-countInvalid)>1 && pause) {
      pause/=(MINFREE-countInvalid);
      display_print("decreasing pause to %d\n", pause);
    }
    */

    while (countInvalid<MINFREE && countData)
      /* make more empty packets available */
      if (setupIndex>-1) {
        /* moving the setup one packet at a time may look inefficient, but this only happens at the start of online acquisition */
        display_print("moving setup from %d to %d\n", setupIndex, minSampleIndex);
        memcpy(&packet[minSampleIndex], &packet[setupIndex], sizeof(ACQ_MessagePacketType));
        countCleared++;
        countInvalid++;
        countData--;
        /* NOTE: don't clear the packet, since Matlab might be reading from it
        display_print("clearing packet %d\n", setupIndex);
        clear_packet(setupIndex);
        */
        packet[setupIndex].message_type = ACQ_MSGQ_INVALID;
        setupIndex     = minSampleIndex;
        minSampleIndex = wrapnumpacket(minSampleIndex+1); /* the next is probably now the oldest */
      }
      else {
        display_print("clearing packet %d\n", minSampleIndex);
        clear_packet(minSampleIndex);
        countCleared++;
        countInvalid++;
        countData--;
        minSampleIndex = wrapnumpacket(minSampleIndex+1); /* the next is probably now the oldest */
      }

    while (setupIndex>-1 && countData && packet[wrapnumpacket(setupIndex+1)].message_type==ACQ_MSGQ_INVALID) {
      /* move the setup to the next empty packet */
      /* moving the setup one packet at a time may look inefficient, but this only happens at the start of online acquisition */
      display_print("moving setup from %d to %d\n", setupIndex, wrapnumpacket(setupIndex+1));
      memcpy(&packet[wrapnumpacket(setupIndex+1)], &packet[setupIndex], sizeof(ACQ_MessagePacketType));
      /* NOTE: don't clear the packet, since Matlab might be reading from it
      clear_packet(setupIndex);
      */
      packet[setupIndex].message_type = ACQ_MSGQ_INVALID;
      countCleared++;
      setupIndex = wrapnumpacket(setupIndex+1);
    }

    display_refresh();
    usleep(pause);

    if (countInvalid>4)
      printf("ok %d\n", countLoop);
    else
      printf("error %d\n", countLoop);
    countLoop++;
  }

  /* detach from the shared memory segment */
  if (shmdt(packet) == -1) {
    perror("shmdt");
    exit(1);
  }

  /* end curses mode */
  display_stop();
  
  /* end of program */
  return 0;
}

void write_setup(int i) {
  /* specify the dataset name in the first packet */
  packet[i].message_type = ACQ_MSGQ_SETUP_COLLECTION;
  packet[i].messageId    = 0;
  packet[i].sampleNumber = 0;
  packet[i].numSamples   = 0;
  packet[i].numChannels  = 0;
  memset(packet[i].data, 0, DATASIZE*sizeof(int));
  strcpy((char *)(packet[i].data), "/home/coherence/roboos/MEG/jansch_roboos_20060209_01.ds");
}

void write_data(int i, int j) {
  /* specify the dataset name in the first packet */
  packet[i].message_type = ACQ_MSGQ_DATA;
  packet[i].messageId    = j;
  packet[i].sampleNumber = j;
  packet[i].numSamples   = 1;
  packet[i].numChannels  = 1;
  memset(packet[i].data, 0, DATASIZE*sizeof(int));
  packet[i].data[0] = j;  /* write something to identify the packet */
}

void clear_packet(int i) {
  /* clear the packet header and data */
  packet[i].message_type = ACQ_MSGQ_INVALID;
  packet[i].messageId    = 0;
  packet[i].sampleNumber = 0;
  packet[i].numSamples   = 0;
  packet[i].numChannels  = 0;
  memset(packet[i].data, 0, DATASIZE*sizeof(int));
}

void clear_all(void) {
  int i;
  for (i=0; i<ACQ_MSGQ_SIZE; i++)
    clear_packet(i);
}

void clear_all_except_setup(void) {
  int i, j = 0;
  char *ptr;
  for (i=0; i<ACQ_MSGQ_SIZE; i++) 
    if (packet[i].message_type==ACQ_MSGQ_SETUP_COLLECTION) {
      /* keep the dataset name intact but clear the rest (including the trigger details) */
      ptr = (char *)packet[i].data;
      while (ptr[j]) j++;
      /* while (j<(DATASIZE*sizeof(int))) {ptr[j]=0; j++;}; */
      while (j<((DATASIZE-9012)*sizeof(int))) {ptr[j]=0; j++;};
      j+=12*sizeof(int);
      while (j<(DATASIZE)*sizeof(int)) {ptr[j]=0; j++;};
      numTrigger = 0;
    }
    else {
      /* clear the packet header and data */
      clear_packet(i);
    }
}

void print_packet(int i) {
  /* print the details of each package */
  display_print("%4d : ", i);
  switch (packet[i].message_type) {
    case ACQ_MSGQ_SETUP_COLLECTION:
      display_print("setup  \t");
      break;
    case ACQ_MSGQ_DATA:
      display_print("data  \t");
      break;
    case ACQ_MSGQ_CLOSE_CONNECTION:
      display_print("close  \t");
      break;
    case ACQ_MSGQ_INVALID:
      display_print("-      \t");
      break;
    default:
      display_print("%d\t", packet[i].message_type);
  }
  display_print("%d\t", packet[i].messageId);
  display_print("%d\t", packet[i].sampleNumber);
  display_print("%d\t", packet[i].numSamples);
  display_print("%d\t", packet[i].numChannels);
  display_print("\n");
}

void display_start(void) {
  initscr();                      /* Start curses mode                */
  cbreak();                       /* Get one character at a time      */
  noecho();                       /* Don't echo() while we do getch   */
  keypad(stdscr, TRUE);           /* We get F1, F2 etc..              */
  nodelay(stdscr, TRUE);          /* Return ERR if no key was pressed */
}

void display_stop(void) {
  endwin();                       /* Stop curses mode                 */
}

void display_refresh(void) {
  if (showDisplay)
    refresh();
}

