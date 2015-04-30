/**********************************************************************
 * MEX wrapper for PortMidi input functions
 *
 * Copyright (C) 2015, Robert Oostenveld
 *
 * This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
 * for the documentation and details.
 *
 *    FieldTrip is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    FieldTrip is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
 *
 * $Id$
 **********************************************************************/

#include <portmidi.h>
#include <porttime.h>
#include <mex.h>

#define FREE(x) {if (x) {free(x); x=NULL;}}
#define WRAP(x,y) ((x) - ((int)((float)(x)/(y)))*(y))
#define INPUT_BUFFER_SIZE 20000

/* these were taken from portmidi/pm_test/mm.c */
#define MIDI_CODE_MASK          0xf0
#define MIDI_CHN_MASK           0x0f
#define MIDI_OFF_NOTE           0x80
#define MIDI_ON_NOTE            0x90
#define MIDI_POLY_TOUCH         0xa0
#define MIDI_CTRL               0xb0
#define MIDI_CH_PROGRAM         0xc0
#define MIDI_TOUCH              0xd0
#define MIDI_BEND               0xe0

#define MIDI_SYSEX              0xf0
#define MIDI_Q_FRAME            0xf1
#define MIDI_SONG_POINTER       0xf2
#define MIDI_SONG_SELECT        0xf3
#define MIDI_TUNE_REQ           0xf6
#define MIDI_EOX                0xf7
#define MIDI_TIME_CLOCK         0xf8
#define MIDI_START              0xfa
#define MIDI_CONTINUE           0xfb
#define MIDI_STOP               0xfc
#define MIDI_ACTIVE_SENSING     0xfe
#define MIDI_SYS_RESET          0xff

#define MIDI_ALL_SOUND_OFF      0x78
#define MIDI_RESET_CONTROLLERS  0x79
#define MIDI_LOCAL              0x7a
#define MIDI_ALL_OFF            0x7b
#define MIDI_OMNI_OFF           0x7c
#define MIDI_OMNI_ON            0x7d
#define MIDI_MONO_ON            0x7e
#define MIDI_POLY_ON            0x7f

PortMidiStream *inStream = NULL;
int isInit = 0, deviceOpen = -1, verbose = 1;
unsigned int numReceived = 0;
double *channel = NULL, *note = NULL, *velocity = NULL, *timestamp = NULL;

void reportPmError(PmError err);
void reportPtError(PtError err);

void receive_poll(PtTimestamp ts, void *userData)
{
  int count, command;
  unsigned int latest;
  PmEvent event;
  
  if (inStream==NULL)
    return;
  
  if (Pm_Poll(inStream)!=TRUE)
    return;
  
  while ((count = Pm_Read(inStream, &event, 1))) {
    if (count == 1) {
      
      /* there seems to be a constant stream of MIDI events, not all of which are interesting */
      /* the status has the command in the highest 4 bits and the channel in the lowest 4 bits */
      command = Pm_MessageStatus(event.message) & MIDI_CODE_MASK;
      if ((command == MIDI_ON_NOTE)         ||
              (command == MIDI_OFF_NOTE )   ||
              (command == MIDI_CH_PROGRAM ) ||
              (command == MIDI_CTRL)        ||
              (command == MIDI_POLY_TOUCH ) ||
              (command == MIDI_TOUCH )      ||
              (command == MIDI_BEND )) {
        
        /* store the latest event at the end of the ring buffer */
        latest = WRAP(numReceived, INPUT_BUFFER_SIZE);
        numReceived++;
        
        channel   [latest] = (double)(Pm_MessageStatus(event.message) & MIDI_CHN_MASK);
        note      [latest] = (double)Pm_MessageData1(event.message);
        velocity  [latest] = (double)Pm_MessageData2(event.message);
        timestamp [latest] = (double)event.timestamp;
        
        if (verbose) {
          mexPrintf("channel = %2d, note = %3d, velocity = %3d, timestamp = %g\n", (int)channel[latest], (int)note[latest], (int)velocity[latest], timestamp[latest]);
          mexEvalString("try, drawnow limitrate nocallbacks; end");
        }
      }
    }
    else
      mexWarnMsgTxt(Pm_GetErrorText(count));
  }
}

void exitFunction() {
  if (isInit) {
    mexPrintf("Terminating PortMidi\n");
    reportPtError(Pt_Stop());
    if (inStream != NULL) {
      reportPmError(Pm_Close(inStream));
      inStream = NULL;
    }
    Pm_Terminate();
    FREE(channel);
    FREE(note);
    FREE(velocity);
    FREE(timestamp);
    isInit = 0;
    deviceOpen = -1;
    numReceived = 0;
  }
}

void init() {
  PmError err1;
  PtError err2;
  
  if (isInit) return;
  
  mexPrintf("Initialising PortMidi\n");
  
  note      = malloc(INPUT_BUFFER_SIZE*sizeof(double));
  channel   = malloc(INPUT_BUFFER_SIZE*sizeof(double));
  velocity  = malloc(INPUT_BUFFER_SIZE*sizeof(double));
  timestamp = malloc(INPUT_BUFFER_SIZE*sizeof(double));
  if (!(channel && note && velocity && timestamp)) {
    FREE(channel);
    FREE(note);
    FREE(velocity);
    FREE(timestamp);
    mexErrMsgTxt("Could not allocate memory");
  }
  
  err1 = Pm_Initialize();
  reportPmError(err1);
  
  err2 = Pt_Start(1, receive_poll, NULL);
  reportPtError(err2);
  
  /* getting here means that PortMidi and PortTime are both fine */
  mexAtExit(exitFunction);
  isInit = 1;
}

void reportPtError(PtError err) {
  switch(err) {
    case pmNoError:
      return;
    case ptHostError:
      mexErrMsgTxt("PortTime error: Host error");
    case ptAlreadyStarted:
      mexErrMsgTxt("PortTime error: Cannot start timer because it is already started");
    case ptAlreadyStopped:
      mexErrMsgTxt("PortTime error: Cannot stop timer because it is already stopped");
    case ptInsufficientMemory:
      mexErrMsgTxt("PortTime error: Memory could not be allocated");
  }
}

void reportPmError(PmError err) {
  switch(err) {
    case pmNoError:
      return;
    case pmGotData:
      return;
    case pmHostError:
      mexErrMsgTxt("PortMidi error: Host error");
    case pmInvalidDeviceId:
      mexErrMsgTxt("PortMidi error: Invalid device ID");
    case pmInsufficientMemory:
      mexErrMsgTxt("PortMidi error: Insufficient memory");
    case pmBufferTooSmall:
      mexErrMsgTxt("PortMidi error: Buffer too small");
    case pmBufferOverflow:
      mexErrMsgTxt("PortMidi error: Buffer overflow");
    case pmBadPtr:
      mexErrMsgTxt("PortMidi error: Bad pointer");
    case pmBadData:
      mexErrMsgTxt("PortMidi error: Bad data");
    case pmInternalError:
      mexErrMsgTxt("PortMidi error: Internal error");
    case pmBufferMaxSize:
      mexErrMsgTxt("PortMidi error: Buffer is already as large as it can be");
  }
}

mxArray *getDevices() {
  const char *fieldNames[] = {"index", "name", "input", "output"};
  int i,numDevs;
  const PmDeviceInfo* info;
  mxArray *R;
  
  numDevs = Pm_CountDevices();
  R = mxCreateStructMatrix(numDevs, 1, 4, fieldNames);
  
  for (i=0;i<numDevs;i++) {
    info = Pm_GetDeviceInfo(i);
    mxSetFieldByNumber(R, i, 0, mxCreateDoubleScalar(i+1));
    mxSetFieldByNumber(R, i, 1, mxCreateString(info->name));
    mxSetFieldByNumber(R, i, 2, mxCreateDoubleScalar(info->input));
    mxSetFieldByNumber(R, i, 3, mxCreateDoubleScalar(info->output));
  }
  return R;
}

mxArray *getEvent() {
  unsigned int i, n;
  double *ptr;
  mxArray *R;
  n = numReceived; /* should not increase while copying */
  R = mxCreateDoubleMatrix(n, 4, mxREAL);
  ptr = (double *)mxGetData(R);
  for (i=0; i<n; i++) {
    ptr[i + 0*n] = channel  [i];
    ptr[i + 1*n] = note     [i];
    ptr[i + 2*n] = velocity [i];
    ptr[i + 3*n] = timestamp[i];
  }
  return R;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  int code = 'L';
  PmError err1;
  PtError err2;
  
  init(); /* returns immediately if already done */
  
  if (nrhs > 0) {
    char arg[2];
    if (!mxIsChar(prhs[0]) || mxGetNumberOfElements(prhs[0]) > 1) {
      mexErrMsgTxt("Bad call\n");
    }
    mxGetString(prhs[0], arg, 2);
    code = arg[0];
  }
  
  switch(code) {
    case 'V':
      if (inStream == NULL) {
        mexErrMsgTxt("No MIDI input device is opened");
      } else {
        if (nrhs < 2) mexErrMsgTxt("Usage for 'verbose': midiIn('V', value)");
        verbose = mxGetScalar(prhs[1]);
      }
      break;
    case 'G': /* Return the buffered events as array */
      if (inStream == NULL) {
        mexErrMsgTxt("No MIDI input device is opened");
      } else {
        plhs[0] = getEvent();
        /* reset the counter to the beginning of the buffer */
        numReceived = 0;
      }
      break;
    case 'F': /* Flush all buffered events */
      if (inStream == NULL) {
        mexErrMsgTxt("No MIDI input device is opened");
      } else {
        /* reset the counter to the beginning of the buffer */
        numReceived = 0;
      }
      break;
    case 'C':   /* Close input stream */
      if (inStream == NULL) {
        mexWarnMsgTxt("No MIDI input device is opened - ignoring 'close' command");
      } else {
        err1 = Pm_Close(inStream);
        inStream = NULL;
        if (err1 != pmNoError) reportPmError(err1);
      }
      break;
    case 'O':	/* Open input stream */
    {
      int device = 0;
      
      if (nrhs<2 || !mxIsNumeric(prhs[1])) mexErrMsgTxt("Bad call\n");
      device = (int) mxGetScalar(prhs[1]);
      if (device < 1 || device > Pm_CountDevices()) mexErrMsgTxt("Device index out of range");
      --device;
      
      if (inStream != NULL) {
        if (deviceOpen == device) {
          mexWarnMsgTxt("MIDI input device is already open - ignoring request");
          return;
        }
        mexWarnMsgTxt("Another MIDI input device is open - closing that one");
        err1 = Pm_Close(inStream);
        inStream = NULL;
        if (err1 != pmNoError) reportPmError(err1);
      }
      
      err1 = Pm_OpenInput(&inStream, device, NULL, INPUT_BUFFER_SIZE, NULL, NULL);
      if (err1 != pmNoError) {
        inStream = NULL;
        reportPmError(err1);
      }
      
      /* remember which device we just opened */
      deviceOpen = device;
    }
    break;
    case 'L':	/* List devices */
      plhs[0] = getDevices();
      break;
    default:
      mexErrMsgTxt("Bad call\n");
  }
}
