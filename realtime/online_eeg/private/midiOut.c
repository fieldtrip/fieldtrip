/**********************************************************************
 * MEX wrapper for PortMidi output functions
 *
 * Copyright (C) 2010, Stefan Klanke
 *
 * This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
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

#define OUTPUT_BUFFER_SIZE 20000

void reportPmError(PmError err);

int isInit = 0, openID = -1;
PortMidiStream *outStream = NULL;

void exitFunction() {
  if (isInit) {
    printf("Terminating PortMidi\n");
    if (outStream != NULL) {
      Pm_Close(outStream);
      outStream = NULL;
    }
    Pt_Stop();
    Pm_Terminate();
    isInit = 0;
  }
}

void init() {
  PmError err;
  
  if (isInit) return;
  
  printf("Initialising PortMidi\n");
  err = Pm_Initialize();
  reportPmError(err);
  
  if (Pt_Start(1,NULL,NULL) == 0) {
    mexAtExit(exitFunction);
    isInit = 1;
  } else {
    Pm_Terminate();
    mexErrMsgTxt("Could start PortMidi, but not PortTime.");
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
      mexErrMsgTxt("PortMidi error: Insufficient memory");
    case pmBufferMaxSize:
      mexErrMsgTxt("PortMidi error: Buffer is already as large as it can be");
  }
}

void sendNoteOnOff(int status, const mxArray *prhs[]) {
  int n,N;
  const double *note, *vel;
  int channel = (int) mxGetScalar(prhs[1]);
  
  if (channel < 1 || channel > 16) mexErrMsgTxt("Channel (2nd arg.) must be 1..16");
  --channel;
  
  if (!mxIsDouble(prhs[2]) || !mxIsDouble(prhs[3])) {
    mexErrMsgTxt("3rd and 4th arguments (note number and velocity) must be 'double' arrays");
  }
  N = mxGetNumberOfElements(prhs[2]);
  if (N != mxGetNumberOfElements(prhs[3])) {
    mexErrMsgTxt("3rd and 4th arguments (note number and velocity) must have same number of elements");
  }
  note = mxGetPr(prhs[2]);
  vel = mxGetPr(prhs[3]);
  
  status = (status & 0xF0) | channel;
  
  for (n=0;n<N;n++) {
    PmMessage msg;
    PmError err;
    
    msg = Pm_Message(status, (int) note[n], (int) vel[n]);
    /* latency=0 => send asap */
    err = Pm_WriteShort(outStream, 0, msg);
    if (err!=pmNoError) reportPmError(err);
  }
}


void sendProgChange(const mxArray *prhs[]) {
  int channel = (int) mxGetScalar(prhs[1]);
  int prog = (int) mxGetScalar(prhs[2]);
  int status;
  PmMessage msg;
  PmError err;
  
  if (channel < 1 || channel > 16) mexErrMsgTxt("Channel (2nd arg.) must be 1..16");
  --channel;
  
  if (!mxIsDouble(prhs[2])) {
    mexErrMsgTxt("3rd argument (program number) must be a 'double'");
  }
  
  status = 0xC0 | channel;
  
  msg = Pm_Message(status, prog-1, 0);
  /* latency=0 => send asap */
  err = Pm_WriteShort(outStream, 0, msg);
  if (err!=pmNoError) reportPmError(err);
}

void playRaw(const mxArray *A) {
  int N = mxGetNumberOfElements(A);
  unsigned char * data = mxGetData(A);
  int n;
  
  PmMessage msg;
  PmError err;
  
  if (N % 3 != 0) mexErrMsgTxt("Raw arguments must be a multiple of 3 elements long.");
  
  N/=3;
  
  for (n=0;n<N;n++) {
    msg = Pm_Message(data[0], data[1], data[2]);
    err = Pm_WriteShort(outStream, 0, msg);
    if (err!=pmNoError) reportPmError(err);
    data+=3;
  }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  int code = 'L';
  PmError err;
  
  init(); /* returns immediately if already done */
  
  if (nrhs > 0) {
    char arg[2];
    
    if (mxIsUint8(prhs[0])) {
      if (outStream == NULL) {
        mexErrMsgTxt("No MIDI output device is opened");
      } else {
        playRaw(prhs[0]);
      }
      return;
    }
    
    if (!mxIsChar(prhs[0]) || mxGetNumberOfElements(prhs[0]) > 1) {
      mexErrMsgTxt("Bad call\n");
    }
    mxGetString(prhs[0], arg, 2);
    code = arg[0];
  }
  
  switch(code) {
    case 'P':
      if (outStream == NULL) {
        mexErrMsgTxt("No MIDI output device is opened");
      } else {
        if (nrhs < 3) mexErrMsgTxt("Usage for 'program change': midiOut('P', channel, program)");
        sendProgChange(prhs);
      }
      break;
    case '+':	/* Note on */
      if (outStream == NULL) {
        mexErrMsgTxt("No MIDI output device is opened");
      } else {
        if (nrhs < 4) mexErrMsgTxt("Usage for 'note on': midiOut('+', channel, notes, velocities)");
        sendNoteOnOff(0x90, prhs);
      }
      break;
    case '-': 	/* Note off */
      if (outStream == NULL) {
        mexErrMsgTxt("No MIDI output device is opened");
      } else {
        if (nrhs < 4) mexErrMsgTxt("Usage for 'note on': midiOut('+', channel, notes, velocities)");
        sendNoteOnOff(0x80, prhs);
      }
      break;
    case '.':	/* All notes off */
      if (outStream == NULL) {
        mexErrMsgTxt("No MIDI output device is opened");
      } else {
        int channel;
        
        if (nrhs < 2) mexErrMsgTxt("Usage for 'all notes off': midiOut('.', channel)");
        channel = (int) mxGetScalar(prhs[1]);
        if (channel < 1 || channel > 16) mexErrMsgTxt("Channel (2nd arg.) must be 1..16");
        --channel;
        err = Pm_WriteShort(outStream, 0, Pm_Message(0xB0 | channel, 123, 0));
        if (err!=pmNoError) reportPmError(err);
      }
      break;
    case 'C':   /* Close output stream */
      if (outStream == NULL) {
        mexWarnMsgTxt("No MIDI output device is opened - ignoring 'close' command");
      } else {
        err = Pm_Close(outStream);
        outStream = NULL;
      }
      break;
    case 'O':	/* Open output stream */
    {
      int device = 0;
      
      
      if (nrhs<2 || !mxIsNumeric(prhs[1])) mexErrMsgTxt("Bad call\n");
      device = (int) mxGetScalar(prhs[1]);
      if (device < 1 || device > Pm_CountDevices()) mexErrMsgTxt("Device index out of range");
      --device;
      
      if (outStream != NULL) {
        if (openID == device) {
          mexWarnMsgTxt("MIDI output device is already open - ignoring request");
          return;
        }
        mexWarnMsgTxt("Another MIDI output device is open - closing that one");
        err = Pm_Close(outStream);
        outStream = NULL;
        if (err != pmNoError) reportPmError(err);
      }
      
      /* last parameter = latency = 0 means that timestamps are ignored */
      err = Pm_OpenOutput(&outStream, device, NULL, OUTPUT_BUFFER_SIZE, NULL, NULL, 0);
      if (err != pmNoError) {
        outStream = NULL;
        reportPmError(err);
      }
      openID = device;
    }
    break;
    case 'L':	/* List devices */
      plhs[0] = getDevices();
      break;
    default:
      mexErrMsgTxt("Bad call\n");
  }
}
