/*
 * Copyright (C) 2008, Robert Oostenveld
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 * $Log: buffer_getevt.c,v $
 * Revision 1.10  2009/06/30 11:35:48  roboos
 * pass more explicit error from buffer onto matlab
 *
 * Revision 1.9  2008/10/29 20:46:12  roboos
 * consistent use of open_connection and close_connection
 * there were some incorrect uses of close() that on windows did not actually close, resulting in the buffer running out of sockets/threads after prolonged use
 *
 * Revision 1.8  2008/07/09 13:34:21  roboos
 * small change in verbose output, using verbose=0|1
 *
 * Revision 1.7  2008/05/22 09:27:22  roboos
 * fixed some issues with Borland compiler, correct pointer arithmetic, declarations at beginning
 *
 * Revision 1.6  2008/03/23 12:47:41  roboos
 * implemented selection of events
 *
 * Revision 1.5  2008/03/09 22:36:38  roboos
 * implemented new buffer mex-function, which acts as dispatcher for all functionality
 * changed buffer_xxx code from mex-function into c-function to be included in dispatcher
 *
 * Revision 1.4  2008/02/26 21:43:25  roboos
 * renamed packet_t structure definition into messagedef_t, added message_t structure (contains def+buf)
 *
 * Revision 1.3  2008/02/20 13:49:16  roboos
 * changed somments into ansi style, needed for matlab
 *
 * Revision 1.2  2008/02/19 10:28:33  roboos
 * added copyright statement and log message
 *
 *
 */

#include "mex.h"
#include "matrix.h"
#include "buffer.h"

#define NUMBER_OF_FIELDS 9

void buffer_getevt(char *hostname, int port, mxArray *plhs[], const mxArray *prhs[])
{
  int server;
  int verbose = 0;
  int i, nevents;
  int offset;
  double *val;
  char msg[512];
  
  mxArray *bufptr;
  
  
  message_t *request  = NULL;
  message_t *response = NULL;
  event_t   *event    = NULL;
  eventsel_t eventsel;
  
  /* this is for the Matlab specific output */
  const char *field_names[] = {
    "type_type",
    "type_numel",
    "value_type",
    "value_numel",
    "sample",
    "offset",
    "duration",
    "bufsize",
    "buf"
  };
  
  /* allocate the elements that will be used in the communication */
  request      = malloc(sizeof(message_t));
  request->def = malloc(sizeof(messagedef_t));
  request->buf = NULL;
  request->def->version = VERSION;
  request->def->command = GET_EVT;
  request->def->bufsize = 0;
  
  if ((prhs[0]!=NULL) && (mxGetNumberOfElements(prhs[0])==2) && (mxIsDouble(prhs[0])) && (!mxIsComplex(prhs[0]))) {
    /* fprintf(stderr, "args OK\n"); */
    val = (double *)mxGetData(prhs[0]);
    eventsel.begevent = (UINT32_T)(val[0]);
    eventsel.endevent = (UINT32_T)(val[1]);
    if (verbose) print_eventsel(&eventsel);
    request->def->bufsize = append(&request->buf, request->def->bufsize, &eventsel, sizeof(eventsel_t));
  }
  
  /* open the TCP socket */
  if ((server = open_connection(hostname, port)) < 0) {
    sprintf(msg, "ERROR: failed to create socket (%d)\n", server);
		mexErrMsgTxt(msg);
  }
  
  if (verbose) print_request(request->def);
  clientrequest(server, request, &response);
  if (verbose) print_response(response->def);
  close_connection(server);
  
  if (response->def->command==GET_OK) {
    event = malloc(sizeof(event_t));
    
    /* first count the number of events */
    nevents = 0;
    offset = 0;
    while (offset<response->def->bufsize) {
      event->def = (char *)response->buf + offset;
      event->buf = (char *)response->buf + offset + sizeof(eventdef_t);
      if (verbose) print_eventdef(event->def);
      offset += sizeof(eventdef_t) + event->def->bufsize;
      nevents++;
    }
    
    /* create a structure array that can hold all events */
    plhs[0] = mxCreateStructMatrix(1, nevents, NUMBER_OF_FIELDS, field_names);
    
    offset = 0;
    for (i=0; i<nevents; i++) {
      event->def = (char *)response->buf + offset;
      event->buf = (char *)response->buf + offset + sizeof(eventdef_t);
      
      bufptr = mxCreateNumericMatrix(1, event->def->bufsize, mxUINT8_CLASS, mxREAL);
      memcpy(mxGetPr(bufptr), event->buf, event->def->bufsize);
      
      mxSetFieldByNumber(plhs[0], i, 0, mxCreateDoubleScalar((double)event->def->type_type));
      mxSetFieldByNumber(plhs[0], i, 1, mxCreateDoubleScalar((double)event->def->type_numel));
      mxSetFieldByNumber(plhs[0], i, 2, mxCreateDoubleScalar((double)event->def->value_type));
      mxSetFieldByNumber(plhs[0], i, 3, mxCreateDoubleScalar((double)event->def->value_numel));
      mxSetFieldByNumber(plhs[0], i, 4, mxCreateDoubleScalar((double)event->def->sample));
      mxSetFieldByNumber(plhs[0], i, 5, mxCreateDoubleScalar((double)event->def->offset));
      mxSetFieldByNumber(plhs[0], i, 6, mxCreateDoubleScalar((double)event->def->duration));
      mxSetFieldByNumber(plhs[0], i, 7, mxCreateDoubleScalar((double)event->def->bufsize));
      mxSetFieldByNumber(plhs[0], i, 8, bufptr);
      offset += sizeof(eventdef_t) + event->def->bufsize;
    }
    FREE(event);
  }
  else {
    sprintf(msg, "ERROR: the buffer returned an error (%d)\n", response->def->command);
		mexErrMsgTxt(msg);
  }
  
  if (request) {
    FREE(request->def);
    FREE(request->buf);
    FREE(request);
  }
  if (response) {
    FREE(response->def);
    FREE(response->buf);
    FREE(response);
  }
  
  return;
}

