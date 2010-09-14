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

#include "buffer_mxutils.h"

#define NUMBER_OF_FIELDS 5



int buffer_getevt(int server, mxArray *plhs[], const mxArray *prhs[])
{
  int verbose = 0;
  int i, nevents;
  int offset;
  double *val;
  int result;
  
  mxArray *bufptr;
  
  
  message_t *request  = NULL;
  message_t *response = NULL;
  eventsel_t eventsel;
  
  /* this is for the Matlab specific output */
  const char *field_names[] = {
    "type",
    "value",
    "sample",
    "offset",
    "duration"
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
  
  if (verbose) print_request(request->def);
  result = clientrequest(server, request, &response);
  if (verbose) print_response(response->def);
  
  if (result == 0) {
    if (response->def->command==GET_OK) {
	  eventdef_t *event_def;
      /* first count the number of events */
      nevents = 0;
      offset = 0;
      while (offset<response->def->bufsize) {
        event_def = (eventdef_t *)((char *)response->buf + offset);
        /* event_buf = (char *)response->buf + offset + sizeof(eventdef_t); */
        if (verbose) print_eventdef(event_def);
        offset += sizeof(eventdef_t) + event_def->bufsize;
        nevents++;
      }
      
      /* create a structure array that can hold all events */
      plhs[0] = mxCreateStructMatrix(1, nevents, NUMBER_OF_FIELDS, field_names);
      
      offset = 0;
      for (i=0; i<nevents; i++) {
	  	char *buf_type,*buf_value;
		
        event_def = (eventdef_t *) ((char *)response->buf + offset);
		buf_type = (char *) response->buf + offset + sizeof(eventdef_t);
		buf_value = buf_type + event_def->type_numel * wordsize_from_type(event_def->type_type);
                
		mxSetFieldByNumber(plhs[0], i, 0, matrix_from_ft_type_data(event_def->type_type, 1, event_def->type_numel, buf_type));
		mxSetFieldByNumber(plhs[0], i, 1, matrix_from_ft_type_data(event_def->value_type, 1, event_def->value_numel, buf_value));
        mxSetFieldByNumber(plhs[0], i, 2, mxCreateDoubleScalar((double)event_def->sample+1)); /* 1-based in Matlab, 0-based in protocol */
        mxSetFieldByNumber(plhs[0], i, 3, mxCreateDoubleScalar((double)event_def->offset));
        mxSetFieldByNumber(plhs[0], i, 4, mxCreateDoubleScalar((double)event_def->duration));
        offset += sizeof(eventdef_t) + event_def->bufsize;
      }
    }
    else {
      result = response->def->command;
    }
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
  
  return result;
}

