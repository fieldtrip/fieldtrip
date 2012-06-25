/*
 * Copyright (C) 2008, Robert Oostenveld
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 * $Id$
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

