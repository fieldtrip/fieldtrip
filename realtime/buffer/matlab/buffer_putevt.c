/*
 * Copyright (C) 2008, Robert Oostenveld
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 * $Log: buffer_putevt.c,v $
 * Revision 1.7  2009/06/30 11:35:48  roboos
 * pass more explicit error from buffer onto matlab
 *
 * Revision 1.6  2008/10/29 20:46:12  roboos
 * consistent use of open_connection and close_connection
 * there were some incorrect uses of close() that on windows did not actually close, resulting in the buffer running out of sockets/threads after prolonged use
 *
 * Revision 1.5  2008/06/19 21:00:36  roboos
 * give more information in the error output messages
 *
 * Revision 1.4  2008/06/17 15:49:25  roboos
 * deal with bufsize and with buf (sofar evt.buf was not added to the message and the code stated a big FIXME)
 *
 * Revision 1.3  2008/03/23 12:51:23  roboos
 * moved command=xxx in code to better location, no functional place
 *
 * Revision 1.2  2008/03/09 22:36:38  roboos
 * implemented new buffer mex-function, which acts as dispatcher for all functionality
 * changed buffer_xxx code from mex-function into c-function to be included in dispatcher
 *
 * Revision 1.1  2008/03/07 14:46:05  roboos
 * new implementation, the events are not yet fully complete but the basic functionality has been tested
 *
 *
 */

#include "mex.h"
#include "matrix.h"
#include "buffer.h"

void buffer_putevt(char *hostname, int port, mxArray * plhs[], const mxArray * prhs[])
{
  size_t n;
  int server, fieldnumber;
  mxArray *field;
  char msg[512];
  
  message_t   *request  = NULL;
  message_t   *response = NULL;
  event_t     *event   = NULL;
  
  /* allocate the event */
  event      = malloc(sizeof(event_t));
  event->def = malloc(sizeof(eventdef_t));
  event->buf = NULL;
  event->def->bufsize   = 0;
  
  /* allocate the request message */
  request      = malloc(sizeof(message_t));
  request->def = malloc(sizeof(messagedef_t));
  request->buf = NULL;
  request->def->version = VERSION;
  request->def->command = PUT_EVT;
  request->def->bufsize = 0;
  
  /* define the event, it has the fields type_type type_numel value_type value_numel sample offset duration */
  
  /* FIXME loop over mutiple events in case of an event-array */
  
  fieldnumber = mxGetFieldNumber(prhs[0], "type_type");
  if (fieldnumber<0) {
    mexErrMsgTxt("field 'type_type' is missing");
    goto cleanup;
  }
  else
  {
    field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
    if (!mxIsNumeric(field) || mxIsEmpty(field)) {
      mexErrMsgTxt("invalid data type for 'type_type'");
      goto cleanup;
    }
    else
      event->def->type_type = (UINT32_T)mxGetScalar(field) ;
  }
  
  fieldnumber = mxGetFieldNumber(prhs[0], "type_numel");
  if (fieldnumber<0) {
    mexErrMsgTxt("field 'type_numel' is missing");
    goto cleanup;
  }
  else
  {
    field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
    if (!mxIsNumeric(field) || mxIsEmpty(field)) {
      mexErrMsgTxt("invalid data type for 'type_numel'");
      goto cleanup;
    }
    else
      event->def->type_numel = (UINT32_T)mxGetScalar(field) ;
  }
  
  fieldnumber = mxGetFieldNumber(prhs[0], "value_type");
  if (fieldnumber<0) {
    mexErrMsgTxt("field 'value_type' is missing");
    goto cleanup;
  }
  else
  {
    field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
    if (!mxIsNumeric(field) || mxIsEmpty(field)) {
      mexErrMsgTxt("invalid data type for 'value_type'");
      goto cleanup;
    }
    else
      event->def->value_type = (UINT32_T)mxGetScalar(field) ;
  }
  
  fieldnumber = mxGetFieldNumber(prhs[0], "value_numel");
  if (fieldnumber<0) {
    mexErrMsgTxt("field 'value_numel' is missing");
    goto cleanup;
  }
  else
  {
    field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
    if (!mxIsNumeric(field) || mxIsEmpty(field)) {
      mexErrMsgTxt("invalid data type for 'value_numel'");
      goto cleanup;
    }
    else
      event->def->value_numel = (UINT32_T)mxGetScalar(field) ;
  }
  
  fieldnumber = mxGetFieldNumber(prhs[0], "sample");
  if (fieldnumber<0) {
    mexErrMsgTxt("field 'sample' is missing");
    goto cleanup;
  }
  else
  {
    field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
    if (!mxIsNumeric(field) || mxIsEmpty(field)) {
      mexErrMsgTxt("invalid data type for 'sample'");
      goto cleanup;
    }
    else
      event->def->sample = (INT32_T)mxGetScalar(field) ;
  }
  
  fieldnumber = mxGetFieldNumber(prhs[0], "offset");
  if (fieldnumber<0) {
    mexErrMsgTxt("field 'offset' is missing");
    goto cleanup;
  }
  else
  {
    field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
    if (!mxIsNumeric(field) || mxIsEmpty(field)) {
      mexErrMsgTxt("invalid data type for 'offset'");
      goto cleanup;
    }
    else
      event->def->offset = (INT32_T)mxGetScalar(field) ;
  }
  
  fieldnumber = mxGetFieldNumber(prhs[0], "duration");
  if (fieldnumber<0) {
    mexErrMsgTxt("field 'duration' is missing");
    goto cleanup;
  }
  else
  {
    field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
    if (!mxIsNumeric(field) || mxIsEmpty(field)) {
      mexErrMsgTxt("invalid data type for 'duration'");
      goto cleanup;
    }
    else
      event->def->duration = (INT32_T)mxGetScalar(field) ;
  }
  
  fieldnumber = mxGetFieldNumber(prhs[0], "bufsize");
  if (fieldnumber<0) {
    mexErrMsgTxt("field 'bufsize' is missing");
    goto cleanup;
  }
  else
  {
    field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
    if (!mxIsNumeric(field) || mxIsEmpty(field)) {
      mexErrMsgTxt("invalid data type for 'bufsize'");
      goto cleanup;
    }
    else
      event->def->bufsize = (INT32_T)mxGetScalar(field) ;
  }
  
  fieldnumber = mxGetFieldNumber(prhs[0], "buf");
  if (fieldnumber<0) {
    mexErrMsgTxt("field 'buf' is missing");
    goto cleanup;
  }
  else
  {
    field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
    if (!mxIsUint8(field)) {
      mexErrMsgTxt("invalid data type for 'buf'");
      goto cleanup;
    }
    else if (mxGetNumberOfElements(field) != event->def->bufsize ) {
      mexErrMsgTxt("invalid number of elements (buf)");
      goto cleanup;
    }
    else {
      /* FIXME check the allocation */
      event->buf = malloc(event->def->bufsize);
      memcpy(event->buf, mxGetPr(field), event->def->bufsize);
    }
  }
  
  /* construct a PUT_EVT request */
  request->def->bufsize = append(&request->buf, request->def->bufsize, event->def, sizeof(eventdef_t));
  request->def->bufsize = append(&request->buf, request->def->bufsize, event->buf, event->def->bufsize);
  
  /* the event structure is not needed any more */
  FREE(event->def);
  FREE(event->buf);
  FREE(event);
  
  /* open the TCP socket */
  if ((server = open_connection(hostname, port)) < 0) {
    sprintf(msg, "ERROR: failed to create socket (%d)\n", server);
		mexErrMsgTxt(msg);
  }
  
  /* write the request, read the response */
  clientrequest(server, request, &response);
  close_connection(server);
  
  /* the request structure is not needed any more */
  if (request) {
    FREE(request->def);
    FREE(request->buf);
    FREE(request);
  }
  
  /* check that the response is PUT_OK */
  if (!response)
    mexErrMsgTxt("unknown error in response\n");
  else if (!response->def)
    mexErrMsgTxt("unknown error in response\n");
  else if (response->def->command!=PUT_OK)
  {
    sprintf(msg, "ERROR: the buffer returned an error (%d)\n", response->def->command);
		mexErrMsgTxt(msg);
  }
  
  /* the response structure is not needed any more */
  if (response) {
    FREE(response->def);
    FREE(response->buf);
    FREE(response);
  }
  
  return;
  
  cleanup:
    FREE(event->def);
    FREE(event->buf);
    FREE(event);
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

