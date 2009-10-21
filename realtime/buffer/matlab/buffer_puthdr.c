/*
 * Copyright (C) 2008, Robert Oostenveld
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 * $Log: buffer_puthdr.c,v $
 * Revision 1.6  2009/06/30 11:35:48  roboos
 * pass more explicit error from buffer onto matlab
 *
 * Revision 1.5  2008/10/29 20:46:12  roboos
 * consistent use of open_connection and close_connection
 * there were some incorrect uses of close() that on windows did not actually close, resulting in the buffer running out of sockets/threads after prolonged use
 *
 * Revision 1.4  2008/06/19 21:00:36  roboos
 * give more information in the error output messages
 *
 * Revision 1.3  2008/03/23 12:51:23  roboos
 * moved command=xxx in code to better location, no functional place
 *
 * Revision 1.2  2008/03/09 22:36:38  roboos
 * implemented new buffer mex-function, which acts as dispatcher for all functionality
 * changed buffer_xxx code from mex-function into c-function to be included in dispatcher
 *
 * Revision 1.1  2008/03/07 14:46:06  roboos
 * new implementation, the events are not yet fully complete but the basic functionality has been tested
 *
 *
 */

#include "mex.h"
#include "matrix.h"
#include "buffer.h"

void buffer_puthdr(char *hostname, int port, mxArray * plhs[], const mxArray * prhs[])
{
  size_t n;
  int server, fieldnumber;
  mxArray *field;
  char msg[512];
  
  message_t    *request  = NULL;
  message_t    *response = NULL;
  header_t     *header   = NULL;
  
  /* allocate the header */
  header      = malloc(sizeof(header_t));
  header->def = malloc(sizeof(headerdef_t));
  header->buf = NULL;
  header->def->bufsize   = 0;
  
  /* allocate the request message */
  request      = malloc(sizeof(message_t));
  request->def = malloc(sizeof(messagedef_t));
  request->buf = NULL;
  request->def->version = VERSION;
  request->def->command = PUT_HDR;
  request->def->bufsize = 0;
  
  /* define the header, it has the fields "nchans", "nsamples", "nevents", "fsample", "data_type" */
  
  fieldnumber = mxGetFieldNumber(prhs[0], "nchans");
  if (fieldnumber<0) {
    mexErrMsgTxt("field 'nchans' is missing");
    goto cleanup; /* FIXME will not be reached */
  }
  else
  {
    field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
    if (!mxIsNumeric(field) || mxIsEmpty(field)) {
      mexErrMsgTxt("invalid data type for 'nchans'");
      goto cleanup;
    }
    else
      header->def->nchans    = (UINT32_T)mxGetScalar(field) ;
  }
  
  fieldnumber = mxGetFieldNumber(prhs[0], "nsamples");
  if (fieldnumber<0) {
    mexErrMsgTxt("field 'nsamples' is missing");
    goto cleanup;
  }
  else
  {
    field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
    if (!mxIsNumeric(field) || mxIsEmpty(field)) {
      mexErrMsgTxt("invalid data type for 'nsamples'");
      goto cleanup;
    }
    else
      header->def->nsamples    = (UINT32_T)mxGetScalar(field) ;
  }
  
  fieldnumber = mxGetFieldNumber(prhs[0], "nevents");
  if (fieldnumber<0) {
    mexErrMsgTxt("field is missing 'nevents'");
    goto cleanup;
  }
  else
  {
    field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
    if (!mxIsNumeric(field) || mxIsEmpty(field)) {
      mexErrMsgTxt("invalid data type for 'nevents'");
      goto cleanup;
    }
    else
      header->def->nevents    = (UINT32_T)mxGetScalar(field) ;
  }
  
  fieldnumber = mxGetFieldNumber(prhs[0], "fsample");
  if (fieldnumber<0) {
    mexErrMsgTxt("field is missing 'fsample'");
    goto cleanup;
  }
  else
  {
    field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
    if (!mxIsNumeric(field) || mxIsEmpty(field)) {
      mexErrMsgTxt("invalid data type for 'fsample'");
      goto cleanup;
    }
    else
      header->def->fsample    = (float)mxGetScalar(field) ;
  }
  
  fieldnumber = mxGetFieldNumber(prhs[0], "data_type");
  if (fieldnumber<0) {
    mexErrMsgTxt("field 'data_type' is missing");
    goto cleanup;
  }
  else
  {
    field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
    if (!mxIsNumeric(field) || mxIsEmpty(field)) {
      mexErrMsgTxt("invalid data type for 'data_type'");
      goto cleanup;
    }
    else
      header->def->data_type    = (UINT32_T)mxGetScalar(field) ;
  }
  
  /* construct a PUT_HDR request */
  request->def->bufsize = append(&request->buf, request->def->bufsize, header->def, sizeof(headerdef_t));
  request->def->bufsize = append(&request->buf, request->def->bufsize, header->buf, header->def->bufsize);
  
  /* the header structure is not needed any more */
  FREE(header->def);
  FREE(header->buf);
  FREE(header);
  
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
    FREE(header->def);
    FREE(header->buf);
    FREE(header);
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

