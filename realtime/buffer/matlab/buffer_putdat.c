/*
 * Copyright (C) 2008, Robert Oostenveld
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 * $Log: buffer_putdat.c,v $
 * Revision 1.3  2009/06/30 11:35:48  roboos
 * pass more explicit error from buffer onto matlab
 *
 * Revision 1.2  2008/10/29 20:46:12  roboos
 * consistent use of open_connection and close_connection
 * there were some incorrect uses of close() that on windows did not actually close, resulting in the buffer running out of sockets/threads after prolonged use
 *
 * Revision 1.1  2008/10/22 10:46:29  roboos
 * created first implementation, based on puthdr
 *
 *
 */

#include "mex.h"
#include "matrix.h"
#include "buffer.h"

void buffer_putdat(char *hostname, int port, mxArray * plhs[], const mxArray * prhs[])
{
  size_t n;
  int server, fieldnumber;
  mxArray *field;
  char msg[512];
  
  message_t    *request  = NULL;
  message_t    *response = NULL;
  data_t       *data     = NULL;
  
    /* allocate the data  */
  data      = malloc(sizeof(data_t));
  data->def = malloc(sizeof(datadef_t));
  
  /* allocate the request message */
  request      = malloc(sizeof(message_t));
  request->def = malloc(sizeof(messagedef_t));
  request->buf = NULL;
  request->def->version = VERSION;
  request->def->command = PUT_DAT;
  request->def->bufsize = 0;
  
  /* define the data, it has the fields "nchans", "nsamples", "data_type" */
  
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
      data->def->nchans    = (UINT32_T)mxGetScalar(field) ;
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
      data->def->nsamples    = (UINT32_T)mxGetScalar(field) ;
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
      data->def->data_type    = (UINT32_T)mxGetScalar(field) ;
  }
  
  fieldnumber = mxGetFieldNumber(prhs[0], "buf");
  if (fieldnumber<0) {
    mexErrMsgTxt("field 'buf' is missing");
    goto cleanup;
  }
  else
  {
    field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
    if (!mxIsNumeric(field) || mxIsEmpty(field)) {
      mexErrMsgTxt("invalid data type for 'buf'");
      goto cleanup;
    }
    else
    {
      if (data->def->data_type == DATATYPE_CHAR)
        data->def->bufsize = data->def->nchans * data->def->nsamples * WORDSIZE_CHAR;
      else if (data->def->data_type == DATATYPE_UINT8)
        data->def->bufsize = data->def->nchans * data->def->nsamples * WORDSIZE_UINT8;
      else if (data->def->data_type == DATATYPE_UINT16)
        data->def->bufsize = data->def->nchans * data->def->nsamples * WORDSIZE_UINT16;
      else if (data->def->data_type == DATATYPE_UINT32)
        data->def->bufsize = data->def->nchans * data->def->nsamples * WORDSIZE_UINT32;
      else if (data->def->data_type == DATATYPE_UINT64)
        data->def->bufsize = data->def->nchans * data->def->nsamples * WORDSIZE_UINT64;
      else if (data->def->data_type == DATATYPE_INT8)
        data->def->bufsize = data->def->nchans * data->def->nsamples * WORDSIZE_INT8;
      else if (data->def->data_type == DATATYPE_INT16)
        data->def->bufsize = data->def->nchans * data->def->nsamples * WORDSIZE_INT16;
      else if (data->def->data_type == DATATYPE_INT32)
        data->def->bufsize = data->def->nchans * data->def->nsamples * WORDSIZE_INT32;
      else if (data->def->data_type == DATATYPE_INT64)
        data->def->bufsize = data->def->nchans * data->def->nsamples * WORDSIZE_INT64;
      else if (data->def->data_type == DATATYPE_FLOAT32)
        data->def->bufsize = data->def->nchans * data->def->nsamples * WORDSIZE_FLOAT32;
      else if (data->def->data_type == DATATYPE_FLOAT64)
        data->def->bufsize = data->def->nchans * data->def->nsamples * WORDSIZE_FLOAT64;
      data->buf = (void *)mxGetPr(field); /* directly point to the memory that is managed by Matlab */
    }
  }
  
  /* construct a PUT_DAT request */
  request->def->bufsize = append(&request->buf, request->def->bufsize, data->def, sizeof(datadef_t));
  request->def->bufsize = append(&request->buf, request->def->bufsize, data->buf, data->def->bufsize);
  
  /* the data structure is not needed any more */
  FREE(data->def);
  /* FREE(data->buf); this should not be freed, since it points to a piece of memory that is managed by Matlab */
  FREE(data);
  
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
    FREE(data->def);
  /* FREE(data->buf); this should not be freed, since it points to a piece of memory that is managed by Matlab */
    FREE(data);
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
