/*
 * Copyright (C) 2008, Robert Oostenveld
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 * $Log: buffer_putprp.c,v $
 * Revision 1.3  2009/06/30 11:35:48  roboos
 * pass more explicit error from buffer onto matlab
 *
 * Revision 1.2  2008/10/29 20:46:12  roboos
 * consistent use of open_connection and close_connection
 * there were some incorrect uses of close() that on windows did not actually close, resulting in the buffer running out of sockets/threads after prolonged use
 *
 * Revision 1.1  2008/03/23 12:45:09  roboos
 * new implementation
 *
 *
 */

#include "mex.h"
#include "matrix.h"
#include "buffer.h"

void buffer_putprp(char *hostname, int port, mxArray * plhs[], const mxArray * prhs[])
{
  size_t n;
  int server, fieldnumber, bufsize;
  mxArray *field;
  char msg[512];
  
  message_t   *request  = NULL;
  message_t   *response = NULL;
  property_t  *property = NULL;
  
  /* allocate the property */
  property      = malloc(sizeof(property_t));
  property->def = malloc(sizeof(propertydef_t));
  property->buf = NULL;
  property->def->bufsize   = 0;
  
  /* allocate the request message */
  request      = malloc(sizeof(message_t));
  request->def = malloc(sizeof(messagedef_t));
  request->buf = NULL;
  request->def->version = VERSION;
  request->def->command = PUT_PRP;
  request->def->bufsize = 0;
  
  /* define the property, it has the fields type_type type_numel value_type value_numel */
  
  /* FIXME loop over mutiple properties in case of a property-array */
  
  fieldnumber = mxGetFieldNumber(prhs[0], "type_type");
  if (fieldnumber<0) {
    mexErrMsgTxt("field is missing");
    goto cleanup;
  }
  else
  {
    field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
    if (!mxIsNumeric(field) || mxIsEmpty(field)) {
      mexErrMsgTxt("invalid data type");
      goto cleanup;
    }
    else
      property->def->type_type = (UINT32_T)mxGetScalar(field) ;
  }
  
  fieldnumber = mxGetFieldNumber(prhs[0], "type_numel");
  if (fieldnumber<0) {
    mexErrMsgTxt("field is missing");
    goto cleanup;
  }
  else
  {
    field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
    if (!mxIsNumeric(field) || mxIsEmpty(field)) {
      mexErrMsgTxt("invalid data type");
      goto cleanup;
    }
    else
      property->def->type_numel = (UINT32_T)mxGetScalar(field) ;
  }
  
  fieldnumber = mxGetFieldNumber(prhs[0], "value_type");
  if (fieldnumber<0) {
    mexErrMsgTxt("field is missing");
    goto cleanup;
  }
  else
  {
    field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
    if (!mxIsNumeric(field) || mxIsEmpty(field)) {
      mexErrMsgTxt("invalid data type");
      goto cleanup;
    }
    else
      property->def->value_type = (UINT32_T)mxGetScalar(field) ;
  }
  
  fieldnumber = mxGetFieldNumber(prhs[0], "value_numel");
  if (fieldnumber<0) {
    mexErrMsgTxt("field is missing");
    goto cleanup;
  }
  else
  {
    field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
    if (!mxIsNumeric(field) || mxIsEmpty(field)) {
      mexErrMsgTxt("invalid data type");
      goto cleanup;
    }
    else
      property->def->value_numel = (UINT32_T)mxGetScalar(field) ;
  }
  
  fieldnumber = mxGetFieldNumber(prhs[0], "bufsize");
  if (fieldnumber<0) {
    mexErrMsgTxt("field is missing");
    goto cleanup;
  }
  else
  {
    field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
    if (!mxIsNumeric(field) || mxIsEmpty(field)) {
      mexErrMsgTxt("invalid data type");
      goto cleanup;
    }
    else
      bufsize = (UINT32_T)mxGetScalar(field) ;
  }
  
  fieldnumber = mxGetFieldNumber(prhs[0], "buf");
  if (fieldnumber<0) {
    mexErrMsgTxt("field is missing");
    goto cleanup;
  }
  else
  {
    field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
    if (!mxIsUint8(field) || mxIsComplex(field)) {
      mexErrMsgTxt("buffer should be real-valued uint8");
      goto cleanup;
    }
    if (mxGetNumberOfElements(field)!=bufsize) {
      mexErrMsgTxt("length of buffer does not correspond to bufsize");
      goto cleanup;
    }
    else
      property->def->bufsize = append(&property->buf, property->def->bufsize, mxGetData(field), bufsize);
  }
  
  /* construct a PUT_PRP request */
  request->def->bufsize = append(&request->buf, request->def->bufsize, property->def, sizeof(propertydef_t));
  request->def->bufsize = append(&request->buf, request->def->bufsize, property->buf, property->def->bufsize);
  
  /* the property structure is not needed any more */
  FREE(property->def);
  FREE(property->buf);
  FREE(property);
  
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
    FREE(property->def);
    FREE(property->buf);
    FREE(property);
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

