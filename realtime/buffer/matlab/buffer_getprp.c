/*
 * Copyright (C) 2008, Robert Oostenveld
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 * $Log: buffer_getprp.c,v $
 * Revision 1.5  2009/06/30 11:35:48  roboos
 * pass more explicit error from buffer onto matlab
 *
 * Revision 1.4  2008/10/29 20:46:12  roboos
 * consistent use of open_connection and close_connection
 * there were some incorrect uses of close() that on windows did not actually close, resulting in the buffer running out of sockets/threads after prolonged use
 *
 * Revision 1.3  2008/07/09 13:34:21  roboos
 * small change in verbose output, using verbose=0|1
 *
 * Revision 1.2  2008/05/22 09:27:22  roboos
 * fixed some issues with Borland compiler, correct pointer arithmetic, declarations at beginning
 *
 * Revision 1.1  2008/03/23 12:45:09  roboos
 * new implementation
 *
 *
 */

#include "mex.h"
#include "matrix.h"
#include "buffer.h"

#define NUMBER_OF_FIELDS 6

void buffer_getprp(char *hostname, int port, mxArray *plhs[], const mxArray *prhs[])
{
  int server;
  int verbose = 0;
  int i, nproperties;
  int offset;
  double *val;
  char msg[512];
  
  mxArray *bufptr;
  
  message_t    *request     = NULL;
  message_t    *response    = NULL;
  property_t   *property    = NULL;
  
  /* this is for the Matlab specific output */
  const char *field_names[] = {
    "type_type",
    "type_numel",
    "value_type",
    "value_numel",
    "bufsize",
    "buf"
  };
  
    /* allocate the elements that will be used in the communication */
  request      = malloc(sizeof(message_t));
  request->def = malloc(sizeof(messagedef_t));
  request->buf = NULL;
  request->def->version = VERSION;
  request->def->command = GET_PRP;
  request->def->bufsize = 0;
  
  #ifdef NONSENSE
  /* properties are selected based on their type and not on number */
  /* this is realised at the moment by sending a propertysel that is equal to a propertydef+buf with value being ignored */
  if ((prhs[0]!=NULL) && (mxGetNumberOfElements(prhs[0])==2) && (mxIsDouble(prhs[0])) && (!mxIsComplex(prhs[0]))) {
    /* fprintf(stderr, "args OK\n"); */
    val = (double *)mxGetData(prhs[0]);
    propertysel.begproperty = (UINT32_T)(val[0]);
    propertysel.endproperty = (UINT32_T)(val[1]);
    print_propertysel(&propertysel);
    request->def->bufsize = append(&request->buf, request->def->bufsize, &propertysel, sizeof(propertysel_t));
  }
  #endif
  
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
    property = malloc(sizeof(property_t));
    
    /* first count the number of propertys */
    nproperties = 0;
    offset = 0;
    while (offset<response->def->bufsize) {
      property->def = (char *)response->buf + offset;
      property->buf = (char *)response->buf + offset + sizeof(propertydef_t);
      print_propertydef(property->def);
      offset += sizeof(propertydef_t) + property->def->bufsize;
      nproperties++;
    }
    
    /* create a structure array that can hold all propertys */
    plhs[0] = mxCreateStructMatrix(1, nproperties, NUMBER_OF_FIELDS, field_names);
    
    offset = 0;
    for (i=0; i<nproperties; i++) {
      property->def = (char *)response->buf + offset;
      property->buf = (char *)response->buf + offset + sizeof(propertydef_t);
      
      bufptr = mxCreateNumericMatrix(1, property->def->bufsize, mxUINT8_CLASS, mxREAL);
      memcpy(mxGetPr(bufptr), property->buf, property->def->bufsize);
      
      mxSetFieldByNumber(plhs[0], i, 0, mxCreateDoubleScalar((double)property->def->type_type));
      mxSetFieldByNumber(plhs[0], i, 1, mxCreateDoubleScalar((double)property->def->type_numel));
      mxSetFieldByNumber(plhs[0], i, 2, mxCreateDoubleScalar((double)property->def->value_type));
      mxSetFieldByNumber(plhs[0], i, 3, mxCreateDoubleScalar((double)property->def->value_numel));
      mxSetFieldByNumber(plhs[0], i, 4, mxCreateDoubleScalar((double)property->def->bufsize));
      mxSetFieldByNumber(plhs[0], i, 5, bufptr);
      offset += sizeof(propertydef_t) + property->def->bufsize;
    }
    FREE(property);
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

