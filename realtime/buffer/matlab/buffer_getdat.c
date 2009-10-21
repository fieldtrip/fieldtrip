/*
 * Copyright (C) 2008, Robert Oostenveld
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 * $Log: buffer_getdat.c,v $
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
 * Revision 1.7  2008/05/22 09:27:21  roboos
 * fixed some issues with Borland compiler, correct pointer arithmetic, declarations at beginning
 *
 * Revision 1.6  2008/03/23 12:47:17  roboos
 * implemented selection of samples
 *
 * Revision 1.5  2008/03/09 22:36:38  roboos
 * implemented new buffer mex-function, which acts as dispatcher for all functionality
 * changed buffer_xxx code from mex-function into c-function to be included in dispatcher
 *
 * Revision 1.4  2008/02/26 21:43:24  roboos
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

#define NUMBER_OF_FIELDS 5

void buffer_getdat(char *hostname, int port, mxArray *plhs[], const mxArray *prhs[])
{
  int server;
  int verbose = 0;
  size_t n;
  double *val;
  char msg[512];
  
  message_t *request  = NULL;
  message_t *response = NULL;
  header_t  *header   = NULL;
  data_t    *data     = NULL;
  datasel_t datasel;
  
  /* this is for the Matlab specific output */
  const char *field_names[NUMBER_OF_FIELDS] = {"nchans", "nsamples", "data_type", "bufsize", "buf"};
  mxArray *datp;
  
  /* allocate the elements that will be used in the communication */
  request      = malloc(sizeof(message_t));
  request->def = malloc(sizeof(messagedef_t));
  request->buf = NULL;
  request->def->version = VERSION;
  request->def->command = GET_DAT;
  request->def->bufsize = 0;
  
  if ((prhs[0]!=NULL) && (mxGetNumberOfElements(prhs[0])==2) && (mxIsDouble(prhs[0])) && (!mxIsComplex(prhs[0]))) {
    /* fprintf(stderr, "args OK\n"); */
    val = (double *)mxGetData(prhs[0]);
    datasel.begsample = (UINT32_T)(val[0]);
    datasel.endsample = (UINT32_T)(val[1]);
    if (verbose) print_datasel(&datasel);
    request->def->bufsize = append(&request->buf, request->def->bufsize, &datasel, sizeof(datasel_t));
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
    data      = malloc(sizeof(data_t));
    data->def = response->buf;
    data->buf = (char *)response->buf + sizeof(datadef_t);
    if (verbose) print_datadef(data->def);
    
    switch (data->def->data_type) {
      case DATATYPE_INT8:
        datp = mxCreateNumericMatrix(data->def->nchans, data->def->nsamples, mxINT8_CLASS, mxREAL);
        memcpy(mxGetPr(datp), data->buf, data->def->nchans*data->def->nsamples*WORDSIZE_INT8);
        break;
      case DATATYPE_INT16:
        datp = mxCreateNumericMatrix(data->def->nchans, data->def->nsamples, mxINT16_CLASS, mxREAL);
        memcpy(mxGetPr(datp), data->buf, data->def->nchans*data->def->nsamples*WORDSIZE_INT16);
        break;
      case DATATYPE_INT32:
        datp = mxCreateNumericMatrix(data->def->nchans, data->def->nsamples, mxINT32_CLASS, mxREAL);
        memcpy(mxGetPr(datp), data->buf, data->def->nchans*data->def->nsamples*WORDSIZE_INT32);
        break;
      case DATATYPE_INT64:
        datp = mxCreateNumericMatrix(data->def->nchans, data->def->nsamples, mxINT64_CLASS, mxREAL);
        memcpy(mxGetPr(datp), data->buf, data->def->nchans*data->def->nsamples*WORDSIZE_INT64);
        break;
      case DATATYPE_FLOAT32:
        datp = mxCreateNumericMatrix(data->def->nchans, data->def->nsamples, mxSINGLE_CLASS, mxREAL);
        memcpy(mxGetPr(datp), data->buf, data->def->nchans*data->def->nsamples*WORDSIZE_FLOAT32);
        break;
      case DATATYPE_FLOAT64:
        datp = mxCreateNumericMatrix(data->def->nchans, data->def->nsamples, mxDOUBLE_CLASS, mxREAL);
        memcpy(mxGetPr(datp), data->buf, data->def->nchans*data->def->nsamples*WORDSIZE_FLOAT64);
        break;
        default:
          mexErrMsgTxt("ERROR; unsupported data type\n");
    }
    
    plhs[0] = mxCreateStructMatrix(1, 1, NUMBER_OF_FIELDS, field_names);
    mxSetFieldByNumber(plhs[0], 0, 0, mxCreateDoubleScalar((double)data->def->nchans));
    mxSetFieldByNumber(plhs[0], 0, 1, mxCreateDoubleScalar((double)(data->def->nsamples)));
    mxSetFieldByNumber(plhs[0], 0, 2, mxCreateDoubleScalar((double)(data->def->data_type)));
    mxSetFieldByNumber(plhs[0], 0, 3, mxCreateDoubleScalar((double)(data->def->bufsize)));
    mxSetFieldByNumber(plhs[0], 0, 4, datp);
    FREE(data);
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

