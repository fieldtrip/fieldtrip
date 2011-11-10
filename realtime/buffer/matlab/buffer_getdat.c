/*
 * Copyright (C) 2008, Robert Oostenveld
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 * $Id$
 */

#include "mex.h"
#include "matrix.h"
#include "buffer.h"

#define NUMBER_OF_FIELDS 5

int buffer_getdat(int server, mxArray *plhs[], const mxArray *prhs[])
{
  int verbose = 0;
  double *val;
  int result = 0;
  
  message_t *request  = NULL;
  message_t *response = NULL;
  datasel_t datasel;
  
  /* this is for the Matlab specific output */
  const char *field_names[NUMBER_OF_FIELDS] = {"nchans", "nsamples", "data_type", "bufsize", "buf"};
  
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
   
  if (verbose) print_request(request->def);
  result = clientrequest(server, request, &response);
  if (verbose) print_response(response->def);
  
  if (result == 0) {
    if (response->def->command==GET_OK) {
      mxArray *datp       = NULL;
      datadef_t *data_def = (datadef_t *) response->buf;
      void *data_buf      = (void *)((char *)response->buf + sizeof(datadef_t));
      
      if (verbose) print_datadef(data_def);

      switch (data_def->data_type) {
        case DATATYPE_UINT8:
        datp = mxCreateNumericMatrix(data_def->nchans, data_def->nsamples, mxUINT8_CLASS, mxREAL);
        memcpy(mxGetPr(datp), data_buf, data_def->nchans*data_def->nsamples*WORDSIZE_UINT8);
        break;
        case DATATYPE_UINT16:
        datp = mxCreateNumericMatrix(data_def->nchans, data_def->nsamples, mxUINT16_CLASS, mxREAL);
        memcpy(mxGetPr(datp), data_buf, data_def->nchans*data_def->nsamples*WORDSIZE_UINT16);
        break;
        case DATATYPE_UINT32:
        datp = mxCreateNumericMatrix(data_def->nchans, data_def->nsamples, mxUINT32_CLASS, mxREAL);
        memcpy(mxGetPr(datp), data_buf, data_def->nchans*data_def->nsamples*WORDSIZE_UINT32);
        break;
        case DATATYPE_UINT64:
        datp = mxCreateNumericMatrix(data_def->nchans, data_def->nsamples, mxUINT64_CLASS, mxREAL);
        memcpy(mxGetPr(datp), data_buf, data_def->nchans*data_def->nsamples*WORDSIZE_UINT64);
        break;	  
        case DATATYPE_INT8:
        datp = mxCreateNumericMatrix(data_def->nchans, data_def->nsamples, mxINT8_CLASS, mxREAL);
        memcpy(mxGetPr(datp), data_buf, data_def->nchans*data_def->nsamples*WORDSIZE_INT8);
        break;
        case DATATYPE_INT16:
        datp = mxCreateNumericMatrix(data_def->nchans, data_def->nsamples, mxINT16_CLASS, mxREAL);
        memcpy(mxGetPr(datp), data_buf, data_def->nchans*data_def->nsamples*WORDSIZE_INT16);
        break;
        case DATATYPE_INT32:
        datp = mxCreateNumericMatrix(data_def->nchans, data_def->nsamples, mxINT32_CLASS, mxREAL);
        memcpy(mxGetPr(datp), data_buf, data_def->nchans*data_def->nsamples*WORDSIZE_INT32);
        break;
        case DATATYPE_INT64:
        datp = mxCreateNumericMatrix(data_def->nchans, data_def->nsamples, mxINT64_CLASS, mxREAL);
        memcpy(mxGetPr(datp), data_buf, data_def->nchans*data_def->nsamples*WORDSIZE_INT64);
        break;
        case DATATYPE_FLOAT32:
        datp = mxCreateNumericMatrix(data_def->nchans, data_def->nsamples, mxSINGLE_CLASS, mxREAL);
        memcpy(mxGetPr(datp), data_buf, data_def->nchans*data_def->nsamples*WORDSIZE_FLOAT32);
        break;
        case DATATYPE_FLOAT64:
        datp = mxCreateNumericMatrix(data_def->nchans, data_def->nsamples, mxDOUBLE_CLASS, mxREAL);
        memcpy(mxGetPr(datp), data_buf, data_def->nchans*data_def->nsamples*WORDSIZE_FLOAT64);
        break;
        default:
          result = -4;  /* mexErrMsgTxt("ERROR; unsupported data type\n"); */
          goto cleanup;
      }

      plhs[0] = mxCreateStructMatrix(1, 1, NUMBER_OF_FIELDS, field_names);
      mxSetFieldByNumber(plhs[0], 0, 0, mxCreateDoubleScalar((double)data_def->nchans));
      mxSetFieldByNumber(plhs[0], 0, 1, mxCreateDoubleScalar((double)(data_def->nsamples)));
      mxSetFieldByNumber(plhs[0], 0, 2, mxCreateDoubleScalar((double)(data_def->data_type)));
      mxSetFieldByNumber(plhs[0], 0, 3, mxCreateDoubleScalar((double)(data_def->bufsize)));
      mxSetFieldByNumber(plhs[0], 0, 4, datp);
    }
    else {
      result = response->def->command;
    }
  }
	 
cleanup:   
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

