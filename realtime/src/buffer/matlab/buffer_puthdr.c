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

#define SIZE_NIFTI_1   348

ft_chunk_t *encodeChannelNames(const mxArray *L, int N) {
  ft_chunk_t *chunk = NULL;
  char *ptr;
  int i,size=0,NL;

  if (!mxIsCell(L)) return NULL;

  NL = mxGetNumberOfElements(L);
  if (NL>N) return NULL;

  for (i=0;i<NL;i++) {
    int li = 0;
    mxArray *S = mxGetCell(L, i);
    if (S == NULL || !mxIsChar(S)) return NULL;
    li = mxGetNumberOfElements(S);
    size += li+1; /* add 1 for trailing 0 */
  }
  /* L contains less names than we've got channels - this is
     for the 0's for the anonymous channels */
  size += N-NL; 

  /* now we (over)allocate the chunk */
  chunk = (ft_chunk_t *) mxMalloc(size + sizeof(ft_chunkdef_t));

  chunk->def.size = size;
  chunk->def.type = FT_CHUNK_CHANNEL_NAMES;

  ptr = chunk->data;

  /* run through the cell-array again */
  for (i=0;i<NL;i++) {
    mxArray *S = mxGetCell(L, i);
    int n = mxGetNumberOfElements(S) + 1;
    mxGetString(S, ptr, n);
    ptr+=n;
  }
  for (i=NL;i<N;i++) {
    *ptr++ = 0;
  }

  /* that's it - the receiver should free this thing using mxFree */
  return chunk;
}


ft_chunk_t *encodeResolutions(const mxArray *R, int N) {
  ft_chunk_t *chunk = NULL;

  /* possible extension: allow other types than 'double' and convert below */
  if (!mxIsDouble(R) || mxIsComplex(R)) return NULL;
  if (!(mxGetM(R)==N && mxGetN(R)==1) && (mxGetM(R)==1 && mxGetN(R)==N)) return NULL;

  /* now we (over)allocate the chunk */
  chunk = (ft_chunk_t *) mxMalloc(N*sizeof(double) + sizeof(ft_chunkdef_t));

  chunk->def.size = N*sizeof(double);
  chunk->def.type = FT_CHUNK_RESOLUTIONS;
  memcpy(chunk->data, mxGetPr(R), N*sizeof(double));

  /* that's it - the receiver should free this thing using mxFree */
  return chunk;
}

int buffer_puthdr(int server, mxArray * plhs[], const mxArray * prhs[])
{
  int fieldnumber;
  mxArray *field;
  int result;

  message_t     request;
  messagedef_t  request_def;
  message_t    *response = NULL;
  headerdef_t   header_def;

  ft_chunkdef_t chunk_def;

  /* allocate the request message */
  request.def = &request_def;
  request.buf = NULL;
  request_def.version = VERSION;
  request_def.command = PUT_HDR;
  request_def.bufsize = 0;

  /* define the header, it has the fields "nchans", "nsamples", "nevents", "fsample", "data_type" */
  if (mxGetNumberOfElements(prhs[0])!=1)
    mexErrMsgTxt("Only one header can be put into the buffer at a time.");

  fieldnumber = mxGetFieldNumber(prhs[0], "nchans");
  if (fieldnumber<0) 
    mexErrMsgTxt("field 'nchans' is missing");
  field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
  if (!mxIsNumeric(field) || mxIsEmpty(field)) 
    mexErrMsgTxt("invalid data type for 'nchans'");
  header_def.nchans    = (UINT32_T)mxGetScalar(field) ;

  fieldnumber = mxGetFieldNumber(prhs[0], "nsamples");
  if (fieldnumber<0) 
    mexErrMsgTxt("field 'nsamples' is missing");
  field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
  if (!mxIsNumeric(field) || mxIsEmpty(field)) 
    mexErrMsgTxt("invalid data type for 'nsamples'");
  header_def.nsamples    = (UINT32_T)mxGetScalar(field) ;

  fieldnumber = mxGetFieldNumber(prhs[0], "nevents");
  if (fieldnumber<0) 
    mexErrMsgTxt("field is missing 'nevents'");
  field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
  if (!mxIsNumeric(field) || mxIsEmpty(field)) 
    mexErrMsgTxt("invalid data type for 'nevents'");
  header_def.nevents    = (UINT32_T)mxGetScalar(field) ;

  fieldnumber = mxGetFieldNumber(prhs[0], "fsample");
  if (fieldnumber<0) 
    mexErrMsgTxt("field is missing 'fsample'");
  field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
  if (!mxIsNumeric(field) || mxIsEmpty(field)) 
    mexErrMsgTxt("invalid data type for 'fsample'");
  header_def.fsample    = (float)mxGetScalar(field) ;

  fieldnumber = mxGetFieldNumber(prhs[0], "data_type");
  if (fieldnumber<0) 
    mexErrMsgTxt("field 'data_type' is missing");
  field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
  if (!mxIsNumeric(field) || mxIsEmpty(field)) 
    mexErrMsgTxt("invalid data type for 'data_type'");
  header_def.data_type    = (UINT32_T)mxGetScalar(field) ;

  /* construct a PUT_HDR request */
  request_def.bufsize = ft_mx_append(&request.buf, request_def.bufsize, &header_def, sizeof(headerdef_t));

  /* append existing chunks to request.buf, set correct header_def.bufsize at the end */
  fieldnumber = mxGetFieldNumber(prhs[0], "nifti_1");
  if (fieldnumber>=0) {
    field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
    if (!mxIsUint8(field) || mxGetNumberOfElements(field)!=SIZE_NIFTI_1) {
      mexWarnMsgTxt("invalid data type for field 'nifti_1' -- ignoring");
    } else {
      chunk_def.size = SIZE_NIFTI_1;
      chunk_def.type = FT_CHUNK_NIFTI1;
      request_def.bufsize = ft_mx_append(&request.buf, request_def.bufsize, &chunk_def, sizeof(chunk_def));
      request_def.bufsize = ft_mx_append(&request.buf, request_def.bufsize, mxGetData(field), chunk_def.size);
    }
  }

  fieldnumber = mxGetFieldNumber(prhs[0], "siemensap");
  if (fieldnumber>=0) {
    field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
    if (!mxIsUint8(field)) {
      mexWarnMsgTxt("invalid data type for field 'siemensap' -- ignoring");
    } else {
      chunk_def.size = mxGetNumberOfElements(field);
      chunk_def.type = FT_CHUNK_SIEMENS_AP;
      request_def.bufsize = ft_mx_append(&request.buf, request_def.bufsize, &chunk_def, sizeof(chunk_def));
      request_def.bufsize = ft_mx_append(&request.buf, request_def.bufsize, mxGetData(field), chunk_def.size);
    }
  }

  fieldnumber = mxGetFieldNumber(prhs[0], "ctf_res4");
  if (fieldnumber>=0) {
    field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
    if (!mxIsUint8(field)) {
      mexWarnMsgTxt("invalid data type for field 'ctf_res4' -- ignoring");
    } else {
      chunk_def.size = mxGetNumberOfElements(field);
      chunk_def.type = FT_CHUNK_CTF_RES4;
      request_def.bufsize = ft_mx_append(&request.buf, request_def.bufsize, &chunk_def, sizeof(chunk_def));
      request_def.bufsize = ft_mx_append(&request.buf, request_def.bufsize, mxGetData(field), chunk_def.size);
    }
  }	

  fieldnumber = mxGetFieldNumber(prhs[0], "neuromag_header");
  if (fieldnumber>=0) {
    field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
    if (!mxIsUint8(field)) {
      mexWarnMsgTxt("invalid data type for field 'neuromag_header' -- ignoring");
    } else {
      chunk_def.size = mxGetNumberOfElements(field);
      chunk_def.type = FT_CHUNK_NEUROMAG_HEADER;
      request_def.bufsize = ft_mx_append(&request.buf, request_def.bufsize, &chunk_def, sizeof(chunk_def));
      request_def.bufsize = ft_mx_append(&request.buf, request_def.bufsize, mxGetData(field), chunk_def.size);
    }
  }	

  fieldnumber = mxGetFieldNumber(prhs[0], "neuromag_isotrak");
  if (fieldnumber>=0) {
    field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
    if (!mxIsUint8(field)) {
      mexWarnMsgTxt("invalid data type for field 'neuromag_isotrak' -- ignoring");
    } else {
      chunk_def.size = mxGetNumberOfElements(field);
      chunk_def.type = FT_CHUNK_NEUROMAG_ISOTRAK;
      request_def.bufsize = ft_mx_append(&request.buf, request_def.bufsize, &chunk_def, sizeof(chunk_def));
      request_def.bufsize = ft_mx_append(&request.buf, request_def.bufsize, mxGetData(field), chunk_def.size);
    }
  }	

  fieldnumber = mxGetFieldNumber(prhs[0], "neuromag_hpiresult");
  if (fieldnumber>=0) {
    field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
    if (!mxIsUint8(field)) {
      mexWarnMsgTxt("invalid data type for field 'neuromag_hpiresult' -- ignoring");
    } else {
      chunk_def.size = mxGetNumberOfElements(field);
      chunk_def.type = FT_CHUNK_NEUROMAG_HPIRESULT;
      request_def.bufsize = ft_mx_append(&request.buf, request_def.bufsize, &chunk_def, sizeof(chunk_def));
      request_def.bufsize = ft_mx_append(&request.buf, request_def.bufsize, mxGetData(field), chunk_def.size);
    }
  }	

  fieldnumber = mxGetFieldNumber(prhs[0], "channel_names");
  if (fieldnumber>=0) {
    ft_chunk_t *chunk = NULL;
    field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
    chunk = encodeChannelNames(field, header_def.nchans);
    if (chunk == NULL) {
      mexWarnMsgTxt("invalid data type for field 'channel_names' -- ignoring.");
    } else {
      request_def.bufsize = ft_mx_append(&request.buf, request_def.bufsize, chunk, sizeof(ft_chunkdef_t) + chunk->def.size);
      mxFree(chunk);
    }
  }

  fieldnumber = mxGetFieldNumber(prhs[0], "resolutions");
  if (fieldnumber>=0) {
    ft_chunk_t *chunk = NULL;
    field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
    chunk = encodeResolutions(field, header_def.nchans);
    if (chunk == NULL) {
      mexWarnMsgTxt("invalid data type for field 'resolutions' -- ignoring.");
    } else {
      request_def.bufsize = ft_mx_append(&request.buf, request_def.bufsize, chunk, sizeof(ft_chunkdef_t) + chunk->def.size);
      mxFree(chunk);
    }
  }		

  /* header->def->bufsize is the request->def->bufsize - sizeof(header->def) */
  ((headerdef_t *) request.buf)->bufsize = request_def.bufsize - sizeof(headerdef_t);

  /* write the request, read the response */
  result = clientrequest(server, &request, &response);

  /* the request structure is not needed any more, but only .buf needs to be free'd */
  if (request.buf != NULL) mxFree(request.buf);

  if (result == 0) {
    /* check that the response is PUT_OK */
    if (!response)
      mexErrMsgTxt("unknown error in response\n");
    if (!response->def) {
      FREE(response->buf);
      FREE(response);
      mexErrMsgTxt("unknown error in response\n");
    }
    if (response->def->command!=PUT_OK) {
      result = response->def->command;
    }
  }
  /* the response structure is not needed any more */
  if (response) {
    FREE(response->def);
    FREE(response->buf);
    FREE(response);
  }
  return result;
}
