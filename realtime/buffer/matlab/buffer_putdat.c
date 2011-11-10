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

int buffer_putdat(int server, mxArray * plhs[], const mxArray * prhs[])
{
	int fieldnumber;
	mxArray *field;
	int result;
  
	messagedef_t  request_def;
	message_t     request;
	datadef_t     data_def;
	data_t        data;  
	message_t    *response = NULL;

	data.def = &data_def;
	request.def = &request_def;
	request.buf = NULL;
	request_def.version = VERSION;	/* this is the same as request.def->version etc.*/
	request_def.command = PUT_DAT;
	request_def.bufsize = 0;
  
  /* define the data, it has the fields "nchans", "nsamples", "data_type" */
  
	fieldnumber = mxGetFieldNumber(prhs[0], "nchans");
	if (fieldnumber<0) 
		mexErrMsgTxt("field 'nchans' is missing");

	field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
	if (!mxIsNumeric(field) || mxIsEmpty(field)) 
		mexErrMsgTxt("invalid data type for 'nchans'");
    data_def.nchans = (UINT32_T)mxGetScalar(field) ;

	fieldnumber = mxGetFieldNumber(prhs[0], "nsamples");
	if (fieldnumber<0) 
		mexErrMsgTxt("field 'nsamples' is missing");
	field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
	if (!mxIsNumeric(field) || mxIsEmpty(field)) 
		mexErrMsgTxt("invalid data type for 'nsamples'");
	data_def.nsamples = (UINT32_T)mxGetScalar(field) ;
  
	fieldnumber = mxGetFieldNumber(prhs[0], "data_type");
	if (fieldnumber<0)
		mexErrMsgTxt("field 'data_type' is missing");

	field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
    if (!mxIsNumeric(field) || mxIsEmpty(field)) 
		mexErrMsgTxt("invalid data type for 'data_type'");
	data_def.data_type = (UINT32_T)mxGetScalar(field) ;
  
	fieldnumber = mxGetFieldNumber(prhs[0], "buf");
	if (fieldnumber<0)
		mexErrMsgTxt("field 'buf' is missing");
	field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
	if (!mxIsNumeric(field) || mxIsEmpty(field))
		mexErrMsgTxt("invalid data type for 'buf'");

	switch(data_def.data_type) {
		case DATATYPE_CHAR:
			data_def.bufsize = data_def.nchans * data_def.nsamples * WORDSIZE_CHAR;
			break;
		case DATATYPE_UINT8:
			data_def.bufsize = data_def.nchans * data_def.nsamples * WORDSIZE_UINT8;
			break;
		case DATATYPE_UINT16:
			data_def.bufsize = data_def.nchans * data_def.nsamples * WORDSIZE_UINT16;
			break;
		case DATATYPE_UINT32:
			data_def.bufsize = data_def.nchans * data_def.nsamples * WORDSIZE_UINT32;
			break;
		case DATATYPE_UINT64:
			data_def.bufsize = data_def.nchans * data_def.nsamples * WORDSIZE_UINT64;
			break;
		case DATATYPE_INT8:
			data_def.bufsize = data_def.nchans * data_def.nsamples * WORDSIZE_INT8;
			break;
		case DATATYPE_INT16:
			data_def.bufsize = data_def.nchans * data_def.nsamples * WORDSIZE_INT16;
			break;
		case DATATYPE_INT32:
			data_def.bufsize = data_def.nchans * data_def.nsamples * WORDSIZE_INT32;
			break;
		case DATATYPE_INT64:
			data_def.bufsize = data_def.nchans * data_def.nsamples * WORDSIZE_INT64;
			break;
		case DATATYPE_FLOAT32:
			data_def.bufsize = data_def.nchans * data_def.nsamples * WORDSIZE_FLOAT32;
			break;
		case DATATYPE_FLOAT64:
			data_def.bufsize = data_def.nchans * data_def.nsamples * WORDSIZE_FLOAT64;
			break;
		default:
			mexErrMsgTxt("Unrecognised data type");
	}
    data.buf = (void *)mxGetPr(field); /* directly point to the memory that is managed by Matlab */
  
  /* construct a PUT_DAT request */
	request_def.bufsize = append(&request.buf, request_def.bufsize, data.def, sizeof(datadef_t));
	request_def.bufsize = append(&request.buf, request_def.bufsize, data.buf, data_def.bufsize);
  
  /* write the request, read the response */
	result = clientrequest(server, &request, &response);
  
  /* the request structure is not needed any more, we free ->buf, the rest is local */
	if (request.buf) {
		FREE(request.buf);
	}
	
	if (result == 0) {
		/* check that the response is PUT_OK */
		if (!response)
			mexErrMsgTxt("unknown error in response\n");
		else if (!response->def) {
			FREE(response->buf);
			FREE(response);
			mexErrMsgTxt("unknown error in response\n");
		}
		else if (response->def->command!=PUT_OK) {
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
