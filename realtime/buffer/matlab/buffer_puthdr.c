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

int buffer_puthdr(int server, mxArray * plhs[], const mxArray * prhs[])
{
	int fieldnumber;
	mxArray *field;
	int result;
  
	message_t     request;
	messagedef_t  request_def;
	message_t    *response = NULL;
	header_t      header;
	headerdef_t   header_def;
  
  /* allocate the header */
	header.def = &header_def;
	header.buf = NULL;
	header_def.bufsize   = 0;
  
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
	request_def.bufsize = append(&request.buf, request_def.bufsize, header.def, sizeof(headerdef_t));
	/* request_def.bufsize = append(&request.buf, request_def.bufsize, header.buf, header_def.bufsize);    -- this is empty */
  
  /* the header structure is not needed any more, but everything's local */
  
  /* write the request, read the response */
	result = clientrequest(server, &request, &response);
  
  /* the request structure is not needed any more, but only .buf needs to be free'd */
    FREE(request.buf);
	
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
