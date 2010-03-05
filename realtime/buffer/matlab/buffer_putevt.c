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

int buffer_putevt(int server, mxArray * plhs[], const mxArray * prhs[])
{
	int fieldnumber;
	mxArray *field;
	int result;
  
	message_t    request;
	messagedef_t request_def;
	message_t   *response = NULL;
	event_t      event;
	eventdef_t   event_def;
  
  /* allocate the event */
	event.def = &event_def;
	event.buf = NULL;
	event_def.bufsize = 0;
  
  /* allocate the request message */
	request.buf = NULL;
	request.def = &request_def;
	request_def.version = VERSION;
	request_def.command = PUT_EVT;
	request_def.bufsize = 0;
  
  /* define the event, it has the fields type_type type_numel value_type value_numel sample offset duration */
  
	if (mxGetNumberOfElements(prhs[0])!=1)
		mexErrMsgTxt("Only one event can be put into the buffer at a time.");
  
	fieldnumber = mxGetFieldNumber(prhs[0], "type_type");
	if (fieldnumber<0) 
		mexErrMsgTxt("field 'type_type' is missing");
	field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
	if (!mxIsNumeric(field) || mxIsEmpty(field)) 
		mexErrMsgTxt("invalid data type for 'type_type'");
	event_def.type_type = (UINT32_T)mxGetScalar(field) ;

	fieldnumber = mxGetFieldNumber(prhs[0], "type_numel");
	if (fieldnumber<0) 
		mexErrMsgTxt("field 'type_numel' is missing");
	field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
	if (!mxIsNumeric(field) || mxIsEmpty(field)) 
		mexErrMsgTxt("invalid data type for 'type_numel'");
	event_def.type_numel = (UINT32_T)mxGetScalar(field) ;
  
	fieldnumber = mxGetFieldNumber(prhs[0], "value_type");
	if (fieldnumber<0) 
		mexErrMsgTxt("field 'value_type' is missing");
		
    field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
	if (!mxIsNumeric(field) || mxIsEmpty(field)) 
		mexErrMsgTxt("invalid data type for 'value_type'");
	event_def.value_type = (UINT32_T)mxGetScalar(field) ;
  
	fieldnumber = mxGetFieldNumber(prhs[0], "value_numel");
	if (fieldnumber<0) 
		mexErrMsgTxt("field 'value_numel' is missing");
	field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
	if (!mxIsNumeric(field) || mxIsEmpty(field)) 
		mexErrMsgTxt("invalid data type for 'value_numel'");
	event_def.value_numel = (UINT32_T)mxGetScalar(field) ;
  
	fieldnumber = mxGetFieldNumber(prhs[0], "sample");
	if (fieldnumber<0) 
		mexErrMsgTxt("field 'sample' is missing");
	field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
    if (!mxIsNumeric(field) || mxIsEmpty(field)) 
		mexErrMsgTxt("invalid data type for 'sample'");
	event_def.sample = (INT32_T)mxGetScalar(field) ;

	fieldnumber = mxGetFieldNumber(prhs[0], "offset");
	if (fieldnumber<0) 
		mexErrMsgTxt("field 'offset' is missing");
    field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
	if (!mxIsNumeric(field) || mxIsEmpty(field)) 
      mexErrMsgTxt("invalid data type for 'offset'");
	event_def.offset = (INT32_T)mxGetScalar(field);
  
	fieldnumber = mxGetFieldNumber(prhs[0], "duration");
	if (fieldnumber<0) 
		mexErrMsgTxt("field 'duration' is missing");
	field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
	if (!mxIsNumeric(field) || mxIsEmpty(field)) 
		mexErrMsgTxt("invalid data type for 'duration'");
	event_def.duration = (INT32_T)mxGetScalar(field) ;
  
	fieldnumber = mxGetFieldNumber(prhs[0], "bufsize");
	if (fieldnumber<0) 
		mexErrMsgTxt("field 'bufsize' is missing");
	field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
	if (!mxIsNumeric(field) || mxIsEmpty(field)) 
		mexErrMsgTxt("invalid data type for 'bufsize'");
	event_def.bufsize = (INT32_T)mxGetScalar(field) ;
  
	fieldnumber = mxGetFieldNumber(prhs[0], "buf");
	if (fieldnumber<0) 
		mexErrMsgTxt("field 'buf' is missing");
    field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
    if (!mxIsUint8(field)) 
      mexErrMsgTxt("invalid data type for 'buf'");
	if (mxGetNumberOfElements(field) != event_def.bufsize ) 
		mexErrMsgTxt("invalid number of elements (buf)");
	event.buf = mxGetData(field);
  
  /* construct a PUT_EVT request */
	request_def.bufsize = append(&request.buf, request_def.bufsize, event.def, sizeof(eventdef_t));
	request_def.bufsize = append(&request.buf, request_def.bufsize, event.buf, event_def.bufsize);
  
  /* the event structure is not needed any more, but everything's local, and .buf is from a Matlab array */
  
  /* write the request, read the response */
	result = clientrequest(server, &request, &response);
  
  /* the request structure is not needed any more, everything apart from request.buf is local */
    FREE(request.buf);
	
	if (result == 0) {	/* no communication errors */
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
