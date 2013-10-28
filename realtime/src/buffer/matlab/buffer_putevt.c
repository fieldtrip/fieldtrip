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
#include "buffer_mxutils.h"

static int fieldnum_type, fieldnum_value, fieldnum_sample, fieldnum_offset, fieldnum_duration;


unsigned int add_event_from_matlab(unsigned int bufsize, void **buf, const mxArray *E, int index) {
	const mxArray *field;
	eventdef_t evdef;
	unsigned int type_bytes, value_bytes;
	unsigned int newsize;
	char *type_buf = NULL, *value_buf = NULL;
	
	field = mxGetFieldByNumber(E, index, fieldnum_type);
	evdef.type_type = ft_type_from_array(field);
	switch(evdef.type_type) {
		case DATATYPE_UNKNOWN:
			mexErrMsgTxt("invalid data type for 'type'");
			break;
		case DATATYPE_CHAR:
			type_buf = mxArrayToString(field);
			if (type_buf == NULL) mexErrMsgTxt("Unexpected error when converting 'type' string");
			type_bytes = evdef.type_numel = strlen(type_buf);
			break;
		default:
			evdef.type_numel = mxGetNumberOfElements(field);
			type_bytes = wordsize_from_type(evdef.type_type) * evdef.type_numel;
			type_buf = (char *) mxGetData(field);
	}
	
	field = mxGetFieldByNumber(E, index, fieldnum_value);
	evdef.value_type = ft_type_from_array(field);
	switch(evdef.value_type) {
		case DATATYPE_UNKNOWN:
			mexErrMsgTxt("invalid data type for 'value'");
			break;
		case DATATYPE_CHAR:
			value_buf = mxArrayToString(field);
			if (value_buf == NULL) mexErrMsgTxt("Unexpected error when converting 'value' string");
			value_bytes = evdef.value_numel = strlen(value_buf);
			break;
		default:
			evdef.value_numel = mxGetNumberOfElements(field);
			value_bytes = wordsize_from_type(evdef.value_type) * evdef.value_numel;
			value_buf = (char *) mxGetData(field);
	}

	evdef.bufsize = type_bytes + value_bytes;
	
	if (fieldnum_sample < 0) {
		/* "sample" field not given? -> will be inserted by server */
		evdef.sample = EVENT_AUTO_SAMPLE;  
	} else {
		field = mxGetFieldByNumber(E, index, fieldnum_sample);
		if (mxIsEmpty(field)) {
			/* also empty "sample" fields will be auto-translated */
			evdef.sample = EVENT_AUTO_SAMPLE;  
		} else {
			if (!mxIsNumeric(field)) {
				mexErrMsgTxt("invalid data type for 'sample'");
			}
			evdef.sample = (UINT32_T) mxGetScalar(field) - 1; /* 0-based index on protocol level */
		}
	}
	
	field = mxGetFieldByNumber(E, index, fieldnum_offset);
	if (!mxIsNumeric(field) || mxIsEmpty(field)) 
		mexErrMsgTxt("invalid data type for 'offset'");
	evdef.offset = (INT32_T) mxGetScalar(field);
	
	field = mxGetFieldByNumber(E, index, fieldnum_duration);
	if (!mxIsNumeric(field) || mxIsEmpty(field)) 
		mexErrMsgTxt("invalid data type for 'duration'");
	evdef.duration = (UINT32_T) mxGetScalar(field);
	
	newsize = ft_mx_append(buf, bufsize, &evdef, sizeof(evdef));
	newsize = ft_mx_append(buf, newsize, type_buf, type_bytes);
	newsize = ft_mx_append(buf, newsize, value_buf, value_bytes);
	
	/* For strings, we need to free the corresponding buffers,
		otherwise they are just pointers to existing Matlab arrays
	*/
	if (evdef.type_type == DATATYPE_CHAR) mxFree(type_buf);
	if (evdef.value_type == DATATYPE_CHAR) mxFree(value_buf);
	
	return newsize;
}

int buffer_putevt(int server, mxArray * plhs[], const mxArray * prhs[])
{
	mxArray *field;
	int i, Nev, result;
  
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
  
	Nev = mxGetNumberOfElements(prhs[0]);
	if (Nev == 0) {
    result = 0; /* negative numbers indicate an error of clientrequest, 0 indicates that it worked fine */
    return result;
  }
	
	if (!mxIsStruct(prhs[0])) mexErrMsgTxt("events must be given as structs.");
		
	fieldnum_type = mxGetFieldNumber(prhs[0], "type");
	if (fieldnum_type < 0) mexErrMsgTxt("field 'type' is missing");
	
	fieldnum_value = mxGetFieldNumber(prhs[0], "value");
	if (fieldnum_value < 0) mexErrMsgTxt("field 'value' is missing");
	
	fieldnum_sample = mxGetFieldNumber(prhs[0], "sample");
	/*
		If no "sample" field is given, we use EVENT_AUTO_SAMPLE and have 
		the server insert the right number automatically.
		
		Old code:
		if (fieldnum_sample < 0) mexErrMsgTxt("field 'sample' is missing");
	*/
	 
	fieldnum_offset = mxGetFieldNumber(prhs[0], "offset");
	if (fieldnum_offset < 0) mexErrMsgTxt("field 'offset' is missing");
	
	fieldnum_duration = mxGetFieldNumber(prhs[0], "duration");
	if (fieldnum_duration < 0) mexErrMsgTxt("field 'duration' is missing");
	
	/* add one event at a time to the request */
	for (i=0;i<Nev;i++) {
		request_def.bufsize = add_event_from_matlab(request_def.bufsize, &request.buf, prhs[0], i);
	}
  
  /* write the request, read the response */
	result = clientrequest(server, &request, &response);
  
  /* the request structure is not needed any more, everything apart from request.buf is local */
	mxFree(request.buf);
	
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
