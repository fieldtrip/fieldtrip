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

int buffer_putprp(int server, mxArray * plhs[], const mxArray * prhs[])
{
  int fieldnumber, bufsize;
  mxArray *field;
  int result;
  
  message_t   *response = NULL;
  
  /* Place the fixed-length property and request structures directly on the stack 
     This avoids a lot of clean-up problems and allows us to leave early using mexErrMsgTxt
  */
  propertydef_t property_def;
  property_t    property;
  message_t     request;
  messagedef_t  request_def;
  
  property.def = &property_def;
  property.buf = NULL;
  property_def.bufsize = 0;

  request.def  = &request_def;
  request.buf = NULL;
  request_def.version = VERSION;
  request_def.command = PUT_PRP;
  request_def.bufsize = 0;  
  
  /* define the property, it has the fields type_type type_numel value_type value_numel */
  
  /* FIXME loop over mutiple properties in case of a property-array */
  
  fieldnumber = mxGetFieldNumber(prhs[0], "type_type");
  if (fieldnumber<0) 
      mexErrMsgTxt("field is missing");

  field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
  if (!mxIsNumeric(field) || mxIsEmpty(field))
      mexErrMsgTxt("invalid data type");
  property_def.type_type = (UINT32_T)mxGetScalar(field) ;

  
  fieldnumber = mxGetFieldNumber(prhs[0], "type_numel");
  if (fieldnumber<0) 
      mexErrMsgTxt("field is missing");
  field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
  if (!mxIsNumeric(field) || mxIsEmpty(field))
      mexErrMsgTxt("invalid data type");
  property_def.type_numel = (UINT32_T)mxGetScalar(field) ;

  
  fieldnumber = mxGetFieldNumber(prhs[0], "value_type");
  if (fieldnumber<0) 
      mexErrMsgTxt("field is missing");
  field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
  if (!mxIsNumeric(field) || mxIsEmpty(field))
      mexErrMsgTxt("invalid data type");
  property_def.value_type = (UINT32_T)mxGetScalar(field) ;
  
  
  fieldnumber = mxGetFieldNumber(prhs[0], "value_numel");
  if (fieldnumber<0) 
      mexErrMsgTxt("field is missing");
  field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
  if (!mxIsNumeric(field) || mxIsEmpty(field)) 
      mexErrMsgTxt("invalid data type");
  property_def.value_numel = (UINT32_T)mxGetScalar(field) ;
  
  
  fieldnumber = mxGetFieldNumber(prhs[0], "bufsize");
  if (fieldnumber<0) 
      mexErrMsgTxt("field is missing");
  field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
  if (!mxIsNumeric(field) || mxIsEmpty(field)) 
      mexErrMsgTxt("invalid data type");
  bufsize = (UINT32_T)mxGetScalar(field) ;
  
  fieldnumber = mxGetFieldNumber(prhs[0], "buf");
  if (fieldnumber<0) 
      mexErrMsgTxt("field is missing");
  field = mxGetFieldByNumber(prhs[0], 0, fieldnumber);
  if (!mxIsUint8(field) || mxIsComplex(field)) 
      mexErrMsgTxt("buffer should be real-valued uint8");
  if (mxGetNumberOfElements(field)!=bufsize) 
      mexErrMsgTxt("length of buffer does not correspond to bufsize");
  property_def.bufsize = append(&property.buf, property_def.bufsize, mxGetData(field), bufsize);
  
  /* construct a PUT_PRP request */
  request_def.bufsize = append(&request.buf, request_def.bufsize, property.def, sizeof(propertydef_t));
  request_def.bufsize = append(&request.buf, request_def.bufsize, property.buf, property_def.bufsize);
  
  /* write the request, read the response */
  result = clientrequest(server, &request, &response);
  
  /* the request structure is not needed any more */
  FREE(request.buf);
  
	if (result == 0) {
		/* check that the response is PUT_OK */
		if (!response) {
			/* nothing to clean up */
			mexErrMsgTxt("unknown error in response\n");
		}
		if (!response->def) {
			FREE(response->buf);
			FREE(response);
			mexErrMsgTxt("unknown error in response\n");
		}
		if (response->def->command!=PUT_OK) {
			result = response->def->command;
		}
	}
	FREE(response->def);
	FREE(response->buf);
	FREE(response);
	return result;
}

