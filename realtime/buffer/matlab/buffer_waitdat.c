/*
 * Copyright (C) 2010, Stefan Klanke
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 */

#include "mex.h"
#include "matrix.h"
#include "buffer.h"

int buffer_waitdat(int server, mxArray * plhs[], const mxArray * prhs[])
{

	int fieldnumber;
	mxArray *field;
	int result;
	double timeout;
	const double *pr;
  
	message_t     request;
	messagedef_t  request_def;
	message_t    *response = NULL;
	waitdef_t 	  waitdef;
  
	request.def = &request_def;
	request.buf = &waitdef;
	request_def.version = VERSION;
	request_def.command = WAIT_DAT;
	request_def.bufsize = sizeof(waitdef_t);
	
	
	if (mxGetNumberOfElements(prhs[0])!=3 || !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) 
		mexErrMsgTxt("3rd parameter must be contain 3 double values [nsamples, nevents, timeout].");
	
	pr = mxGetPr(prhs[0]);
	waitdef.threshold.nsamples = (pr[0]<0) ? 0xFFFFFFFF : (UINT32_T) pr[0];
	waitdef.threshold.nevents  = (pr[1]<0) ? 0xFFFFFFFF : (UINT32_T) pr[1];
	waitdef.milliseconds = (pr[2]<0) ? 0 : (UINT32_T) pr[2];
	
	/* write the request, read the response */
	result = clientrequest(server, &request, &response);
  	
	if (result == 0) {
		/* check that the response is WAIT_OK */
		if (!response)
			mexErrMsgTxt("unknown error in response\n");
		if (!response->def) {
			FREE(response->buf);
			FREE(response);
			mexErrMsgTxt("unknown error in response\n");
		}
		if (response->def->command!=WAIT_OK || response->buf == NULL || response->def->bufsize != sizeof(samples_events_t)) {
			result = response->def->command;
		} else {
			const char *field_names[2] = {"nsamples","nevents"};
			samples_events_t *nes = (samples_events_t *) response->buf;
			
			plhs[0] = mxCreateStructMatrix(1, 1, 2, field_names);
			mxSetFieldByNumber(plhs[0], 0, 0, mxCreateDoubleScalar((double)(nes->nsamples)));
			mxSetFieldByNumber(plhs[0], 0, 1, mxCreateDoubleScalar((double)(nes->nevents)));
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
