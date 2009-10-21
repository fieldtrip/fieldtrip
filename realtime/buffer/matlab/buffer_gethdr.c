/*
 * Copyright (C) 2008, Robert Oostenveld
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 * $Log: buffer_gethdr.c,v $
 * Revision 1.9  2009/06/30 11:35:48  roboos
 * pass more explicit error from buffer onto matlab
 *
 * Revision 1.8  2008/10/29 20:46:12  roboos
 * consistent use of open_connection and close_connection
 * there were some incorrect uses of close() that on windows did not actually close, resulting in the buffer running out of sockets/threads after prolonged use
 *
 * Revision 1.7  2008/07/09 13:34:21  roboos
 * small change in verbose output, using verbose=0|1
 *
 * Revision 1.6  2008/05/22 09:27:22  roboos
 * fixed some issues with Borland compiler, correct pointer arithmetic, declarations at beginning
 *
 * Revision 1.5  2008/03/09 22:36:38  roboos
 * implemented new buffer mex-function, which acts as dispatcher for all functionality
 * changed buffer_xxx code from mex-function into c-function to be included in dispatcher
 *
 * Revision 1.4  2008/02/26 21:43:25  roboos
 * renamed packet_t structure definition into messagedef_t, added message_t structure (contains def+buf)
 *
 * Revision 1.3  2008/02/20 13:49:17  roboos
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

#define NUMBER_OF_FIELDS 6

void buffer_gethdr(char *hostname, int port, mxArray *plhs[], const mxArray *prhs[])
{
	int server;
	int verbose = 0;
	size_t n;
  char msg[512];

	message_t *request  = NULL;
	message_t *response = NULL;
	header_t  *header   = NULL;

	/* this is for the Matlab specific output */
	const char *field_names[NUMBER_OF_FIELDS] = {"nchans", "nsamples", "nevents", "fsample", "data_type", "bufsize"};

	/* allocate the elements that will be used in the communication */
	request      = malloc(sizeof(message_t));
	request->def = malloc(sizeof(messagedef_t));
	request->buf = NULL;
	request->def->version = VERSION;
	request->def->command = GET_HDR;
	request->def->bufsize = 0;

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
		header      = malloc(sizeof(header_t));
		header->def = response->buf;
		header->buf = (char *)response->buf + sizeof(headerdef_t);
		if (verbose) print_headerdef(header->def);

		plhs[0] = mxCreateStructMatrix(1, 1, NUMBER_OF_FIELDS, field_names);
		mxSetFieldByNumber(plhs[0], 0, 0, mxCreateDoubleScalar((double)(header->def->nchans)));
		mxSetFieldByNumber(plhs[0], 0, 1, mxCreateDoubleScalar((double)(header->def->nsamples)));
		mxSetFieldByNumber(plhs[0], 0, 2, mxCreateDoubleScalar((double)(header->def->nevents)));
		mxSetFieldByNumber(plhs[0], 0, 3, mxCreateDoubleScalar((double)(header->def->fsample)));
		mxSetFieldByNumber(plhs[0], 0, 4, mxCreateDoubleScalar((double)(header->def->data_type)));
		mxSetFieldByNumber(plhs[0], 0, 5, mxCreateDoubleScalar((double)(header->def->bufsize)));
		FREE(header);
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

