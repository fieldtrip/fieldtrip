/*
 * Copyright (C) 2008, Robert Oostenveld
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 * $Log: buffer_flushdat.c,v $
 * Revision 1.3  2008/10/29 20:46:12  roboos
 * consistent use of open_connection and close_connection
 * there were some incorrect uses of close() that on windows did not actually close, resulting in the buffer running out of sockets/threads after prolonged use
 *
 * Revision 1.2  2008/07/09 13:34:21  roboos
 * small change in verbose output, using verbose=0|1
 *
 * Revision 1.1  2008/06/19 20:48:32  roboos
 * added support for flushing header, data and events
 *
 *
 */

#include "mex.h"
#include "matrix.h"
#include "buffer.h"

int buffer_flushdat(int server, mxArray *plhs[], const mxArray *prhs[])
{
	int verbose = 0;
	int result = 0;
	message_t *request  = NULL;
	message_t *response = NULL;

	/* allocate the elements that will be used in the communication */
	request      = malloc(sizeof(message_t));
	request->def = malloc(sizeof(messagedef_t));
	request->buf = NULL;
	request->def->version = VERSION;
	request->def->command = FLUSH_DAT;
	request->def->bufsize = 0;

	if (verbose) print_request(request->def);
	result = clientrequest(server, request, &response);
	if (verbose) print_response(response->def);

	if (result==0) {
		if (response->def->command!=FLUSH_OK) {
			result = response->def->command;
		}
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

	return result;
}

