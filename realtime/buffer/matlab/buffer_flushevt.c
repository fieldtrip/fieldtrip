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

int buffer_flushevt(int server, mxArray *plhs[], const mxArray *prhs[])
{
	int verbose = 0;
	int result;
  
	message_t *request  = NULL;
	message_t *response = NULL;

	/* allocate the elements that will be used in the communication */
	request      = malloc(sizeof(message_t));
	request->def = malloc(sizeof(messagedef_t));
	request->buf = NULL;
	request->def->version = VERSION;
	request->def->command = FLUSH_EVT;
	request->def->bufsize = 0;

	if (verbose) print_request(request->def);
	result = clientrequest(server, request, &response);
	if (verbose) print_response(response->def);

	if (result == 0) {
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

