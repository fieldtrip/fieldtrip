/*
 * Copyright (C) 2008, Robert Oostenveld
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 */

#include <stdio.h>
#include <stdlib.h>

#include "buffer.h"

/*******************************************************************************
 * this function is called by the client
 * it takes care that the request is processed by the buffer
 *******************************************************************************/
int clientrequest(int server, const message_t *request, message_t **response_ptr) {
    int verbose = 0;

	if (verbose>0) fprintf(stderr, "clientrequest: server = %d\n", server);
	if (verbose>0) print_request(request->def);

	if (server<0) {
		fprintf(stderr, "clientrequest: invalid value for server (%d)\n", server);
		return -1;
	}

	else if (server==0) {
		/* use direct memory acces to the buffer */
		if (dmarequest(request, response_ptr)!=0)
			return -2;
	}

	else if (server>0) {
		/* use TCP connection to the buffer */
		if (tcprequest(server, request, response_ptr)!=0)
			return -3;
	}

	if (verbose>0) print_response((*response_ptr)->def);

	/* everything went fine */
	return 0;
}
