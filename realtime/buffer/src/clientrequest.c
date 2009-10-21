/*
 * Copyright (C) 2008, Robert Oostenveld
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 * $Log: clientrequest.c,v $
 * Revision 1.4  2009/01/23 19:47:50  roboos
 * changed verbosity
 *
 * Revision 1.3  2008/05/22 09:49:32  roboos
 * some small changes by Jurgen Mellinger related to compatibility with Borland
 * moved linux and socket specific includes and defines into seperate header files
 *
 * Revision 1.2  2008/04/24 15:39:44  roboos
 * more variability in return values, changed verbosity
 *
 * Revision 1.1  2008/03/08 10:30:20  roboos
 * renamed write_request to clientrequest, moved tcp specific code to seperate tcprequest function
 *
 * Revision 1.1  2008/03/07 14:49:39  roboos
 * new implementation, the function itself decides whether to use sockets or direct memcpy access
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include "unix_includes.h"
#include <pthread.h>

#include "buffer.h"
#include "message.h"

/*******************************************************************************
 * this function is called by the client
 * it takes care that the request is processed by the buffer
 *******************************************************************************/
int clientrequest(int server, message_t *request, message_t **response_ptr) {
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

