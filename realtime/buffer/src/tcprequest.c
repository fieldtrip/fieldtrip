/*
 * Copyright (C) 2008, Robert Oostenveld
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 * $Log: tcprequest.c,v $
 * Revision 1.3  2008/05/22 09:49:32  roboos
 * some small changes by Jurgen Mellinger related to compatibility with Borland
 * moved linux and socket specific includes and defines into seperate header files
 *
 * Revision 1.2  2008/03/13 13:37:55  roboos
 * remove include pthread
 *
 * Revision 1.1  2008/03/08 10:38:56  roboos
 * moved code to standalone function, seperate from other ways of handling the request (e.g. dma and udp)
 * no functional changes
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include "unix_includes.h"

#include "buffer.h"
#include "message.h"

/*******************************************************************************
 * commnuicate with the buffer through TCP
 *******************************************************************************/
int tcprequest(int server, message_t *request, message_t **response_ptr) {
	int n;

	// this will hold the response
	message_t *response;
	response      = (message_t*)malloc(sizeof(message_t));
	response->def = (messagedef_t*)malloc(sizeof(messagedef_t));
	response->buf = NULL;
	// the response should be passed to the calling function, where it should be freed
	*response_ptr = response;

	/* send the request to the server, first the message definition */
	if ((n = bufwrite(server, request->def, sizeof(messagedef_t)))!=sizeof(messagedef_t)) {
		fprintf(stderr, "write size = %d, should be %d\n", n, sizeof(messagedef_t));
		goto cleanup;
	}

	/* send the request to the server, then the message payload */
	if ((n = bufwrite(server, request->buf, request->def->bufsize))!=request->def->bufsize) {
		fprintf(stderr, "write size = %d, should be %d\n", n, request->def->bufsize);
		goto cleanup;
	}

	/* read the response from the server, first the message definition */
	if ((n = bufread(server, response->def, sizeof(messagedef_t))) != sizeof(messagedef_t)) {
		fprintf(stderr, "packet size = %d, should be %d\n", n, sizeof(messagedef_t));
		goto cleanup;
	}

	if (response->def->version!=VERSION) {
		fprintf(stderr, "incorrect version\n");
		goto cleanup;
	}

	/* read the response from the server, then the message payload */
	if (response->def->bufsize>0) {
		response->buf = malloc(response->def->bufsize);
		if ((n = bufread(server, response->buf, response->def->bufsize)) != response->def->bufsize) {
			fprintf(stderr, "read size = %d, should be %d\n", n, response->def->bufsize);
			goto cleanup;
		}
	}

	/* everything went fine, return with the response */
	//print_response(response->def);
	return 0;

cleanup:
	/* there was a problem, clear the response and return */
	FREE(response->def);
	FREE(response->buf);
	FREE(response);
	response_ptr = NULL;
	return -1;
}



