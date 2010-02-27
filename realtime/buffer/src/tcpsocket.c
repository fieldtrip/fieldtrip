/*
 * Copyright (C) 2008, Robert Oostenveld
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 *
 */

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "buffer.h"

#ifdef ENABLE_POLLING
  #include <poll.h>
#endif

#define THREADSLEEP 1000000  // in microseconds
#define POLLSLEEP   100     // in microseconds

extern pthread_mutex_t mutexsocketcount;
extern int socketcount;
extern pthread_mutex_t mutexthreadcount;
extern int threadcount;

/* this function deals with the incoming client request */
void *tcpsocket(void *arg) {
	int n;
	int status = 0, verbose = 0;
	int oldcancelstate, oldcanceltype;

#ifdef ENABLE_POLLING
	struct pollfd fds;
#endif

	// these are used for communication over the TCP socket
	int client = 0;
	message_t *request = NULL, *response = NULL;

	/* the connection to the client has been made by the server */
	client = (int)arg;

    /* this is for debugging */
    pthread_mutex_lock(&mutexsocketcount);
    socketcount++;
    pthread_mutex_unlock(&mutexsocketcount);

	/* this is to prevent closing the thread at an unwanted moment and memory from leaking */
	pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &oldcancelstate);
	pthread_setcanceltype(PTHREAD_CANCEL_DEFERRED, &oldcanceltype);
	pthread_cleanup_push(cleanup_message, &request);
	pthread_cleanup_push(cleanup_message, &response)
	pthread_cleanup_push(cleanup_socket, &client);

    if (verbose>1) fprintf(stderr, "tcpsocket: client = %d, socketcount = %d, threadcount = %d\n", client, socketcount, threadcount);

	/* keep processing messages untill the connection is closed */
	while (1) {

	request       = (message_t*)malloc(sizeof(message_t));
	request->def  = (messagedef_t*)malloc(sizeof(messagedef_t));
	request->buf  = NULL;

#ifdef ENABLE_POLLING
		/* wait for data to become available or until the connection is closed */
		/* thohar: i think this is not neccessary as we dont need a timeout. */
		/* roboos: we need it to detect when the socket is closed by the client */
		while (1) {
			fds.fd      = client;
			fds.events  = POLLIN | POLLRDNORM | POLLRDBAND | POLLPRI | POLLOUT | POLLWRNORM | POLLWRBAND | POLLERR | POLLNVAL;
			fds.revents = 0;

			if (poll(&fds, 1, 1)==-1) {
				perror("poll");
				goto cleanup;
			}

			if (fds.revents & POLLHUP)
				goto cleanup;				// the connection has been closed
			else if (fds.revents & POLLERR)
				goto cleanup;				// the connection has been closed
			else if (fds.revents & POLLIN)
				break;						// data is available, process the message
			else
				usleep(POLLSLEEP);			// wait for data or closed connection
		}
#endif

		if ((n = bufread(client, request->def, sizeof(messagedef_t))) != sizeof(messagedef_t)) {
			if (verbose>0) fprintf(stderr, "tcpsocket: packet size = %d, should be %d\n", n, sizeof(messagedef_t));
			goto cleanup;
		}

		if (request->def->version!=VERSION) {
			if (verbose>0) fprintf(stderr, "tcpsocket: incorrect request version\n");
			goto cleanup;
		}

		if (request->def->bufsize>0) {
			request->buf = malloc(request->def->bufsize);
			if ((n = bufread(client, request->buf, request->def->bufsize)) != request->def->bufsize) {
				if (verbose>0) fprintf(stderr, "tcpsocket: read size = %d, should be %d\n", n, request->def->bufsize);
				goto cleanup;
			}
		}

		if (verbose>1) print_request(request->def);
		if (verbose>1) print_buf(request->buf, request->def->bufsize);

		if ((status = dmarequest(request, &response)) != 0) {
			if (verbose>0) fprintf(stderr, "tcpsocket: an unexpected error occurred\n");
			goto cleanup;
		}

		if (verbose>1) print_response(response->def);
		if (verbose>1) print_buf(request->buf, request->def->bufsize);

		if ((n = bufwrite(client, response->def, sizeof(messagedef_t)))!=sizeof(messagedef_t)) {
			if (verbose>0) fprintf(stderr, "tcpsocket: write size = %d, should be %d\n", n, sizeof(messagedef_t));
			goto cleanup;
		}
		if ((n = bufwrite(client, response->buf, response->def->bufsize))!=response->def->bufsize) {
			if (verbose>0) fprintf(stderr, "tcpsocket: write size = %d, should be %d\n", n, response->def->bufsize);
			goto cleanup;
		}

		cleanup_message(&request);
		cleanup_message(&response);
        request = NULL;
        response = NULL;

	} /* while (1) */

cleanup:
	/* from now on it is safe to cancel the thread */
	pthread_setcancelstate(oldcancelstate, NULL);
	pthread_setcanceltype(oldcanceltype, NULL);

	pthread_cleanup_pop(1); // request
	pthread_cleanup_pop(1); // response
	pthread_cleanup_pop(1); // socket

    /* this is for debugging */
    pthread_mutex_lock(&mutexsocketcount);
    socketcount--;
    pthread_mutex_unlock(&mutexsocketcount);

    /* this is for debugging */
    pthread_mutex_lock(&mutexthreadcount);
    threadcount--;
    pthread_mutex_unlock(&mutexthreadcount);

	pthread_exit(NULL);
	return NULL;
}

