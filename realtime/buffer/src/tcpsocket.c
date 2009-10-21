/*
 * Copyright (C) 2008, Robert Oostenveld
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 * $Log: tcpsocket.c,v $
 * Revision 1.16  2009/01/23 19:47:50  roboos
 * changed verbosity
 *
 * Revision 1.15  2009/01/23 17:55:34  roboos
 * set request and response to NULL after cleanup, to ensure that they are not free'd twice
 *
 * Revision 1.14  2009/01/23 08:26:44  roboos
 * fixed a serious bug that caused a lot of memory to leak (in fact all packets that were sent over the socket would eventually leack away), both on the client and server side
 *
 * Revision 1.13  2009/01/21 20:54:35  roboos
 * added some debugging code to keep track of number of sockets and threads (both seem ok)
 * cleaned up the debugging: more consistent use of verbose flag and function name in fprintf
 *
 * Revision 1.12  2009/01/21 08:47:22  roboos
 * fixed bug in thread counting
 *
 * Revision 1.11  2009/01/21 08:33:01  roboos
 * added some code for debugging the leaking of sockets and threads
 *
 * Revision 1.10  2008/07/09 10:38:56  roboos
 * changed verbosity to 0, some changes in whitespace
 *
 * Revision 1.9  2008/05/22 13:39:01  roboos
 * moved declarations to top
 *
 * Revision 1.8  2008/05/22 09:49:32  roboos
 * some small changes by Jurgen Mellinger related to compatibility with Borland
 * moved linux and socket specific includes and defines into seperate header files
 *
 * Revision 1.7  2008/04/24 15:52:20  roboos
 * disabled the section on polling using a define, changed verbosity
 *
 * Revision 1.6  2008/03/26 14:36:52  thohar
 * minor portability stuff
 *
 * Revision 1.5  2008/03/23 16:53:21  roboos
 * implemented thread cancelation, reactivated polling code, needed to detect closed socked
 *
 * Revision 1.4  2008/03/13 13:40:54  roboos
 * changed handling input argument, added thread cancelation, changed cleanup (not yet final, since that would mean that it moves to the cleanup function)
 *
 * Revision 1.3  2008/03/13 12:35:23  thohar
 * added headers for win32 build
 * removed polling as it is not neccessary and just eats up cpu-time
 *
 * Revision 1.2  2008/03/08 10:33:48  roboos
 * changed some code to reflect the renamed low-level functions, no functional changes
 *
 * Revision 1.1  2008/03/08 09:41:02  roboos
 * renamed buffer_thread function to tcpserver, idem for the c-file
 * renamed socket_thread function to tcpsocket, idem for the c-file
 *
 * Revision 1.12  2008/03/02 13:29:23  roboos
 * seperated socket and buffer-handing code
 * adde polling to the socket section
 * allow connections to remain open (i.e. multiple request/response pairs)
 *
 * Revision 1.11  2008/02/27 10:17:30  roboos
 * Use append() and message_t and messagedef_t in constructing the response in memory, delay the writing of the response to after the switch-ladder. All responses are now written the same way, equal to reading them. This allows for the next step, which should be to seperate this code into a buffer_request section and a socket/network specific section.
 *
 * Revision 1.10  2008/02/26 21:43:26  roboos
 * renamed packet_t structure definition into messagedef_t, added message_t structure (contains def+buf)
 *
 * Revision 1.9  2008/02/20 13:41:06  roboos
 * malloc the datasel and eventsel, also if specified in buf (free buf immediately)
 * implemented delayed writing on socket for GET_DAT
 *
 * Revision 1.8  2008/02/20 07:14:02  roboos
 * added support for the other data types
 *
 * Revision 1.7  2008/02/19 17:22:38  roboos
 * added support for int32 data buffer
 * changed the verbosity
 *
 * Revision 1.6  2008/02/19 10:41:28  roboos
 * re-enabled the test for reading outside of the buffer
 *
 * Revision 1.5  2008/02/19 10:22:56  roboos
 * added consistent copyrights and log message to each of the files
 *
 * Revision 1.4  2008/02/19 09:29:25  roboos
 * replaced all read/write calls with bufread/bufwrite calls
 * fixed problem in offset calculation (accidentaly removed the (int) cast, which was meant for rounding off)
 *
 * Revision 1.3  2008/02/19 09:08:32  roboos
 * moved some mutexes around
 * fixed a nasty bug in GET_EVT (offset beyond end of event buffer)
 * many small changes
 *
 * Revision 1.2  2008/02/18 12:13:46  roboos
 * moved executable from buffer to demo
 * fixed bugs in sinewave and socket for events
 * stripped down the eventdef_t fields
 * many small changes
 *
 * Revision 1.1  2008/02/18 10:05:27  roboos
 * restructured the directory layout, copied all code into src, added directory for external code
 *
 * Revision 1.7  2008/02/14 16:52:05  roboos
 * implemented GET_EVT in analogy with GET_DAT
 *
 * Revision 1.6  2008/02/13 13:59:03  roboos
 * made a start with implementing GET_EVT
 *
 * Revision 1.5  2008/02/13 13:07:56  roboos
 * fixed numerous bugs
 * implemented append function, which after all does not seem to be neccessary (since the problem is in the TCP buffer size)
 *
 * Revision 1.4  2008/02/11 21:16:23  roboos
 * after testing with a variety of messages, this version seems to be the first one that is fully functional
 * reimplemented many of the GET/PUT implementations based on the new header/data/event structures wheer each one ahs a def and a buf section
 * ensure that all pointers are set to NULL after a free()
 * print more debugging information
 * still to be done:
 * - add channel labels to header
 * - add channel selection to GET_DAT
 * - add event selection to GET_EVT
 * - implement mutexes around global memory access
 *
 * Revision 1.3  2008/02/10 20:28:05  roboos
 * implemented put_evt, get_dat
 * other small changes
 *
* Revision 1.2  2008/02/10 10:50:38  roboos
* made helper functions for freeing buffer
* made helper functions for initializing buffer
* implemented some commands
* at the moment probably not functional due to work-in-progress in GET_DAT, but for the rest it is tested and functional
*
* Revision 1.1  2008/02/09 16:44:26  roboos
* renamed network.c into socket.c
* made first skeleton implementation, seems partially functional
*
*/

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "unix_includes.h"
#ifdef ENABLE_POLLING
  #include <poll.h>
#endif

#include "socket_includes.h"

#include "buffer.h"
#include "message.h"

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

