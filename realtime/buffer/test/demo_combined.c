/*
 * demo - a multithreaded and networked buffer for EEG/MEG data and events
 *
 * Copyright (C) 2008, Christian Hesse & Robert Oostenveld
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 * $Log: demo_combined.c,v $
 * Revision 1.3  2009/01/23 17:50:32  roboos
 * removed some unnecessary code
 *
 * Revision 1.2  2008/12/16 20:35:55  roboos
 * changed c++ style comments into ansi c
 *
 * Revision 1.1  2008/07/09 10:08:25  roboos
 * renamed the various demo applications into demo_xxx
 *
 * Revision 1.1  2008/07/08 18:44:48  release
 * moved from src to test, renamed into sinewave_buffer
 *
 * Revision 1.7  2008/05/29 07:54:59  roboos
 * moved checks for the packing of structs and wordsize to seperate function
 *
 * Revision 1.6  2008/05/22 13:37:28  roboos
 * some small changes to make it compile with borland
 *
 * Revision 1.5  2008/03/13 13:34:48  roboos
 * use varargin for host and port or defauls
 *
 * Revision 1.4  2008/03/10 10:04:57  roboos
 * use host argument in making threads
 *
 * Revision 1.3  2008/03/08 10:33:48  roboos
 * changed some code to reflect the renamed low-level functions, no functional changes
 *
 * Revision 1.2  2008/02/26 21:43:25  roboos
 * renamed packet_t structure definition into messagedef_t, added message_t structure (contains def+buf)
 *
 * Revision 1.1  2008/02/21 12:10:24  roboos
 * this used to be the "demo" application
 *
 * Revision 1.4  2008/02/19 17:26:20  roboos
 * some small changes, nothing functional
 *
 * Revision 1.3  2008/02/19 10:22:56  roboos
 * added consistent copyrights and log message to each of the files
 *
 * Revision 1.2  2008/02/18 17:04:08  roboos
 * lots of small changes, debugging for brainamp
 *
 * Revision 1.1  2008/02/18 12:13:46  roboos
 * moved executable from buffer to demo
 * fixed bugs in sinewave and socket for events
 * stripped down the eventdef_t fields
 * many small changes
 *
 */

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "buffer.h"
#include "message.h"
#include "unix_includes.h"
#include "socket_includes.h"

/***********************************************************************
  this is the main thread
 ***********************************************************************/
int main(int argc, char *argv[]) {
	host_t host;

	/* these variables are for the threading */
	int rc;
	pthread_t tid;

    check_datatypes();

	if (argc>2) {
		sprintf(host.name, argv[1]);
		host.port = atoi(argv[2]);
	}
	else {
		sprintf(host.name, DEFAULT_HOSTNAME);
		host.port = DEFAULT_PORT;
	}

	/* start the buffer in a seperate thread */
	rc = pthread_create(&tid, NULL, tcpserver, (void *)(&host));
	if (rc) {
		fprintf(stderr, "main err1: return code from pthread_create() is %d\n", rc);
		exit(-1);
	}
    pthread_detach(tid);
	usleep(1000000);

	/* start the acquisition */
    sinewave_thread((void *)(&host));

	return 0;
}

