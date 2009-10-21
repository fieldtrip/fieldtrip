/*
 * Copyright (C) 2008, Christian Hesse & Robert Oostenveld
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 * $Log: demo_event.c,v $
 * Revision 1.2  2008/07/09 10:40:08  roboos
 * some minor cleanup
 *
 * Revision 1.1  2008/07/09 10:08:25  roboos
 * renamed the various demo applications into demo_xxx
 *
 * Revision 1.1  2008/07/08 20:23:23  roboos
 * new demo code: based on old code from sinewave demo, but more simple
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "message.h"
#include "buffer.h"
#include "unix_includes.h"
#include "socket_includes.h"

int main(int argc, char *argv[]) {
	host_t host;
    check_datatypes();

	if (argc>2) {
		sprintf(host.name, argv[1]);
		host.port = atoi(argv[2]);
	}
	else {
		sprintf(host.name, DEFAULT_HOSTNAME);
		host.port = DEFAULT_PORT;
	}

    fprintf(stderr, "demo_event: host.name =  %s\n", host.name);
    fprintf(stderr, "demo_event: host.port =  %d\n", host.port);

	/* start the acquisition */
	event_thread((void *)(&host));

	return 0;
}

