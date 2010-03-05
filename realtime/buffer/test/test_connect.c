/*
 * Use as
 *    ./test_connect localhost 1972
 *
 * Copyright (C) 2008, Christian Hesse & Robert Oostenveld
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include "buffer.h"

int main(int argc, char *argv[]) {
		host_t host;
		int server = -1;

		/* these variables are for writing the data */
		int status = 0, verbose = 0, numconnect, numloop;
		struct timeval tic, toc;
		float elapsed;

		/* start with defaults */
		sprintf(host.name, DEFAULT_HOSTNAME);
		host.port = DEFAULT_PORT;

		if (argc>1)
				sprintf(host.name, argv[1]);

		if (argc>2)
				host.port = atoi(argv[2]);

		gettimeofday(&tic, NULL);

		numconnect = 0;
		numloop    = 0;

		while (1) {
				numconnect++;
				numloop++;

				server = open_connection(host.name, host.port);
				status = close_connection(server);

				if (numloop==1000) {
						gettimeofday(&toc, NULL);
						elapsed = 1.0 * (toc.tv_sec-tic.tv_sec) + 0.000001 * (toc.tv_usec-tic.tv_usec);
						fprintf(stderr, "numconnect = %d, elapsed = %f, numconnect/sec = %f\n", numconnect, elapsed, ((float)(numconnect))/((float)elapsed));
						numloop = 0;
				}

		} /* while(1) */

		return 0;
}

