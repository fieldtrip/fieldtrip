#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "peer.h"
#include "extern.h"

int main(int argc, char *argv[]) {
		int verbose = 0;

		/* these variables are for the threading */
		int rc;
		pthread_t tid;

		peerinit(NULL);

		sleep(1);
		host->status = 2;
		sleep(1);

		if ((rc = pthread_create(&tid, NULL, announce, (void *)&host))>0) {
				/* the code should never arrive here */
				if (verbose>0) fprintf(stderr, "test_recv: return code from pthread_create() is %d\n", rc);
				goto cleanup;
		}

		if ((rc = pthread_create(&tid, NULL, tcpserver, (void *)&host))>0) {
				/* the code should never arrive here */
				if (verbose>0) fprintf(stderr, "test_recv: return code from pthread_create() is %d\n", rc);
				goto cleanup;
		}

		for (;;) {
				sleep(1);
		}

		peerexit(NULL);
		exit(0);

cleanup:
		fprintf(stderr, "test_recv: cleanup\n");
		exit(1);
}

