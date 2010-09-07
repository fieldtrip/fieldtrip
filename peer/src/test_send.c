#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "peer.h"
#include "extern.h"

int main(int argc, char *argv[]) {
		int verbose = 0, server, i;
		peerlist_t *peer;

		jobdef_t *def;
		void *buf;

		/* these variables are for the threading */
		int rc;
		pthread_t tid;

		peerinit(NULL);

		if (verbose>2) {
				pthread_mutex_lock(&mutexhost);
				fprintf(stderr, "test_send: host.name =  %s\n", host->name);
				fprintf(stderr, "test_send: host.port =  %u\n", host->port);
				pthread_mutex_unlock(&mutexhost);
		}

		if ((rc = pthread_create(&tid, NULL, discover, (void *)NULL))>0) {
				/* the code should never arrive here */
				if (verbose>0) fprintf(stderr, "peer: return code from pthread_create() is %d\n", rc);
				goto cleanup;
		}

		if ((rc = pthread_create(&tid, NULL, expire, (void *)NULL))>0) {
				/* the code should never arrive here */
				if (verbose>0) fprintf(stderr, "peer: return code from pthread_create() is %d\n", rc);
				goto cleanup;
		}

		/* sleep for two seconds */
		sleep(2);

		def = (jobdef_t *)malloc(sizeof(jobdef_t));
		def->version = VERSION;
		def->id      = random();
		def->argsize = 256;
		def->optsize = 256;

		buf = malloc(256);
		for (i=0; i<256; i++)
				((UINT8_T *)buf)[i] = i;

		/* send all available peers a message */
		while (1) {
				pthread_mutex_lock(&mutexpeerlist);
				peer = peerlist;
				i = 0;
				while(peer) {
						fprintf(stderr, "-------------------------------------------------\n");
						fprintf(stderr, "test_send: peerlist[%d].host.name = %s\n", i, peer->host->name);
						fprintf(stderr, "test_send: peerlist[%d].host.port = %u\n", i, peer->host->port);

						/* open the TCP socket */
						server = open_tcp_connection(peer->host->name, peer->host->port);

						if (server<0) {
								fprintf(stderr, "test_send: failed to create socket\n");
						}
						else {
								write(server, def, sizeof(jobdef_t));
								pthread_mutex_lock(&mutexhost);
								write(server, host, sizeof(hostdef_t));
								pthread_mutex_unlock(&mutexhost);
								write(server, buf, 256); /* for the job arguments */
								write(server, buf, 256); /* for the job options   */
								close(server);
						}

						peer = peer->next ;
						i++;

				} /* while (peer) */
				pthread_mutex_unlock(&mutexpeerlist);
		} /* while (1) */

		peerexit(NULL);
		exit(0);

cleanup:
		fprintf(stderr, "test_send: cleanup\n");
		exit(1);

} /* main */
