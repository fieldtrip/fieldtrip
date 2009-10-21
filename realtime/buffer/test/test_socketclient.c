#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <strings.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>

#include "buffer.h"

void main(int argc, char *argv[]) {
	int server;
	int n = 0, status;
	void *buf = NULL;
	time_t t0, t1, t2;

	/* open the TCP socket */
	if ((server = open_connection(argv[1], atoi(argv[2]))) < 0) {
		fprintf(stderr, "ERROR; failed to create socket\n");
		exit(1);
	}

	time(&t0);
	status = 1;
	while (status) {
		time(&t1);
		n++;
		buf = malloc(BUFSIZE);
		status = (bufwrite(server, buf, BUFSIZE)==BUFSIZE);
		FREE(buf);
		time(&t2);
		fprintf(stderr, "n = %d, ct = %d, dt = %d\n", n, t2-t0, t2-t1);
	}

	close(server);
}
