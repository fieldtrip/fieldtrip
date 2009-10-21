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
	/* these variables are for the socket */
	struct sockaddr_in sa;
	int s, c, b, n = 0;
	int optval, status;
	void *buf = NULL;

	/* setup socket */
	if ((s = socket(PF_INET, SOCK_STREAM, 0)) < 0) {
		perror("socket");
		return 1;
	}

	/* prevend "bind: address already in use" */
	optval = 1;
	if (setsockopt(s, SOL_SOCKET, SO_REUSEADDR, &optval, sizeof(optval)) < 0) {
		perror("setsockopt");
		return 1;
	}

	bzero(&sa, sizeof sa);

	sa.sin_family = AF_INET;
	sa.sin_port   = htons(atoi(argv[1]));

	if (INADDR_ANY)
		sa.sin_addr.s_addr = htonl(INADDR_ANY);

	if (bind(s, (struct sockaddr *)&sa, sizeof sa) < 0) {
		perror("bind");
		return 2;
	}

	listen(s, BACKLOG);

	for (;;) {
		b = sizeof sa;

		if ((c = accept(s, (struct sockaddr *)&sa, &b)) < 0) {
			perror("buffer accept");
			return 4;
		}
		else {
			fprintf(stderr, "opened connection to client on socket %d (%d)\n", c, n);

			/* set larger buffer */
			optval = SO_RCVBUF_SIZE;
			if (setsockopt(c, SOL_SOCKET, SO_RCVBUF, &optval, sizeof(optval)) < 0) {
				perror("setsockopt");
				return 1;
			}

			/* set larger buffer */
			optval = SO_SNDBUF_SIZE;
			if (setsockopt(c, SOL_SOCKET, SO_SNDBUF, &optval, sizeof(optval)) < 0) {
				perror("setsockopt");
				return 1;
			}

			/* disable the Nagle buffering algorithm */
			/*
			optval = 1;
			if (setsockopt(c, IPPROTO_TCP, TCP_NODELAY, &optval, sizeof(optval)) < 0) {
				perror("setsockopt");
				return 1;
			}
			*/

			status = 1;
			while (status) {
				fprintf(stderr, "reading from socket\n");
				buf = malloc(BUFSIZE);
				status = (bufread(c, buf, BUFSIZE)==BUFSIZE);
				FREE(buf);
			}

		}
	}

	return 0;
}
