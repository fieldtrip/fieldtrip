/*
 *
 * Use as
 *    ./test_nslookup <hostname>
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include "buffer.h"

int main(int argc, char *argv[]) {
		int server = -1;
		struct hostent *host;

		/* these variables are for writing the data */
		int status = 0, verbose = 0, numlookup, numloop;
		struct timeval tic, toc;
		float elapsed;
		
		#ifdef WIN32
		WSADATA wsaData;
		if(WSAStartup(MAKEWORD(1, 1), &wsaData)) {
			fprintf(stderr, "Cannot start Windows sockets.\n");
			return 1;
		}
		#endif

		gettimeofday(&tic, NULL);

		numlookup = 0;
		numloop    = 0;

		while (1) {
				numlookup++;
				numloop++;

				if ((host = gethostbyname(argv[1])) == NULL) {
						fprintf(stderr, "nslookup1 failed on '%s'\n", argv[1]);
						return -1;
				}

				if (host->h_length == 0) {
						fprintf(stderr, "nslookup2 failed on '%s'\n", argv[1]);
						return -1;
				}

				if (numloop==10000) {
						gettimeofday(&toc, NULL);
						elapsed = 1.0 * (toc.tv_sec-tic.tv_sec) + 0.000001 * (toc.tv_usec-tic.tv_usec);
						fprintf(stderr, "numlookup = %d, elapsed = %f, numlookup/sec = %f\n", numlookup, elapsed, ((float)(numlookup))/((float)elapsed));
						numloop = 0;
				}

		} /* while(1) */

		return 0;
}


