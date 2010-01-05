/*
 * Copyright (C) 2008, Robert Oostenveld & Christian Hesse
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 * $Log: util.c,v $
 * Revision 1.33  2009/06/29 10:48:02  roboos
 * retry opening the connection 10 times
 *
 * Revision 1.32  2009/01/23 19:47:50  roboos
 * changed verbosity
 *
 * Revision 1.31  2009/01/23 08:26:44  roboos
 * fixed a serious bug that caused a lot of memory to leak (in fact all packets that were sent over the socket would eventually leack away), both on the client and server side
 *
 * Revision 1.30  2009/01/21 20:54:35  roboos
 * added some debugging code to keep track of number of sockets and threads (both seem ok)
 * cleaned up the debugging: more consistent use of verbose flag and function name in fprintf
 *
 * Revision 1.29  2008/12/16 20:38:23  roboos
 * removed message_t, header_t, data_t and event_t from check_datatypes()
 * function, because they contain a void pointer which varies on linux32
 * and linux64, and these types are not used themselves in the communication
 * (i.e. never sent over the network)
 *
 * Revision 1.28  2008/11/14 15:48:34  roboos
 * number of small changes, nothing significant
 *
 * Revision 1.27  2008/10/29 20:12:27  roboos
 * renamed open_remotehost into open_connection and added function close_connection
 * added some fprintf statements, no functional changes
 *
 * Revision 1.26  2008/07/09 10:37:52  roboos
 * moved some printing functions to seperate file
 *
 * Revision 1.25  2008/06/20 07:53:52  roboos
 * addec check on sizeof(xxx_t) to ensure correct packed representation of network transparent structures
 *
 * Revision 1.24  2008/05/29 11:55:52  roboos
 * small change in debugging feedback only
 *
 * Revision 1.23  2008/05/29 10:33:16  roboos
 * open_remotehost: moved WSA declaration to top, now also works on borland
 *
 * Revision 1.22  2008/05/29 07:55:00  roboos
 * moved checks for the packing of structs and wordsize to seperate function
 *
 * Revision 1.21  2008/05/22 09:55:16  roboos
 * some changes for compatibility wioth Borland, thanks to Jurgen
 *
 * Revision 1.19  2008/04/24 15:41:59  roboos
 * changed verbosity
 *
 * Revision 1.18  2008/04/14 14:13:37  thohar
 * moved win32 socket init code
 *
 * Revision 1.17  2008/03/26 14:37:37  thohar
 * bufread and bufwrite now also break when connection has been closed
 *
 * Revision 1.16  2008/03/23 17:01:12  roboos
 * removed polling code, which was already commented out
 *
 * Revision 1.15  2008/03/17 13:47:55  roboos
 * print_buf also detects NULL and prints debug info accordingly
 *
 * Revision 1.14  2008/03/13 13:45:55  roboos
 * removed open_localhost
 * some other small changes
 *
 * Revision 1.13  2008/03/13 12:38:01  thohar
 * added headers for win32 build
 * removed polling as it is not neccessary
 * replaced read with recv for portability
 * replaced write with send for portability
 * bufread does not exit app now when recv return < 0
 * open_remothost now starts the sockets on win32 build
 *
 * Revision 1.12  2008/03/10 09:44:27  roboos
 * added property functions
 *
 * Revision 1.11  2008/03/07 14:51:06  roboos
* more strict input checks in append
* added inteligence in open_remotehost so that "localhost:0" will result in direct memcpy
*
* Revision 1.10  2008/03/02 13:22:34  roboos
* added polling for closed connection to bufread add write
*
* Revision 1.9  2008/02/27 10:13:46  roboos
* added print_buf
*
* Revision 1.8  2008/02/26 21:43:26  roboos
* renamed packet_t structure definition into messagedef_t, added message_t structure (contains def+buf)
*
* Revision 1.7  2008/02/20 13:39:45  roboos
* modified append function
* deleted append3 function
* do a retry if no connection can be made (pause 1 sec)
*
* Revision 1.6  2008/02/19 17:24:47  roboos
* also print fsample
*
* Revision 1.5  2008/02/19 10:22:56  roboos
* added consistent copyrights and log message to each of the files
*
*
*/

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>       /* for strerror */

#include "socket_includes.h"
#include "unix_includes.h"

#include "buffer.h"
#include "message.h"

/* these are for debugging */
pthread_mutex_t mutexappendcount = PTHREAD_MUTEX_INITIALIZER;
int appendcount = 0;
pthread_mutex_t mutexthreadcount = PTHREAD_MUTEX_INITIALIZER;
int threadcount = 0;
pthread_mutex_t mutexsocketcount = PTHREAD_MUTEX_INITIALIZER;
int socketcount = 0;

int bufread(int s, void *buf, int numel) {
		int numcall = 0, numthis = 0, numread = 0, verbose = 0;

		while (numread<numel) {

				numthis = recv(s, (char*)buf+numread, numel-numread, 0);
				if (numthis<0) {
						perror("bufread");
						break;
				}
				else if (numthis == 0)
						break;

				if (verbose>0) fprintf(stderr, "bufread: read %d bytes\n", numthis);
				numread += numthis;
				numcall ++;
				usleep(1000);
		}
		if (verbose>1)
				fprintf(stderr, "bufread: reading the complete buffer required %d calls\n", numcall);
		return numread;
}

int bufwrite(int s, void *buf, int numel) {
		int numcall = 0, numthis = 0, numwrite = 0, verbose = 0;

		while (numwrite<numel) {

				numthis = send(s, (char*)buf+numwrite, numel-numwrite, 0);
				if (numthis<0) {
						perror("bufwrite");
						break;
				}
				else if(numthis == 0)
						break;

				if (verbose) fprintf(stderr, "bufwrite: wrote %d bytes\n", numthis);
				numwrite += numthis;
				numcall ++;
				usleep(1000);
		}
		if (verbose>1)
				fprintf(stderr, "bufwrite: writing the complete buffer required %d calls\n", numcall);
		return numwrite;
}

int append(void **buf1, int bufsize1, void *buf2, int bufsize2) {
		int verbose = 0;

		if (verbose>1) {
				pthread_mutex_lock(&mutexappendcount);
				appendcount++;
				fprintf(stderr, "append: appendcount = %d\n", appendcount);
				pthread_mutex_unlock(&mutexappendcount);
		}

		if (((*buf1)==NULL) && (!bufsize1)) {
				if (verbose>0) fprintf(stderr, "append: allocating %d bytes\n", bufsize1+bufsize2);
				(*buf1) = malloc(bufsize1+bufsize2);
		}
		else if (((*buf1)!=NULL) && (bufsize1)) {
				if (verbose>0) fprintf(stderr, "append: reallocating from %d to %d bytes\n", bufsize1, bufsize1+bufsize2);
				(*buf1) = realloc((*buf1), bufsize1+bufsize2);
		}
		else {
				perror("append err1");
				exit(1);
		}
		if ((*buf1)==NULL) {
				perror("append err2");
				exit(1);
		}
		memcpy((char*)(*buf1)+bufsize1, buf2, bufsize2);
		return (bufsize1+bufsize2);
}

int close_connection(int s) {
		int status = 0, verbose = 0;
		if (verbose>0) fprintf(stderr, "close_connection: socket = %d\n", s);
		if (s>0)
				status = closesocket(s);	/* it is a TCP connection */
		if (status!=0)
				perror("close_connection");
		return status;
}

int open_connection(const char *hostname, int port) {
		int verbose = 0;
        int retry = 10;
		int s;
		struct sockaddr_in sa;
		struct hostent *host;

#ifdef WIN32
		WSADATA wsa;
#endif

		if (port==0) {
				if (verbose>0) fprintf(stderr, "open_connection: using direct memory copy\n");
				return 0;
		}
		else {
				if (verbose>0) fprintf(stderr, "open_connection: server = %s, port = %d\n", hostname, port);
		}

#ifdef WIN32
		if(WSAStartup(MAKEWORD(1, 1), &wsa)) {
				fprintf(stderr, "open_connection: cannot start sockets\n");
				/* FIXME should this exception be handled more explicitely?  */
		}
#endif

		if ((s = socket(PF_INET, SOCK_STREAM, 0)) < 0) {
				if (verbose>0) fprintf(stderr, "open_connection: socket = %d\n", s);
				perror("open_connection");
				return -1;
		}

		if ((host = gethostbyname(hostname)) == NULL) {
				fprintf(stderr, "open_connection: nslookup1 failed on '%s'\n", hostname);
				return -1;
		}

		if (host->h_length == 0) {
				fprintf(stderr, "open_connection: nslookup2 failed on '%s'\n", hostname);
				return -1;
		}

		bzero(&sa, sizeof sa);
		sa.sin_family = AF_INET;
		sa.sin_port = htons(port);
		memcpy(&(sa.sin_addr.s_addr), host->h_addr_list[0], sizeof(sa.sin_addr.s_addr));

		while (retry>0) {
			if (connect(s, (struct sockaddr *)&sa, sizeof sa)<0) {
					perror("open_connection");
					usleep(5000); /* wait 5 miliseconds */
					retry = retry--;
			}
			else {
				/* this signals that the connection has been made */
				retry = -1;
			}
		}
		if (retry==0) {
			/* it failed on mutliple attempts, give up */
			return -2;
		}

		/*
		   while (connect(s, (struct sockaddr *)&sa, sizeof sa) < 0) {
		   perror("open_connection connect");
		   usleep(1000000);
		   }
		 */

		if (verbose>0) fprintf(stderr, "open_connection: connected to %s:%d on socket %d\n", hostname, port, s);
		return s;
}

void check_datatypes() {
		/* check datatypes */
		if (WORDSIZE_CHAR    !=1) { fprintf(stderr, "invalid size of CHAR    (%d)\n", WORDSIZE_CHAR   ); exit(-1); }
		if (WORDSIZE_UINT8   !=1) { fprintf(stderr, "invalid size of UINT8   (%d)\n", WORDSIZE_UINT8  ); exit(-1); }
		if (WORDSIZE_UINT16  !=2) { fprintf(stderr, "invalid size of UINT16  (%d)\n", WORDSIZE_UINT16 ); exit(-1); }
		if (WORDSIZE_UINT32  !=4) { fprintf(stderr, "invalid size of UINT32  (%d)\n", WORDSIZE_UINT32 ); exit(-1); }
		if (WORDSIZE_UINT64  !=8) { fprintf(stderr, "invalid size of UINT64  (%d)\n", WORDSIZE_UINT64 ); exit(-1); }
		if (WORDSIZE_INT8    !=1) { fprintf(stderr, "invalid size of INT8    (%d)\n", WORDSIZE_INT8   ); exit(-1); }
		if (WORDSIZE_INT16   !=2) { fprintf(stderr, "invalid size of INT16   (%d)\n", WORDSIZE_INT16  ); exit(-1); }
		if (WORDSIZE_INT32   !=4) { fprintf(stderr, "invalid size of INT32   (%d)\n", WORDSIZE_INT32  ); exit(-1); }
		if (WORDSIZE_INT64   !=8) { fprintf(stderr, "invalid size of INT64   (%d)\n", WORDSIZE_INT64  ); exit(-1); }
		if (WORDSIZE_FLOAT32 !=4) { fprintf(stderr, "invalid size of FLOAT32 (%d)\n", WORDSIZE_FLOAT32); exit(-1); }
		if (WORDSIZE_FLOAT64 !=8) { fprintf(stderr, "invalid size of FLOAT64 (%d)\n", WORDSIZE_FLOAT64); exit(-1); }
		if (sizeof(messagedef_t)  !=8 ) { fprintf(stderr, "invalid size of messagedef_t\n"); exit(-1); }
		if (sizeof(headerdef_t)   !=24) { fprintf(stderr, "invalid size of headerdef_t \n"); exit(-1); }
		if (sizeof(propertydef_t) !=20) { fprintf(stderr, "invalid size of propertydef_t \n"); exit(-1); }
		if (sizeof(datadef_t)     !=16) { fprintf(stderr, "invalid size of datadef_t   \n"); exit(-1); }
		if (sizeof(eventdef_t)    !=32) { fprintf(stderr, "invalid size of eventdef_t  \n"); exit(-1); }
		if (sizeof(datasel_t)     !=8 ) { fprintf(stderr, "invalid size of datasel_t   \n"); exit(-1); }
		if (sizeof(eventsel_t)    !=8 ) { fprintf(stderr, "invalid size of eventsel_t  \n"); exit(-1); }
}

