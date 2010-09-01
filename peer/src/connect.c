/*
 * Copyright (C) 2010, Robert Oostenveld
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/
 *
 */

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>       /* for strerror */

#include "peer.h"
#include "extern.h"
#include "platform_includes.h"

int open_uds_connection(const char *socketname) {
#if defined (PLATFORM_WIN32) || defined (PLATFORM_WIN64)
		/* not yet implemented */
#else
		int s, len, verbose = 0;
		struct sockaddr_un remote;

		if ((s = socket(AF_UNIX, SOCK_STREAM, 0)) == -1) {
				perror("open_uds_connection socket");
				return -1;
		}

		remote.sun_family = AF_UNIX;
		strcpy(remote.sun_path, socketname);
		len = strlen(remote.sun_path) + sizeof(remote.sun_family);
#ifdef USE_ABSTRACT_UDS_NAMES
		remote.sun_path[0] = 0;
#endif
		if (connect(s, (struct sockaddr *)&remote, len) == -1) {
				perror("open_uds_connection connect");
				return -1;
		}

		pthread_mutex_lock(&mutexconnectioncount);
		connectioncount++;
		pthread_mutex_unlock(&mutexconnectioncount);

		if (verbose>0)
				fprintf(stderr, "open_uds_connection: connected to %s on socket %d\n", socketname, s);
		return s;
#endif
}


int open_tcp_connection(const char *hostname, int port) {
		int verbose = 0;
		int s, retry;
		struct sockaddr_in sa;
		struct hostent *host = NULL;

#ifdef WIN32
		WSADATA wsa;
#endif

		if (port==0) {
				if (verbose>0)
						fprintf(stderr, "open_tcp_connection: using direct memory copy\n");
				return 0;
		}
		else {
				if (verbose>0)
						fprintf(stderr, "open_tcp_connection: server = %s, port = %d\n", hostname, port);
		}

#ifdef WIN32
		if(WSAStartup(MAKEWORD(1, 1), &wsa)) {
				fprintf(stderr, "open_tcp_connection: cannot start sockets\n");
				/* FIXME should this exception be handled more explicitely?  */
		}
#endif

		if ((host = gethostbyname(hostname)) == NULL) {
				fprintf(stderr, "open_tcp_connection: nslookup1 failed on '%s'\n", hostname);
				return -1;
		}

		if (host->h_length == 0) {
				fprintf(stderr, "open_tcp_connection: nslookup2 failed on '%s'\n", hostname);
				return -1;
		}

		bzero(&sa, sizeof sa);
		sa.sin_family = AF_INET;
		sa.sin_port   = htons(port);
		memcpy(&(sa.sin_addr.s_addr), host->h_addr_list[0], sizeof(sa.sin_addr.s_addr));

		s = socket(PF_INET, SOCK_STREAM, 0);
		if (verbose>0)
				fprintf(stderr, "open_tcp_connection: socket = %d\n", s);
		if (s<0) {
				perror("open_tcp_connection");
				return -1;
		}

		retry = 3;
		while (retry>0) {
				if (connect(s, (struct sockaddr *)&sa, sizeof sa)<0) {
						/* wait 5 miliseconds and try again */
						perror("open_tcp_connection");
						usleep(5000);
						retry--;
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
		   perror("open_tcp_connection connect");
		   usleep(1000000);
		   }
		 */

		pthread_mutex_lock(&mutexconnectioncount);
		connectioncount++;
		pthread_mutex_unlock(&mutexconnectioncount);

		if (verbose>1)
				fprintf(stderr, "open_tcp_connection: connectioncount = %d\n", connectioncount);

		if (verbose>0)
				fprintf(stderr, "open_tcopen_tcpnnection: connected to %s:%d on socket %d\n", hostname, port, s);
		return s;
}

int close_connection(int s) {
		int retval, verbose = 0;

		if (verbose>0)
				fprintf(stderr, "close_connection: socket = %d\n", s);

		if (s==0)
				/* for a DMA connection */
				retval = 0;	
		else if (s>0)
				/* for a TPC or UDS socket connection */
				retval = closesocket(s);
		else
				/* unkown connection */
				retval = -1;

		if (retval<0)
				perror("close_connection");

		pthread_mutex_lock(&mutexconnectioncount);
		connectioncount--;
		pthread_mutex_unlock(&mutexconnectioncount);

		if (verbose>1)
				fprintf(stderr, "close_connection: connectioncount = %d\n", connectioncount);

		return retval;
}

