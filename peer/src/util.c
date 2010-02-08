#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>       /* for strerror */

#include "peer.h"
#include "extern.h"
#include "platform_includes.h"

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

				if (verbose>0)
						fprintf(stderr, "bufread: read %d bytes\n", numthis);
				numread += numthis;
				numcall++;
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

				if (verbose)
						fprintf(stderr, "bufwrite: wrote %d bytes\n", numthis);
				numwrite += numthis;
				numcall++;
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


		if (((*buf1)!=NULL) && (bufsize1==0)) {
				perror("append err1");
				return -1;
		}
		else if (((*buf1)==NULL) && (bufsize1!=0)) {
				perror("append err2");
				return -1;
		}

		if ((*buf1)==NULL) {
				if (verbose>0)
						fprintf(stderr, "append: allocating %d bytes\n", bufsize2);
				(*buf1) = malloc(bufsize2);
		}
		else if ((*buf1)!=NULL) {
				if (verbose>0)
						fprintf(stderr, "append: reallocating from %d to %d bytes\n", bufsize1, bufsize1+bufsize2);
				(*buf1) = realloc((*buf1), bufsize1+bufsize2);
		}

		memcpy((char*)(*buf1)+bufsize1, buf2, bufsize2);
		return (bufsize1+bufsize2);
}

int close_connection(int s) {
		int status = 0, verbose = 0;
		if (verbose>0)
				fprintf(stderr, "close_connection: socket = %d\n", s);
		if (s>0)
				status = closesocket(s);	/* it is a TCP connection */
		if (status!=0)
				perror("close_connection");
		return status;
}

int open_connection(const char *hostname, int port) {
		int verbose = 0;
		int s, retry;
		struct sockaddr_in sa;
		struct hostent *host = NULL;

#ifdef WIN32
		WSADATA wsa;

		if(WSAStartup(MAKEWORD(1, 1), &wsa)) {
				fprintf(stderr, "open_connection: cannot start sockets\n");
				/* FIXME should this exception be handled more explicitely?  */
		}
#endif

		if (verbose>0)
				fprintf(stderr, "open_connection: server = %s, port = %d\n", hostname, port);

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
		sa.sin_port   = htons(port);
		memcpy(&(sa.sin_addr.s_addr), host->h_addr_list[0], sizeof(sa.sin_addr.s_addr));

		if ((s = socket(PF_INET, SOCK_STREAM, 0)) < 0) {
				if (verbose>0) fprintf(stderr, "open_connection: socket = %d\n", s);
				perror("open_connection");
				return -1;
		}

		retry = 3;
		while (retry>0) {
				if (connect(s, (struct sockaddr *)&sa, sizeof sa)<0) {
						/* wait 5 miliseconds and try again */
						perror("open_connection");
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
		   perror("open_connection connect");
		   usleep(1000000);
		   }
		 */

		if (verbose>0)
				fprintf(stderr, "open_connection: connected to %s:%d on socket %d\n", hostname, port, s);
		return s;
}

int jobcount(void) {
		int jobcount = 0;
		joblist_t *job;
		pthread_mutex_lock(&mutexjoblist);
		job = joblist;
		while (job) {
				jobcount++;
				job = job->next;
		}
		pthread_mutex_unlock(&mutexjoblist);
		return jobcount;
}

int peercount(void) {
		int peercount = 0;
		peerlist_t *peer;
		pthread_mutex_lock(&mutexpeerlist);
		peer = peerlist;
		while (peer) {
				peercount++;
				peer = peer->next;
		}
		pthread_mutex_unlock(&mutexpeerlist);
		return peercount;
}

int hoststatus(void) {
		int status = -1;
		pthread_mutex_lock(&mutexhost);
		if (host)
				status = host->status;
		pthread_mutex_unlock(&mutexhost);
		return status;
}

int ismember_userlist(char *str) {
  userlist_t *user = NULL;
  int ismember = 1;
  pthread_mutex_lock(&mutexuserlist);
  if (userlist) {
    ismember = 0;
    user = userlist;
    while (user) {
      if (strncmp(user->name, str, STRLEN)==0)
        ismember = 1;
      user = user->next;
    }
  }
  pthread_mutex_unlock(&mutexuserlist);
		return ismember;
}

int ismember_grouplist(char *str) {
		grouplist_t *group = NULL;
		int ismember = 1; 
		pthread_mutex_lock(&mutexgrouplist);
		if (grouplist) {
				ismember = 0;
				group = grouplist;
				while (group) {
						if (strncmp(group->name, str, STRLEN)==0)
								ismember = 1;
						group = group->next;
				}
		}
		pthread_mutex_unlock(&mutexgrouplist);
		return ismember;
}

int ismember_hostlist(char *str) {
		hostlist_t *host = NULL;
		int ismember = 1; 
		pthread_mutex_lock(&mutexhostlist);
		if (hostlist) {
				ismember = 0;
				host = hostlist;
				while (host) {
						if (strncmp(host->name, str, STRLEN)==0)
								ismember = 1;
						host = host->next;
				}
		}
		pthread_mutex_unlock(&mutexhostlist);
		return ismember;
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
}

