/*
 * Copyright (C) 2008, Robert Oostenveld & Christian Hesse
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>       /* for strerror */

#include "buffer.h"
#include <pthread.h>
#include "extern.h"

unsigned int bufread(int s, void *buf, unsigned int numel) {
		unsigned int numcall = 0, numread = 0, verbose = 0;
		int numthis = 0;
		
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
				numcall ++;
#ifndef PLATFORM_WIN32			/* SK: I think this shouldn't be necessary on any platform: the sockets are blocking */
				if (numread<numel)
						usleep(1000);
#endif
		}
		if (verbose>1)
				fprintf(stderr, "bufread: reading the complete buffer required %d calls\n", numcall);
		return numread;
}

unsigned int bufwrite(int s, const void *buf, unsigned int numel) {
		int numthis = 0;
		unsigned int numcall = 0, numwrite = 0, verbose = 0;

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
				numcall ++;
#ifndef PLATFORM_WIN32			/* SK: I think this shouldn't be necessary on any platform: the sockets are blocking */				
				if (numwrite<numel)
						usleep(1000);
#endif
		}
		if (verbose>1)
				fprintf(stderr, "bufwrite: writing the complete buffer required %d calls\n", numcall);
		return numwrite;
}

unsigned int append(void **buf1, unsigned int bufsize1, void *buf2, unsigned int bufsize2) {
		int verbose = 0;

		if (verbose>1) {
				pthread_mutex_lock(&mutexappendcount);
				appendcount++;
				fprintf(stderr, "append: appendcount = %d\n", appendcount);
				pthread_mutex_unlock(&mutexappendcount);
		}

		if (((*buf1)!=NULL) && (bufsize1==0)) {
				perror("append err1");
				return 0;	/* was -1, but this is never checked anyway */
		}
		else if (((*buf1)==NULL) && (bufsize1!=0)) {
				perror("append err2");
				return 0;	/* was -1, but this is never checked anyway */
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
		struct hostent *host;
#ifdef WIN32
		static WSADATA wsa = {0,0}; /* check version fields to only initialise once */
#endif

		if (port==0) {
				if (verbose>0)
						fprintf(stderr, "open_connection: using direct memory copy\n");
				return 0;
		}
		else {
				if (verbose>0)
						fprintf(stderr, "open_connection: server = %s, port = %d\n", hostname, port);
		}

#ifdef WIN32
		if (wsa.wVersion == 0) {
			/* We only need to do this once ... and actually have a corresponding WSACleanup call somewhere */
			if(WSAStartup(MAKEWORD(1, 1), &wsa)) {
				fprintf(stderr, "open_connection: cannot start sockets\n");
				/* FIXME should this exception be handled more explicitely?  */
			}
		}
#endif

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

		if ((s = socket(PF_INET, SOCK_STREAM, 0)) < 0) {
				if (verbose>0)
						fprintf(stderr, "open_connection: socket = %d\n", s);
				perror("open_connection");
				return -1;
		}
    
		retry = 10;
		while (retry>0) {
				if (connect(s, (struct sockaddr *)&sa, sizeof sa)<0) {
						/* wait 5 miliseconds and try again */
						usleep(5000);
						retry--;
				}
				else {
						/* this signals that the connection has been made */
						retry = -1;
				}
		}
		if (retry==0) {
				/* close the socket */
				closesocket(s);
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

#ifdef DISABLE_NAGLE
		{
			int optval = 1;
			setsockopt(s, IPPROTO_TCP, TCP_NODELAY, &optval, sizeof(optval));
		}
#endif
        
		return s;
}

void check_datatypes() {
		/* check datatypes */
		if (WORDSIZE_CHAR    !=1) { fprintf(stderr, "invalid size of CHAR    (%d)\n", (int) WORDSIZE_CHAR   ); exit(-1); }
		if (WORDSIZE_UINT8   !=1) { fprintf(stderr, "invalid size of UINT8   (%d)\n", (int) WORDSIZE_UINT8  ); exit(-1); }
		if (WORDSIZE_UINT16  !=2) { fprintf(stderr, "invalid size of UINT16  (%d)\n", (int) WORDSIZE_UINT16 ); exit(-1); }
		if (WORDSIZE_UINT32  !=4) { fprintf(stderr, "invalid size of UINT32  (%d)\n", (int) WORDSIZE_UINT32 ); exit(-1); }
		if (WORDSIZE_UINT64  !=8) { fprintf(stderr, "invalid size of UINT64  (%d)\n", (int) WORDSIZE_UINT64 ); exit(-1); }
		if (WORDSIZE_INT8    !=1) { fprintf(stderr, "invalid size of INT8    (%d)\n", (int) WORDSIZE_INT8   ); exit(-1); }
		if (WORDSIZE_INT16   !=2) { fprintf(stderr, "invalid size of INT16   (%d)\n", (int) WORDSIZE_INT16  ); exit(-1); }
		if (WORDSIZE_INT32   !=4) { fprintf(stderr, "invalid size of INT32   (%d)\n", (int) WORDSIZE_INT32  ); exit(-1); }
		if (WORDSIZE_INT64   !=8) { fprintf(stderr, "invalid size of INT64   (%d)\n", (int) WORDSIZE_INT64  ); exit(-1); }
		if (WORDSIZE_FLOAT32 !=4) { fprintf(stderr, "invalid size of FLOAT32 (%d)\n", (int) WORDSIZE_FLOAT32); exit(-1); }
		if (WORDSIZE_FLOAT64 !=8) { fprintf(stderr, "invalid size of FLOAT64 (%d)\n", (int) WORDSIZE_FLOAT64); exit(-1); }
		if (sizeof(messagedef_t)  !=8 ) { fprintf(stderr, "invalid size of messagedef_t\n"); exit(-1); }
		if (sizeof(headerdef_t)   !=24) { fprintf(stderr, "invalid size of headerdef_t \n"); exit(-1); }
		if (sizeof(datadef_t)     !=16) { fprintf(stderr, "invalid size of datadef_t   \n"); exit(-1); }
		if (sizeof(eventdef_t)    !=32) { fprintf(stderr, "invalid size of eventdef_t  \n"); exit(-1); }
		if (sizeof(datasel_t)     !=8 ) { fprintf(stderr, "invalid size of datasel_t   \n"); exit(-1); }
		if (sizeof(eventsel_t)    !=8 ) { fprintf(stderr, "invalid size of eventsel_t  \n"); exit(-1); }
}



unsigned int wordsize_from_type(UINT32_T data_type) {
	switch(data_type) {
		case DATATYPE_CHAR:
			return WORDSIZE_CHAR;
		case DATATYPE_UINT8:		
		case DATATYPE_INT8:
			return WORDSIZE_INT8;
		case DATATYPE_UINT16:
		case DATATYPE_INT16:
			return WORDSIZE_INT16;
		case DATATYPE_UINT32:				
		case DATATYPE_INT32:
			return WORDSIZE_INT32;
		case DATATYPE_UINT64:
		case DATATYPE_INT64:
			return WORDSIZE_INT64;
		case DATATYPE_FLOAT32:
			return WORDSIZE_FLOAT32;
		case DATATYPE_FLOAT64:
			return WORDSIZE_FLOAT64;
	}
	return 0;
}

const ft_chunk_t *find_chunk(const void *buf, unsigned int offset0, unsigned int size, UINT32_T chunk_type) {
	unsigned int bufpos = offset0;
	while (bufpos + sizeof(ft_chunkdef_t) <= size) {
		const ft_chunk_t *chunk = (ft_chunk_t *) ((char *)buf + bufpos);
		if (chunk->def.type == chunk_type) return chunk;
		
		bufpos += sizeof(ft_chunkdef_t) + chunk->def.size;
	}
	return NULL;
}


/** Iterate through an array of events and check whether all of them are properly defined,
	that is, whether the "type" and "value" fields are of valid type and size, and whether the
	"bufsize" fields are correct (that is, fully contained in the passed buffer, and big enough
	to hold "type" and "value".
	Returns the number of events on success (might also be 0), or
	a negative number that indicates in which event definition an error happend,
	for example, a return value of -2 means that the first event was ok, but the second event
	definition was invalid. This function returns at the first error.
*/
int check_event_array(unsigned int size, const void *buf) {
	unsigned int offset=0;
	int numEvents = 0;
	
	while (offset + sizeof(eventdef_t) <= size) {
		const eventdef_t *E;	
		unsigned int wsType, wsValue;
		
		/* Set our event pointer to the current location within the array */
		E = (const eventdef_t *) ((char *) buf + offset);
		
		/* Increase the offset by the size of this event, and check whether it's fully
		   contained within the given array.
		*/
		offset += sizeof(eventdef_t) + E->bufsize;
		if (offset > size) goto error;
		
		/* Check whether "type" and "value" are of known type */
		wsType = wordsize_from_type(E->type_type);
		if (wsType == 0) goto error;
		wsValue = wordsize_from_type(E->value_type);
		if (wsValue == 0) goto error;
		
		/* Check whether "type" and "value" are contained in this event's "buf" */
		if (wsType * E->type_numel + wsValue * E->value_numel > E->bufsize) goto error;

		/* all checks passed, continue at next offset */
		++numEvents;
	}
	return numEvents;
error:
	return -(1+numEvents);
}





#ifdef WIN32
int open_unix_connection(const char *name) {
	return -1;
}
#else
int open_unix_connection(const char *name) {
	int verbose = 0;
	int s, retry;
	struct sockaddr_un sa;

	bzero(&sa, sizeof(sa));
	sa.sun_family = AF_UNIX;
	strncpy(sa.sun_path, name, sizeof(sa.sun_path));
		
	s = socket(AF_UNIX, SOCK_STREAM, 0);
	if (s < 0) {
		perror("open_unix_connection, socket");
		return -1;
	}

	retry = 10;
	while (retry>0) {
		if (connect(s, (struct sockaddr *)&sa, sizeof(sa))<0) {
			/* wait 5 miliseconds and try again */
			perror("open_connection");
			usleep(5000);
			retry--;
		} else {
			/* this signals that the connection has been made */
			retry = -1;
			}
		}
	if (retry==0) {
		/* it failed on mutliple attempts, give up */
		return -2;
	}

	if (verbose>0)
		fprintf(stderr, "open_unix_connection: connected to %s on socket %d\n", name,  s);
        
	return s;
}
#endif



#ifdef WIN32
#ifndef COMPILER_MINGW

/*
 * timeval.h    1.0 01/12/19
 *
 * Defines gettimeofday, timeval, etc. for Win32
 *
 * By Wu Yongwei
 *
 */

//#define EPOCHFILETIME (116444736000000000i64)
#define EPOCHFILETIME ((INT64_T) 116444736000000000LL)

#ifdef COMPILER_LCC
VOID STDCALL GetSystemTimeAsFileTime(LPFILETIME);
#endif

int gettimeofday(struct timeval *tv, struct timezone *tz)
{
    FILETIME        ft;
    LARGE_INTEGER   li;
    INT64_T         t;
    static int      tzflag;
    
    if (tv) {
        GetSystemTimeAsFileTime(&ft);
        li.LowPart  = ft.dwLowDateTime;
        li.HighPart = ft.dwHighDateTime;
        t  = li.QuadPart;       /* In 100-nanosecond intervals */
        t -= EPOCHFILETIME;     /* Offset to the Epoch time */
        t /= 10;                /* In microseconds */
        tv->tv_sec  = (long)(t / 1000000);
        tv->tv_usec = (long)(t % 1000000);
    }
    
	#ifndef COMPILER_LCC
	/* LCC that comes with Matlab has problems with _timezone and _daylight, 
		and we don't need it anyway */
    if (tz) {
        if (!tzflag) {
            _tzset();
            tzflag++;
        }
        tz->tz_minuteswest = _timezone / 60;
        tz->tz_dsttime = _daylight;
    }
	#endif

    return 0;
}

#endif
#endif
