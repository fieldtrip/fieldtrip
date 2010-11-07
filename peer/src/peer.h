/* prevent double include */
#ifndef PEER_H
#define PEER_H

#include <pthread.h>
#include "platform_includes.h"
#if defined(PLATFORM_WIN32) || defined(PLATFORM_WIN64)
#include "win32/stdint.h"
#else
#include <stdint.h>
#include <syslog.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef POLLRDNORM
#define POLLRDNORM POLLIN
#endif

#ifndef POLLRDBAND
#define POLLRDBAND POLLPRI
#endif

#ifndef POLLWRNORM
#define POLLWRNORM POLLOUT
#endif

#ifndef POLLWRBAND
#define POLLWRBAND POLLOUT
#endif

#ifndef SO_REUSEPORT
/* SO_REUSEADDR is used on linux */
/* SO_REUSEPORT is used on osx   */
#define SO_REUSEPORT      SO_REUSEADDR
#endif

/* this is UNIX only, and may not work on all flavours */
#define USE_ABSTRACT_UDS_NAMES

#define STATUS_ZOMBIE            0			/* status = 0 means zombie mode, don't accept anything   */
#define STATUS_MASTER            1			/* status = 1 means master mode, accept everything       */
#define STATUS_IDLE              2			/* status = 2 means idle slave, accept only a single job */
#define STATUS_BUSY              3 			/* status = 3 means busy slave, don't accept a new job   */

#define VERSION                  0x0012
#define ANNOUNCE_GROUP           "225.0.0.88"
#define ANNOUNCE_PORT 	         1700		/* it will auto-increment if the port is not available */
#define DEFAULT_GROUP            "unknown"
#define DEFAULT_USER             "unknown"
#define DEFAULT_HOST             "localhost"
#define DEFAULT_PORT             1701
#define DEFAULT_STATUS           STATUS_ZOMBIE
#define DEFAULT_MEMAVAIL         UINT32_MAX
#define DEFAULT_CPUAVAIL         0
#define DEFAULT_TIMAVAIL         (24*3600)

#define ACCEPTSLEEP              10000		/* in usec */
#define ANNOUNCESLEEP            1000000	/* in usec */
#define ANNOUNCEJITTER           10000 		/* in usec */
#define EXPIRESLEEP              1500000	/* in usec, should be longer than ANNOUNCESLEEP+ANNOUNCEJITTER */
#define BACKLOG                  16
#define EXPIRATION               3			/* in sec  */
#define SMARTMEM_MINIMUM         104857600	/* 100 MB  */
#define SMARTSHARE_HISTORY       2			/* number if history items peer peer */
#define SMARTSHARE_PREVHOSTCOUNT 3			/* number of times that a host has to "knock" */
#define SMARTSHARE_TIMEOUT       5			/* in sec  */
#define SMARTSHARE_TIMER         3			/* idle time in seconds after which smartshare is disabled */
#define SMARTCPU_TOLERANCE       0.05		/* the ideal load of a computer is N+0.05, with N the number of CPUs */
#define SO_RCVBUF_SIZE           16384
#define SO_SNDBUF_SIZE           16384

#define MAXPWDSIZE		 		 16384
#define MAXPATHSIZE		 		 16384
#define MAXARGSIZE				 INT32_MAX
#define STRLEN			 		 128

#ifndef SYSLOG
#define SYSLOG 1 /* should be 1 for command-line peerslave, 0 or 3 for MATLAB */
#endif

#if   SYSLOG == 0
#define PANIC(format, args...)			{exit(-1);}
#define DEBUG(level, format, args...)	{ }
#elif SYSLOG == 1
#define PANIC(format, args...)			{syslog(LOG_ERR, format, ## args); exit(-1);}
#define DEBUG(level, format, args...)	{syslog(level, format, ## args);}
#elif SYSLOG == 2
#define PANIC(format, args...)			{fprintf(stderr, format"\n", ## args); exit(-1);}
#define DEBUG(level, format, args...)	{if (level<=syslog_level) fprintf(stderr, format"\n", ## args);}
#elif SYSLOG == 3
#define PANIC(format, args...)			{mexPrintf(format"\n", ## args); exit(-1);}
#define DEBUG(level, format, args...)	{if (level<=syslog_level) mexPrintf(format"\n", ## args);}
#endif

#define FREE(x)							{if (x) {free(x); x=NULL;}}

/* FIXME these are obvious at the moment, but should be formally defined */
typedef char      CHAR_T;
typedef float  FLOAT32_T;
typedef double FLOAT64_T;

/* the following types should be according to "ISO C99: 7.18 Integer types" (see /usr/include/stdint.h on OSX and Linux) */
/* FIXME different endianness between client/server is not supported at the moment */

#ifndef INT8_T
typedef   int8_t   INT8_T;
#endif

#ifndef INT16_T
typedef  int16_t  INT16_T;
#endif

#ifndef INT32_T
typedef  int32_t  INT32_T;
#endif

#ifndef INT64_T
typedef  int64_t  INT64_T;
#endif

#ifndef UINT8_T
typedef  uint8_t  UINT8_T;
#endif

#ifndef UINT16_T
typedef uint16_t UINT16_T;
#endif

#ifndef UINT32_T
typedef uint32_t UINT32_T;
#endif

#ifndef UINT64_T
typedef uint64_t UINT64_T;
#endif

/* these are used for type checking */
#define WORDSIZE_CHAR    sizeof(CHAR_T   )
#define WORDSIZE_UINT8   sizeof(UINT8_T  ) 
#define WORDSIZE_UINT16  sizeof(UINT16_T )
#define WORDSIZE_UINT32  sizeof(UINT32_T )
#define WORDSIZE_UINT64  sizeof(UINT64_T )
#define WORDSIZE_INT8    sizeof(INT8_T   )
#define WORDSIZE_INT16   sizeof(INT16_T  )
#define WORDSIZE_INT32   sizeof(INT32_T  )
#define WORDSIZE_INT64   sizeof(INT64_T  )
#define WORDSIZE_FLOAT32 sizeof(FLOAT32_T)
#define WORDSIZE_FLOAT64 sizeof(FLOAT64_T)

/* this describes the details of the current job (applies only to busy slave) */
typedef struct {
		UINT32_T pid;			/* UNIX process identifier of the peerslave */
		UINT32_T id;            /* identifier of the peer where the job originates from */
		char name[STRLEN];      /* hostname of the peer where the job originates from */
		char user[STRLEN];
		char group[STRLEN];
		UINT64_T timreq; 
		UINT64_T memreq; 
		UINT64_T cpureq; 
		UINT32_T argsize;   	/* size of the job arguments in bytes */
		UINT32_T optsize;   	/* size of the job options in bytes */
} current_t;

/* this is the binary packet that is announced over the network */
typedef struct {
		UINT32_T version;	   
		UINT32_T id;	   
		char name[STRLEN]; 
		char user[STRLEN]; 
		char group[STRLEN];
		char socket[STRLEN];  /* name of the unix domain socket, or empty if not available */
		UINT32_T port;
		UINT32_T status;
		UINT64_T timavail; 
		UINT64_T memavail; 
		UINT64_T cpuavail; 
		current_t current;
} hostdef_t;

typedef struct {
		UINT32_T version;
		UINT32_T id;
		UINT64_T timreq; 
		UINT64_T memreq; 
		UINT64_T cpureq; 
		UINT32_T argsize;   	/* size of the job arguments in bytes */
		UINT32_T optsize;   	/* size of the job options in bytes */
} jobdef_t;

typedef struct joblist_s {
		jobdef_t  *job;
		hostdef_t *host;
		void      *arg;
		void      *opt;
		struct joblist_s *next;
} joblist_t;

typedef struct peerlist_s {
		hostdef_t *host;
		time_t time;        /* time in seconds since January 1, 1970, Coordinated Universal Time */
		char ipaddr[INET_ADDRSTRLEN];
		struct peerlist_s *next;
} peerlist_t;

typedef struct smartsharelist_s {
		UINT64_T timreq; 
		struct smartsharelist_s *next;
} smartsharelist_t;

/* this is for restricting the access to a peer to one user or a list of users */
typedef struct userlist_s {
		char *name;
		struct userlist_s *next;
} userlist_t;

/* this is for restricting the access to a peer to one group or a list of groups */
typedef struct grouplist_s {
		char *name;
		struct grouplist_s *next;
} grouplist_t;

/* this is for restricting the access to a peer to one host or a list of hosts */
typedef struct hostlist_s {
		char *name;
		struct hostlist_s *next;
} hostlist_t;

typedef struct {
		hostdef_t *host; /* this defines the host details                 */
		jobdef_t  *job;  /* this defines the job contents                 */
		void      *arg;  /* this contains the input or output arguments   */
		void      *opt;  /* this contains the options for the job         */
} message_t;

#ifdef __cplusplus
extern "C" {
#endif

/* core function definitions */
void *udsserver (void *);
void *tcpserver (void *);
void *tcpsocket (void *);
void *announce  (void *);
void *discover  (void *);
void *expire    (void *);
void  peerinit (void *);
void  peerexit (void *);
int   announce_once(void);

/* functions from smartshare.c */
void smartshare_reset   (void);
int  smartshare_check   (float timreq, int hostid);
void smartshare_history (jobdef_t *job);

/* fnuctions from smartmem.c */
int smartmem_update(void);
int smartcpu_update(void);

/* functions from security.c */
int security_check(hostdef_t *host);
int ismember_userlist(char *);
int ismember_grouplist(char *);
int ismember_hostlist(char *);

/* functions from util.c and elsewhere */
int  append(void **buf1, int bufsize1, void *buf2, int bufsize2);
int  bufread(int s, void *buf, int numel);
int  bufwrite(int s, void *buf, int numel);
int  close_connection(int s);
int  hoststatus(void);
int  jobcount(void);
int  open_tcp_connection(const char *hostname, int port);
int  open_uds_connection(const char *socketname);
int  peercount(void);
int  check_localhost(const char *ipaddr);
int  check_killswitch(void);
void check_datatypes(void);
void clear_grouplist(void);
void clear_hostlist(void);
void clear_joblist(void);
void clear_peerlist(void);
void clear_smartsharelist(void);
void clear_userlist(void);

#ifdef __cplusplus
}
#endif

#endif /* PEER_H */
