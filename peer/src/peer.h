/* prevent double include */
#ifndef PEER_H
#define PEER_H

#include <pthread.h>
#include "platform_includes.h"
#if defined(PLATFORM_WIN32) || defined(PLATFORM_WIN64)
#include "win32/stdint.h"
#else
#include <stdint.h>
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

#define VERSION           0x0009
#define ANNOUNCE_GROUP    "225.0.0.88"
#define ANNOUNCE_PORT 	  1700				/* it will auto-increment if the port is not available */
#define DEFAULT_GROUP     "unknown"
#define DEFAULT_USER      "unknown"
#define DEFAULT_HOST      "localhost"
#define DEFAULT_PORT      1701
#define DEFAULT_STATUS    0
#define DEFAULT_MEMAVAIL  UINT32_MAX
#define DEFAULT_CPUAVAIL  UINT32_MAX
#define DEFAULT_TIMAVAIL  UINT32_MAX
#define FAIRSHARE_PREVHOSTCOUNT 3	/* number of times that a host has to "knock" */
#define FAIRSHARE_TIMER         3	/* idle time in seconds after which fairshare is disabled */
#define FAIRSHARE_HISTORY       2	/* number if history items peer peer */
#define EXPIRATION        5
#define BACKLOG           16
#define SO_RCVBUF_SIZE    16384
#define SO_SNDBUF_SIZE    16384
#define ACCEPTSLEEP       10000    /* in usec */
#define ANNOUNCESLEEP     1000000  /* in usec */
#define EXPIRESLEEP       1000000  /* in usec */
#define FAIRSHARE_TIMEOUT 5        /* in sec  */

#define MAXPWDSIZE		  16384
#define MAXPATHSIZE		  16384
#define MAXARGSIZE		  INT32_MAX
#define STRLEN			  256

#define FREE(x) {if (x) {free(x); x=NULL;}}

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

typedef struct {
		UINT32_T version;
		UINT32_T id;
		char name[STRLEN];
		char user[STRLEN];
		char group[STRLEN];
		UINT32_T port;
		UINT32_T status;
		UINT64_T timavail; 
		UINT64_T memavail; 
		UINT64_T cpuavail; 
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

typedef struct fairsharelist_s {
		UINT32_T timreq; 
		struct fairsharelist_s *next;
} fairsharelist_t;

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
void *tcpserver (void *);
void *tcpsocket (void *);
void *announce  (void *);
void *discover  (void *);
void *expire    (void *);
void  peerinit (void *);
void  peerexit (void *);
int   announce_once(void);

/* functions from fairshare.c */
void fairshare_reset   (void);
int  fairshare_check   (float timreq, int hostid);
void fairshare_history (jobdef_t *job);

/* fnuctions from smartmem.c */
int meminfo(UINT64_T *MemTotal, UINT64_T *MemFree);
int smartmem_update(void);
UINT64_T smartmem_avail(void);

/* functions from security.c */
int security_check(hostdef_t *host);
int ismember_userlist(char *);
int ismember_grouplist(char *);
int ismember_hostlist(char *);

/* functions from util.c */
int  append(void **buf1, int bufsize1, void *buf2, int bufsize2);
int  bufread(int s, void *buf, int numel);
int  bufwrite(int s, void *buf, int numel);
int  open_connection(const char *hostname, int port);
int  close_connection(int s);
void check_datatypes(void);
int  jobcount(void);
int  peercount(void);
int  hoststatus(void);
int  localhost(const char *ipaddr);
void clear_peerlist(void);
void clear_joblist(void);
void clear_userlist(void);
void clear_grouplist(void);
void clear_hostlist(void);
void clear_fairsharelist(void);

#ifdef __cplusplus
}
#endif

#endif /* PEER_H */
