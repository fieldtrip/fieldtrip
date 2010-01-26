
/* prevent double include */
#ifndef PEER_H
#define PEER_H

#include <pthread.h>
#include <stdint.h>
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

#define VERSION           0x0005
#define ANNOUNCE_GROUP    "225.0.0.88"
#define ANNOUNCE_PORT 	  1700				/* it will auto-increment if the port is not available */
#define DEFAULT_GROUP     "all"
#define DEFAULT_USER      "unknown"
#define DEFAULT_HOST      "localhost"
#define DEFAULT_PORT      1701
#define DEFAULT_STATUS    0
#define EXPIRATION        5
#define BACKLOG           16
#define SO_RCVBUF_SIZE    16384
#define SO_SNDBUF_SIZE    16384
#define ACCEPTSLEEP       10000    /* in usec */
#define ANNOUNCESLEEP     1000000  /* in usec */
#define EXPIRESLEEP       1000000  /* in usec */

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
		char name[STRLEN];
		char addr[STRLEN];
		char user[STRLEN];
		char group[STRLEN];
		UINT32_T id;
		UINT32_T port;
		UINT32_T status;
} hostdef_t;

typedef struct {
		UINT32_T version;
		UINT32_T id;
		UINT32_T argsize;   /* size of the job arguments in bytes */
		UINT32_T optsize;   /* size of the job options in bytes */
} jobdef_t;

struct joblist_s {
		jobdef_t  *job;
		hostdef_t *host;
		void      *arg;
		void      *opt;
		struct joblist_s *next;
};

struct peerlist_s {
		hostdef_t *host;
		time_t time;        /* time in seconds since January 1, 1970, Coordinated Universal Time */
		struct peerlist_s *next;
};

/* this is for restricting the access to a peer to one user or a list of users */
struct userlist_s {
  char *name;
  struct userlist_s *next;
};

/* this is for restricting the access to a peer to one group or a list of groups */
struct grouplist_s {
  char *name;
  struct grouplist_s *next;
};

/* this is for restricting the access to a peer to one host or a list of hosts */
struct hostlist_s {
  char *name;
  struct hostlist_s *next;
};

typedef struct joblist_s joblist_t;
typedef struct peerlist_s peerlist_t;
typedef struct userlist_s userlist_t;
typedef struct grouplist_s grouplist_t;
typedef struct hostlist_s hostlist_t;

typedef struct {
		jobdef_t  *job;  /* this defines the job contents                 */
		hostdef_t *host; /* this defines where it originates from         */
		void      *arg;  /* this contains the input or output arguments   */
		void      *opt;  /* this contains the options for the job         */
} message_t;

#ifdef __cplusplus
extern "C" {
#endif

		/* function definitions */
		void *tcpserver (void *);
		void *tcpsocket (void *);
		void *announce  (void *);
		void *discover  (void *);
		void *expire    (void *);
		void  peerinit (void *);
		void  peerexit (void *);

		/* functions from util.c */
		int bufread(int s, void *buf, int numel);
		int bufwrite(int s, void *buf, int numel);
		int append(void **buf1, int bufsize1, void *buf2, int bufsize2);
		int close_connection(int s);
		int open_connection(const char *hostname, int port);
		void check_datatypes(void);
		int jobcount(void);
		int peercount(void);
		int hoststatus(void);
		int ismember_userlist(char *);
		int ismember_grouplist(char *);
		int ismember_hostlist(char *);

#ifdef __cplusplus
}
#endif

#endif /* PEER_H */
