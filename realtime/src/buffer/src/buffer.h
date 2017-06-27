/*
 * Copyright (C) 2008, Robert Oostenveld & Christian Hesse
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 */

#ifndef BUFFER_H
#define BUFFER_H

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "platform_includes.h"
#include "message.h"

#include "cleanup.h"
#include "endianutil.h"
#include "interface.h"
#include "message.h"
#include "printstruct.h"
#include "swapbytes.h"
#include "util.h"

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

#define BACKLOG           16
#define DEFAULT_HOSTNAME "localhost"
#define DEFAULT_PORT      1972
#define HOSTNAME_LENGTH   256

#define SO_RCVBUF_SIZE 16384
#define SO_SNDBUF_SIZE 16384

/* this is because the function has been renamed, but is perhaps already in use in other software */
#define open_remotehost open_connection

/* FIXME these should be variable */
#define MAXNUMBYTE      (512*1024*1024)
#define MAXNUMSAMPLE    600000
#define MAXNUMEVENT     100

#define WRAP(x,y) ((x) - ((int)((float)(x)/(y)))*(y))
#define FREE(x) {if (x) {free(x); x=NULL;}}
#define DIE_BAD_MALLOC(ptr) if ((ptr)==NULL) { fprintf(stderr, "Out of memory in line %d", __LINE__); exit(1); }

typedef struct {
	char name[HOSTNAME_LENGTH];
	int  port;
} host_t;

#ifdef __cplusplus
extern "C" {
#endif

	/* definition of the functions that implement the network transparent server */
	void *tcpserver(void *);
	void *tcpsocket(void *);
	int clientrequest(int, const message_t *, message_t**);
	int dmarequest(const message_t *, message_t**);
	int tcprequest(int, const message_t *, message_t**);

#ifdef __cplusplus
}
#endif

#endif
