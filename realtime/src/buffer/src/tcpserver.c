/*
 * Copyright (C) 2008, Robert Oostenveld & Christian Hesse
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* these are for changing the socket to non-blocking mode */
#include <fcntl.h>
#include <errno.h>

#include <pthread.h>
#include "buffer.h"
#include "extern.h"

#define ACCEPTSLEEP 1000

typedef struct {
  int fd;
} threadlocal_t;

void cleanup_tcpserver(void *arg) {
  threadlocal_t *threadlocal;
  threadlocal = (threadlocal_t *)arg;
  if (threadlocal && threadlocal->fd>0) {
    closesocket(threadlocal->fd);
    threadlocal->fd = -1;
  }

  pthread_mutex_lock(&mutexstatus);
  tcpserverStatus = 0;
  pthread_mutex_unlock(&mutexstatus);
}

/* pthread_attr_t attr; */ /* this one would be passed to the thread */

/***********************************************************************
 * this thread listens to incoming TCP connections
 * if a connection is made by a client, it starts the tcpsocket function
 ***********************************************************************/
void *tcpserver(void *arg) {
  int verbose = 0;
  host_t *host;

  /* these variables are for the socket */
  struct sockaddr_in sa;
  int s, c;
  int b;
  int optval;
  /* struct timeval timeout; */
  int oldcancelstate, oldcanceltype;

#if defined(PLATFORM_WIN32) || defined(PLATFORM_WIN64)
  unsigned long enable = 0;
  static WSADATA wsa = {0,0};
#endif

  /* these variables are for the threading */
  int rc;
  pthread_t tid;

  threadlocal_t threadlocal;
  threadlocal.fd = -1;

  /* this determines the port on which the server will listen */
  if (!arg)
    exit(1);
  else
    host = (host_t *)arg;

  if (verbose>0) fprintf(stderr, "tcpserver: host.name =  %s\n", host->name);
  if (verbose>0) fprintf(stderr, "tcpserver: host.port =  %d\n", host->port);

  pthread_cleanup_push(cleanup_tcpserver, &threadlocal);

  /* the status contains the thread id when running, or zero when not running */
  pthread_mutex_lock(&mutexstatus);
  if (tcpserverStatus==0) {
    tcpserverStatus = 1;
    pthread_mutex_unlock(&mutexstatus);
  }
  else {
    pthread_mutex_unlock(&mutexstatus);
    goto cleanup;
  }

#if defined(PLATFORM_WIN32) || defined(PLATFORM_WIN64)
  /* We only need to do this once ... and actually have a corresponding WSACleanup call somewhere */
  if (wsa.wVersion == 0) {
    if(WSAStartup(MAKEWORD(1, 1), &wsa)) {
      if (verbose>0) fprintf(stderr, "tcpserver: cannot start sockets\n");
    }
  }
#endif

  /* setup socket */
  if ((s = socket(PF_INET, SOCK_STREAM, 0)) < 0) {
    perror("tcpserver socket");
    goto cleanup;
  }

  /* this will be closed at cleanup */
  threadlocal.fd = s;

  /* place the socket in non-blocking mode, required to do thread cancelation */
#if defined(PLATFORM_WIN32) || defined(PLATFORM_WIN64)
  enable = 0;
  ioctlsocket(s, FIONBIO, &enable);
#else
  optval = fcntl(s, F_GETFL, NULL);
  optval = optval | O_NONBLOCK;
  if (fcntl(s, F_SETFL, optval)<0) {
    perror("tcpserver fcntl");
    goto cleanup;
  }
#endif

  /* change the receive timeout */
  /*
     timeout.tv_sec  = 1;
     timeout.tv_usec = 1;
     if (setsockopt(s, SOL_SOCKET, SO_RCVTIMEO, &timeout, sizeof(optval)) < 0) {
     perror("tcpserver setsockopt");
     if (verbose>0) fprintf(stderr, "tcpserver: errno = %d\n", errno);
     goto cleanup;
     }
   */

  /* prevend "bind: address already in use" */
  optval = 1;
  if (setsockopt(s, SOL_SOCKET, SO_REUSEADDR, (const char*)&optval, sizeof(optval)) < 0) {
    perror("tcpserver setsockopt");
    goto cleanup;
  }

  bzero(&sa, sizeof sa);
  sa.sin_family = AF_INET;
  sa.sin_port   = htons(host->port);

  if (INADDR_ANY)
    sa.sin_addr.s_addr = htonl(INADDR_ANY);

  if (bind(s, (struct sockaddr *)&sa, sizeof sa) < 0) {
    perror("tcpserver bind");
    goto cleanup;
  }

  if (listen(s, BACKLOG)<0) {
    perror("tcpserver listen");
    goto cleanup;
  }

  for (;;) {
    /*
     * If no pending connections are present on the queue, and the socket
     * is not marked as non-blocking, accept() blocks the caller until a
     * connection is present.  If the socket is marked non-blocking and
     * no pending connections are present on the queue, accept() returns
     * an error as described below.
     */

    b = sizeof sa;
    c = accept(s, (struct sockaddr *)&sa, &b);

    if (verbose>1)
      fprintf(stderr, "tcpserver: accept\n");

    if (c<0) {
#if defined(PLATFORM_WIN32) || defined(PLATFORM_WIN64)
      if(errno == 0) {
        pthread_testcancel();
        usleep(ACCEPTSLEEP);
      }
      else {
        perror("tcpserver accept");
        goto cleanup;
      }
#else
      if (errno==EWOULDBLOCK) {
        pthread_testcancel();
        usleep(ACCEPTSLEEP);
      }
      else {
        perror("tcpserver accept");
        goto cleanup;
      }
#endif
    }

    else {
      if (verbose>0) fprintf(stderr, "tcpserver: opened connection to client on socket %d\n", c);

      /* SK: we could set socket option "SO_LINGER" so resources are freed up immediately
         but we leave this at the default for now
       */
      if (0) {
        struct linger lg;

        lg.l_onoff = 1;
        lg.l_linger = 0;

        setsockopt(s, SOL_SOCKET, SO_LINGER, (char *) &lg, sizeof(lg));
      }


      /* set larger buffer */
      optval = SO_RCVBUF_SIZE;
      if (setsockopt(c, SOL_SOCKET, SO_RCVBUF, (const char*)&optval, sizeof(optval)) < 0) {
        perror("tcpserver setsockopt");
        goto cleanup;
      }

      /* set larger buffer */
      optval = SO_SNDBUF_SIZE;
      if (setsockopt(c, SOL_SOCKET, SO_SNDBUF, (const char*)&optval, sizeof(optval)) < 0) {
        perror("tcpserver setsockopt");
        goto cleanup;
      }

      /* place the socket back in blocking mode, this is needed for tcpsocket  */
#if defined(PLATFORM_WIN32) || defined(PLATFORM_WIN64) 
      enable = 0;
      ioctlsocket(c, FIONBIO, &enable);
#else
      optval = fcntl(c, F_GETFL, NULL);
      optval = optval & (!O_NONBLOCK);
      if (fcntl(c, F_SETFL, optval)<0) {
        perror("tcpserver fcntl");
        goto cleanup;
      }
#endif

#ifdef DISABLE_NAGLE
      /* disable the Nagle buffering algorithm */
      optval = 1;
      if (setsockopt(c, IPPROTO_TCP, TCP_NODELAY, &optval, sizeof(optval)) < 0) {
        perror("tcpserver setsockopt");
        goto cleanup;
      }
#endif

      /* deal with the incoming connection on the TCP socket in a seperate thread */
      /* rc = pthread_create(&tid, &attr, tcpsocket, (void *)c); */
      rc = pthread_create(&tid, NULL, tcpsocket, (void *)c);

      if (rc) {
        if (verbose>0) fprintf(stderr, "tcpserver: return code from pthread_create() is %d\n", rc);
        goto cleanup;
      }
      else {
        /* this is for debugging */
        pthread_mutex_lock(&mutexthreadcount);
        threadcount++;
        pthread_mutex_unlock(&mutexthreadcount);
        if (verbose>0) fprintf(stderr, "tcpserver: c = %d, threadcount = %d\n", c, threadcount);
        pthread_detach(tid);
      }
    }
  }

cleanup:
  printf("");
  pthread_cleanup_pop(1);
  return NULL;
}
