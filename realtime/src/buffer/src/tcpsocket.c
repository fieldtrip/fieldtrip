/*
 * Copyright (C) 2008, Robert Oostenveld
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>

#include "buffer.h"
#include "extern.h"

#ifdef ENABLE_POLLING
#include <poll.h>
#endif

#define THREADSLEEP      1000000  /* in microseconds */
#define POLLSLEEP        100      /* in microseconds */
#define MERGE_THRESHOLD  4096     /* TODO: optimize this value? Maybe look at MTU size */

typedef struct {
  void *message;
  int fd;
} threadlocal_t;

void cleanup_tcpsocket(void *arg) {
  threadlocal_t *threadlocal;
  threadlocal = (threadlocal_t *)arg;
  if (threadlocal && threadlocal->message) {
    FREE(threadlocal->message);
  }
  if (threadlocal && threadlocal->fd>0) {
    closesocket(threadlocal->fd);
    threadlocal->fd = -1;
  }

  pthread_mutex_lock(&mutexsocketcount);
  socketcount--;
  pthread_mutex_unlock(&mutexsocketcount);

  pthread_mutex_lock(&mutexthreadcount);
  threadcount--;
  pthread_mutex_unlock(&mutexthreadcount);
}

/* this function deals with the incoming client request */
void *tcpsocket(void *arg) {
  int n;
  int status = 0, verbose = 0, noresponse = 0;

#ifdef ENABLE_POLLING
  struct pollfd fds;
#endif

  /* these are used for communication over the TCP socket */
  int client = 0;
  message_t *request = NULL, *response = NULL;

  threadlocal_t threadlocal;
  threadlocal.message = NULL;
  threadlocal.fd = -1;

  /* the connection to the client has been made by the server */
  client = (int)arg;

  /* this will be closed at cleanup */
  threadlocal.fd = client;

  pthread_cleanup_push(cleanup_tcpsocket, &threadlocal);

  /* this is for debugging */
  pthread_mutex_lock(&mutexsocketcount);
  socketcount++;
  pthread_mutex_unlock(&mutexsocketcount);

  if (verbose>1) fprintf(stderr, "tcpsocket: client = %d, socketcount = %d, threadcount = %d\n", client, socketcount, threadcount);

  /* keep processing messages untill the connection is closed */
  while (1) {
    int swap = 0;
    UINT16_T reqCommand;
    UINT32_T respBufSize;

    request       = (message_t*)malloc(sizeof(message_t));
    DIE_BAD_MALLOC(request);

    request->def  = (messagedef_t*)malloc(sizeof(messagedef_t));
    DIE_BAD_MALLOC(request->def);
    request->buf  = NULL;

#ifdef ENABLE_POLLING
    /* wait for data to become available or until the connection is closed */
    /* thohar: i think this is not necessary as we dont need a timeout. */
    /* roboos: we need it to detect when the socket is closed by the client */
    while (1) {
      fds.fd      = client;
      fds.events  = POLLIN | POLLRDNORM | POLLRDBAND | POLLPRI | POLLOUT | POLLWRNORM | POLLWRBAND | POLLERR | POLLNVAL;
      fds.revents = 0;

      if (poll(&fds, 1, 1)==-1) {
        perror("poll");
        goto cleanup;
      }

      if (fds.revents & POLLHUP)
        goto cleanup; /* the connection has been closed */
      else if (fds.revents & POLLERR)
        goto cleanup; /* the connection has been closed */
      else if (fds.revents & POLLIN)
        break; /* data is available, process the message */
      else
        usleep(POLLSLEEP); /* wait for data or closed connection */
    }
#endif

    if ((n = bufread(client, request->def, sizeof(messagedef_t))) != sizeof(messagedef_t)) {
      if (verbose>0) fprintf(stderr, "tcpsocket: packet size = %d, should be %lu\n", n, sizeof(messagedef_t));
      goto cleanup;
    }

    if (request->def->version==VERSION_OE) {
      swap = 1;
      ft_swap16(2, &request->def->version); /* version + command */
      ft_swap32(1, &request->def->bufsize);
      reqCommand = request->def->command;
    }

    if (request->def->version!=VERSION) {
      if (verbose>0) fprintf(stderr, "tcpsocket: incorrect request version\n");
      goto cleanup;
    }

    if (request->def->bufsize>0) {
      request->buf = malloc(request->def->bufsize);
      DIE_BAD_MALLOC(request->buf);
      if ((n = bufread(client, request->buf, request->def->bufsize)) != request->def->bufsize) {
        if (verbose>0) fprintf(stderr, "tcpsocket: read size = %d, should be %d\n", n, request->def->bufsize);
        goto cleanup;
      }
    }

    if (swap && request->def->bufsize > 0) ft_swap_buf_to_native(reqCommand, request->def->bufsize, request->buf);

    if (verbose>1) print_request(request->def);
    if (verbose>1) print_buf(request->buf, request->def->bufsize);

    if ((status = dmarequest(request, &response)) != 0) {
      if (verbose>0) fprintf(stderr, "tcpsocket: an unexpected error occurred\n");
      goto cleanup;
    }

    DIE_BAD_MALLOC(response);
    DIE_BAD_MALLOC(response->def);

    if (verbose>1) print_response(response->def);
    if (verbose>1) print_buf(request->buf, request->def->bufsize);

    respBufSize = response->def->bufsize;
    if (swap) ft_swap_from_native(reqCommand, response);

    /* check whether a response is requested */
    switch (request->def->command) {
      case PUT_HDR_NORESPONSE:
      case PUT_DAT_NORESPONSE:
      case PUT_EVT_NORESPONSE:
        noresponse = 1;
        break;
      default:
        noresponse = 0;
    }

    /* we don't need the request anymore */
    cleanup_message(&request);
    request = NULL;

    if (!noresponse) {
      /* the request should be answered with a response */

      if (respBufSize + sizeof(messagedef_t) <= MERGE_THRESHOLD) {
        /* merge response->def and response->buf if they are small, so we can send it in one go over TCP */
        int msize = respBufSize + sizeof(messagedef_t);
        void *merged = NULL;

        append(&merged, 0, response->def, sizeof(messagedef_t));
        DIE_BAD_MALLOC(merged);
        append(&merged, sizeof(messagedef_t), response->buf, respBufSize);
        DIE_BAD_MALLOC(merged);

        if ((n=bufwrite(client, merged, msize) != msize)) {
          if (verbose>0) fprintf(stderr, "tcpsocket: write size = %d, should be %d\n", n, msize);
          FREE(merged);
          goto cleanup;
        }
        FREE(merged);
      } else {
        if ((n = bufwrite(client, response->def, sizeof(messagedef_t)))!=sizeof(messagedef_t)) {
          if (verbose>0) fprintf(stderr, "tcpsocket: write size = %d, should be %lu\n", n, sizeof(messagedef_t));
          goto cleanup;
        }
        if ((n = bufwrite(client, response->buf, respBufSize))!=respBufSize) {
          if (verbose>0) fprintf(stderr, "tcpsocket: write size = %d, should be %u\n", n, respBufSize);
          goto cleanup;
        }
      }
    }

    cleanup_message(&response);
    response = NULL;

  } /* while (1) */

cleanup:
  printf(""); /* otherwise the pthread_cleanup_pop won't compile */

  if (response!=NULL)
    cleanup_message(&response);
  response = NULL; /* SK: prevent double free in following pthread_cleanup_pop */

  pthread_cleanup_pop(1);

  /* this is for debugging */
  pthread_mutex_lock(&mutexsocketcount);
  socketcount--;
  pthread_mutex_unlock(&mutexsocketcount);

  /* this is for debugging */
  pthread_mutex_lock(&mutexthreadcount);
  threadcount--;
  pthread_mutex_unlock(&mutexthreadcount);

  pthread_exit(NULL);
  return NULL;
}
