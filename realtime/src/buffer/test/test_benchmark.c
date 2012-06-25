/*
 * Use as
 *    ./test_benchmark localhost 1972 32 512 1
 *
 * Copyright (C) 2008, Christian Hesse & Robert Oostenveld
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "buffer.h"
#include <pthread.h>

int main(int argc, char *argv[]) {
  host_t host;
  
    /* these variables are for the threading */
  int rc;
  pthread_t tid;
  
    /* these variables are for writing the data */
  int i, j, k, status = 0, verbose = 0, samplecount = 0;
  time_t tic, toc;
  time_t elapsed;
  
  /* these represent the acquisition system properties */
  int fsample        = 512; /* this is just for the header */
  int nchans         = 32;
  int nsamples       = 64;
  int stateless      = 0;	  /* boolean */
  
  /* these are used in the communication and represent statefull information */
  int server             = -1;
  message_t    *request  = NULL;
  message_t    *response = NULL;
  header_t     *header   = NULL;
  data_t       *data     = NULL;
  event_t      *event    = NULL;
    
  
    /* start with defaults */
  sprintf(host.name, DEFAULT_HOSTNAME);
  host.port = DEFAULT_PORT;
  
  if (argc>1)
    sprintf(host.name, argv[1]);
  
  if (argc>2)
    host.port = atoi(argv[2]);
  
  if (argc>3)
    nchans = atoi(argv[3]);
  
  if (argc>4)
    nsamples = atoi(argv[4]);
  
  if (argc>5)
    stateless = atoi(argv[5]);
  

  check_datatypes();
  
  if (verbose>0) fprintf(stderr, "test_benchmark: host.name =  %s\n", host.name);
  if (verbose>0) fprintf(stderr, "test_benchmark: host.port =  %d\n", host.port);

  /* allocate the elements that will be used in the communication */
  request      = malloc(sizeof(message_t));
  request->def = malloc(sizeof(messagedef_t));
  request->buf = NULL;
  request->def->version = VERSION;
  request->def->bufsize = 0;
  
  header      = malloc(sizeof(header_t));
  header->def = malloc(sizeof(headerdef_t));
  header->buf = NULL;
  
  data      = malloc(sizeof(data_t));
  data->def = malloc(sizeof(datadef_t));
  data->buf = NULL;
  
  event      = malloc(sizeof(event_t));
  event->def = malloc(sizeof(eventdef_t));
  event->buf = NULL;
  
  /* define the header */
  header->def->nchans    = nchans;
  header->def->nsamples  = 0;
  header->def->nevents   = 0;
  header->def->fsample   = fsample;
  header->def->data_type = DATATYPE_FLOAT32;
  header->def->bufsize   = 0;
  FREE(header->buf);
  
  /* define the constant part of the data and allocate space for the variable part */
  data->def->nchans    = nchans;
  data->def->nsamples  = nsamples;
  data->def->data_type = DATATYPE_FLOAT32;
  data->def->bufsize   = WORDSIZE_FLOAT32*nchans*nsamples;
  FREE(data->buf);
  data->buf            = malloc(WORDSIZE_FLOAT32*nchans*nsamples);
    
    /* create the random data */
    for (j=0; j<nsamples; j++)
      for (i=0; i<nchans;    i++)
        ((FLOAT32_T *)(data->buf))[j*nchans+i] = 2.0*((FLOAT32_T)rand())/RAND_MAX - 1.0;
  
  /* initialization phase, send the header */
  request->def->command = PUT_HDR;
  request->def->bufsize = append(&request->buf, request->def->bufsize, header->def, sizeof(headerdef_t));
  request->def->bufsize = append(&request->buf, request->def->bufsize, header->buf, header->def->bufsize);
  
  server = open_connection(host.name, host.port);
  status = clientrequest(server, request, &response);
  if (verbose>0) fprintf(stderr, "test_benchmark: clientrequest returned %d\n", status);
  if (status) {
    fprintf(stderr, "random: err1\n");
    goto cleanup;
  }
  
  if (stateless) {
    status = close_connection(server);
    if (status) {
      fprintf(stderr, "random: err2\n");
      goto cleanup;
    }
  }
  
  cleanup_message(&request);
  cleanup_message(&response);
  request = NULL;
  response = NULL;
  
  tic = time(NULL);
  
  while (1) {
    
    /* create the request */
    request      = malloc(sizeof(message_t));
    request->def = malloc(sizeof(messagedef_t));
    request->buf = NULL;
    request->def->version = VERSION;
    request->def->bufsize = 0;
    request->def->command = PUT_DAT;
    request->def->bufsize = append(&request->buf, request->def->bufsize, data->def, sizeof(datadef_t));
    request->def->bufsize = append(&request->buf, request->def->bufsize, data->buf, data->def->bufsize);
    
    if (stateless)
      server = open_connection(host.name, host.port);
    
    status = clientrequest(server, request, &response);
    if (verbose>0) fprintf(stderr, "random: clientrequest returned %d\n", status);
    if (status) {
      fprintf(stderr, "random: err3\n");
      goto cleanup;
    }

    if (stateless) {
      status = close_connection(server);
      if (status) {
        fprintf(stderr, "random: err4\n");
        goto cleanup;
      }
    }

    cleanup_message(&request);
    cleanup_message(&response);
    request = NULL;
    response = NULL;
    
    samplecount += nsamples*nchans;
    toc = time(NULL);
    elapsed = toc-tic;
    
   fprintf(stderr, "samplecount = %d, elapsed = %d, samples/sec = %f\n", samplecount, elapsed, ((float)(samplecount))/((float)elapsed));
    
  } /* while(1) */
  
  cleanup:
    cleanup_event(&event);
    cleanup_data(&data);
    cleanup_header(&header);
    cleanup_message(&request);
    cleanup_message(&response);
    
    pthread_exit(0);
    return 0;
}

