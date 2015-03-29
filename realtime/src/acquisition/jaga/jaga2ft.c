/*
 * Use as
 *   jaga2ft
 *   jaga2ft <hostname>
 *   jaga2ft <hostname> <port>
 *
 * Copyright (C) 2015 Robert Oostenveld
 * Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
 * Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 */

#include <arpa/inet.h>
#include <netinet/in.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "buffer.h"

#define BUFLEN   1388
#define JAGAPORT 55000
#define FTHOST   "localhost"
#define FTPORT   1972

void diep(char *s)
{
  perror(s);
  exit(1);
}

int main(int argc, char *argv[])
{
  struct sockaddr_in si_me, si_other;
  socklen_t slen = sizeof(struct sockaddr_in);
  int udpsocket, n;
  char buf[BUFLEN];

  struct {
    uint16_t version;
    uint16_t nchans;
    uint16_t nbit;
    uint16_t fsample;
    uint16_t sec;
    uint16_t smp;
  } packet;

  int sample = 0, status = 0, verbose = 0;
  host_t host;

  /* these represent the acquisition system properties */
  int nchans         = 16;   /* will be updated later on */
  int fsample        = 1000; /* will be updated later on */
  int nbit           = 16;
  int blocksize      = 43;

  /* these are used in the communication with the FT buffer and represent statefull information */
  int ftserver           = -1;
  message_t    *request  = NULL;
  message_t    *response = NULL;
  header_t     *header   = NULL;
  data_t       *data     = NULL;

  if (argc<3) {
    printf("Usage: %s <hostname> <port>\n", argv[0]);
  }

  if (argc>1)
    strcpy(host.name, argv[1]);
  else {
    strcpy(host.name, FTHOST);
  }

  if (argc>2)
    host.port = atoi(argv[2]);
  else {
    host.port = FTPORT;
  }

  fprintf(stderr, "jaga2ft: hostname     =  %s\n", host.name);
  fprintf(stderr, "jaga2ft: port         =  %d\n", host.port);
  ftserver = open_connection(host.name, host.port);
  if (ftserver<0) {
    fprintf(stderr, "jaga2ft: could not open connection to buffer\n");
    exit(1);
  }

  /* open the UDP server */
  if ((udpsocket=socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP))==-1)
    diep("socket");
  memset((char *) &si_me, 0, sizeof(si_me));
  si_me.sin_family      = AF_INET;
  si_me.sin_port        = htons(JAGAPORT);
  si_me.sin_addr.s_addr = htonl(INADDR_ANY);
  if (bind(udpsocket, &si_me, sizeof(si_me))==-1)
    diep("bind");

  /* allocate the elements that will be used in the communication to the FT buffer */
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

  /* read the first packet to get some information */
  if ((n=recvfrom(udpsocket, buf, BUFLEN, 0, &si_other, &slen))==-1)
    diep("recvfrom()");
  if (verbose>0)
    printf("jaga2ft: received %d byte packet from %s:%d\n", n, inet_ntoa(si_other.sin_addr), ntohs(si_other.sin_port));

  /* parse the UDP package */
  packet.version = *(uint16_t *)(buf+0);
  packet.nchans  = *(uint16_t *)(buf+2);
  packet.nbit    = *(uint16_t *)(buf+4);
  packet.fsample = *(uint16_t *)(buf+6);
  packet.sec     = *(uint16_t *)(buf+8);
  packet.smp     = *(uint16_t *)(buf+10);

  /* update the defaults */
  nchans  = packet.nchans;
  fsample = packet.fsample;

  /* define the header */
  header->def->nchans    = nchans;
  header->def->nsamples  = 0;
  header->def->nevents   = 0;
  header->def->fsample   = fsample;
  header->def->data_type = DATATYPE_UINT16;
  header->def->bufsize   = 0;

  /* define the constant part of the data and allocate space for the variable part */
  data->def->nchans    = nchans;
  data->def->nsamples  = blocksize;
  data->def->data_type = DATATYPE_UINT16;
  data->def->bufsize   = WORDSIZE_UINT16*nchans*blocksize;

  /* initialization phase, send the header */
  request->def->command = PUT_HDR;
  request->def->bufsize = append(&request->buf, request->def->bufsize, header->def, sizeof(headerdef_t));
  request->def->bufsize = append(&request->buf, request->def->bufsize, header->buf, header->def->bufsize);

  /* this is not needed any more */
  cleanup_header(&header);

  status = clientrequest(ftserver, request, &response);
  if (verbose>0) fprintf(stderr, "jaga2ft: clientrequest returned %d\n", status);
  if (status) {
    fprintf(stderr, "jaga2ft: could not send request to buffer\n");
    exit(1);
  }

  if (status || response==NULL || response->def == NULL) {
    fprintf(stderr, "jaga2ft: err2\n");
    exit(1);
  }

  cleanup_message(&request);

  if (response->def->command != PUT_OK) {
    fprintf(stderr, "jaga2ft: error in 'put header' request.\n");
    exit(1);
  }

  cleanup_message(&response);

  /* add a small pause between writing header + first data block */
  usleep(200000);

  while (1) 
  {

    if (verbose>1) 
      for (unsigned int i=0; i<12; i++) 
        printf("buf[%2u] = %hhu\n", i, buf[i]);

    /* parse the UDP package */
    packet.version = *(uint16_t *)(buf+0);
    packet.nchans  = *(uint16_t *)(buf+2);
    packet.nbit    = *(uint16_t *)(buf+4);
    packet.fsample = *(uint16_t *)(buf+6);
    packet.sec     = *(uint16_t *)(buf+8);
    packet.smp     = *(uint16_t *)(buf+10);

    /* point to the data */
    data->buf = (buf+12);

    /* do some sanity checks */
    if (packet.version!=0) {
      fprintf(stderr, "jaga2ft: inconsistent version %hu\n", packet.version);
      exit(1);
    }
    if (packet.nchans!=nchans) {
      fprintf(stderr, "jaga2ft: inconsistent number of channels %hu\n", packet.nchans);
      exit(1);
    }
    if (packet.nbit!=nbit) {
      fprintf(stderr, "jaga2ft: inconsistent number of bits %hu\n", packet.nbit);
      exit(1);
    }
    if (packet.fsample!=fsample) {
      fprintf(stderr, "jaga2ft: inconsistent samling rate %hu\n", packet.fsample);
      exit(1);
    }

    /* create the request */
    request      = malloc(sizeof(message_t));
    request->def = malloc(sizeof(messagedef_t));
    request->buf = NULL;
    request->def->version = VERSION;
    request->def->bufsize = 0;
    request->def->command = PUT_DAT;
    request->def->bufsize = append(&request->buf, request->def->bufsize, data->def, sizeof(datadef_t));
    request->def->bufsize = append(&request->buf, request->def->bufsize, data->buf, data->def->bufsize);

    status = clientrequest(ftserver, request, &response);
    if (verbose>0) fprintf(stderr, "jaga2ft: clientrequest returned %d\n", status);
    if (status) {
      fprintf(stderr, "jaga2ft: err3\n");
      exit(1);
    }

    if (status) {
      fprintf(stderr, "jaga2ft: err4\n");
      exit(1);
    }

    sample += blocksize;
    printf("jaga2ft: sample count = %i\n", sample);

    /* FIXME do someting with the response, i.e. check that it is OK */
    cleanup_message(&request);

    if (response == NULL || response->def == NULL || response->def->command!=PUT_OK) {
      fprintf(stderr, "Error when writing samples.\n");
    }
    cleanup_message(&response);

    /* read the next packet */
    if ((n=recvfrom(udpsocket, buf, BUFLEN, 0, &si_other, &slen))==-1)
      diep("recvfrom()");
    if (verbose>0)
      printf("jaga2ft: received %d byte packet from %s:%d\n", n, inet_ntoa(si_other.sin_addr), ntohs(si_other.sin_port));

  } /* while(1) */

  close_connection(ftserver);
  close(udpsocket);
  return 0;
}
