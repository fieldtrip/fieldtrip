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
#include <signal.h>

#include "buffer.h"
#include "socketserver.h"

#define BUFLEN   1388
#define JAGAPORT 55000
#define FTHOST   "-"
#define FTPORT   1972

static char usage[] =
"\n" \
    "Usage: jaga2ft [ftHost] [ftPort]\n" \
    "\n" \
    "When ftPort is omitted, it will default to 1972.\n" \
    "When ftHost is omitted, it will default to '-'.\n" \
    "Using '-' for the buffer hostname (ftHost) starts a local buffer on the given port (ftPort).\n" \
    "\n" \
    "Example use:\n" \
    "   jaga2ft                 # start a local buffer on port 1972\n" \
    "   jaga2ft - 1234          # start a local buffer on port 1234\n" \
    "   jaga2ft mentat002 1234  # connect to remote buffer\n" \
    "\n" \
    ;  

int keepRunning = 1;

void abortHandler(int sig) {
  printf("Ctrl-C pressed -- stopping...\n");
  keepRunning = 0;
}

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
  } packet_v0;

  struct {
    uint8_t version;
    uint8_t nchans;
    uint16_t diagnostic_word;
    uint16_t mode_word;
    uint16_t fsample;
    uint32_t smp;
  } packet_v3;

  /* this is the common denominator of packet format v0 and v3 */
  struct {
    uint16_t version;
    uint16_t nchans;
    uint16_t nbit;
    uint16_t fsample;
    uint32_t smp;
  } packet;

  int sample = 0, status = 0, verbose = 0;
  host_t host;

  /* these represent the acquisition system properties */
  int nchans         = 16;   /* will be updated later on */
  int fsample        = 1000; /* will be updated later on */
  int nbit           = 16;
  int blocksize      = 43;

  /* these are used in the communication with the FT buffer and represent statefull information */
  int ftSocket           = -1;
  ft_buffer_server_t *ftServer;
  message_t    *request  = NULL;
  message_t    *response = NULL;
  header_t     *header   = NULL;
  data_t       *data     = NULL;

  printf(usage);

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

  /* Spawn tcpserver or connect to remote buffer */
  if (strcmp(host.name, "-") == 0) {
    ftServer = ft_start_buffer_server(host.port, NULL, NULL, NULL);
    if (ftServer==NULL) {
      fprintf(stderr, "jaga2ft: could not start up a local buffer serving at port %i\n", host.port);
      return 1;
    }
    ftSocket = 0;
    printf("jaga2ft: streaming to local buffer on port %i\n", host.port);
  }
  else {
    ftSocket = open_connection(host.name, host.port);

    if (ftSocket < 0) {
      fprintf(stderr, "jaga2ft: could not connect to remote buffer at %s:%i\n", host.name, host.port);
      return 1;
    }
    printf("jaga2ft: streaming to remote buffer at %s:%i\n", host.name, host.port);
  }  

  /* open the UDP server */
  if ((udpsocket=socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP))==-1)
    diep("socket udp");
  int enable = 1;
  if (setsockopt(udpsocket, SOL_SOCKET, SO_REUSEADDR, &enable, sizeof(int)) < 0)
    diep("setsockopt");
  memset((char *) &si_me, 0, sizeof(si_me));
  si_me.sin_family      = AF_INET;
  si_me.sin_port        = htons(JAGAPORT);
  si_me.sin_addr.s_addr = htonl(INADDR_ANY);
  if (bind(udpsocket, &si_me, sizeof(si_me))==-1)
    diep("bind udp");

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
  if (buf[0]==0) {
    packet_v0.version = *(uint16_t *)(buf+0);
    packet_v0.nchans  = *(uint16_t *)(buf+2);
    packet_v0.nbit    = *(uint16_t *)(buf+4);
    packet_v0.fsample = *(uint16_t *)(buf+6);
    packet_v0.sec     = *(uint16_t *)(buf+8);
    packet_v0.smp     = *(uint16_t *)(buf+10);
    /* the packets are quite similar, the data starts at the same location */
    packet.version = packet_v0.version;
    packet.nchans  = packet_v0.nchans;
    packet.nbit    = packet_v0.nbit;
    packet.fsample = packet_v0.fsample;
    packet.smp     = packet_v0.smp;
  }
  else if (buf[0]==3) {
    packet_v3.version         = *(uint8_t  *)(buf+0);
    packet_v3.nchans          = *(uint8_t  *)(buf+1);
    packet_v3.diagnostic_word = *(uint16_t *)(buf+2);
    packet_v3.mode_word       = *(uint16_t *)(buf+4);
    packet_v3.fsample         = *(uint16_t *)(buf+6);
    packet_v3.smp             = *(uint32_t *)(buf+8);
    /* the packets are quite similar, the data starts at the same location */
    packet.version = packet_v3.version;
    packet.nchans  = packet_v3.nchans;
    packet.nbit    = 16;
    packet.fsample = packet_v3.fsample;
    packet.smp     = packet_v3.smp;
  }
  else {
    fprintf(stderr, "invalid packet version");
    exit(1);
  }

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

  status = clientrequest(ftSocket, request, &response);
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

  /* register CTRL-C handler */
  signal(SIGINT, abortHandler);

  printf("Starting to listen - press CTRL-C to quit\n");

  /* add a small pause between writing header + first data block */
  usleep(200000);

  while (keepRunning) 
  {

    if (verbose>1) 
      for (n=0; n<12; n++) 
	printf("buf[%2u] = %hhu\n", n, buf[n]);

    /* parse the UDP package */
    if (buf[0]==0) {
      packet_v0.version = *(uint16_t *)(buf+0);
      packet_v0.nchans  = *(uint16_t *)(buf+2);
      packet_v0.nbit    = *(uint16_t *)(buf+4);
      packet_v0.fsample = *(uint16_t *)(buf+6);
      packet_v0.sec     = *(uint16_t *)(buf+8);
      packet_v0.smp     = *(uint16_t *)(buf+10);
      /* the packets are quite similar, the data starts at the same location */
      packet.version = packet_v0.version;
      packet.nchans  = packet_v0.nchans;
      packet.nbit    = packet_v0.nbit;
      packet.fsample = packet_v0.fsample;
      packet.smp     = packet_v0.smp;
    }
    else if (buf[0]==3) {
      packet_v3.version         = *(uint8_t  *)(buf+0);
      packet_v3.nchans          = *(uint8_t  *)(buf+1);
      packet_v3.diagnostic_word = *(uint16_t *)(buf+2);
      packet_v3.mode_word       = *(uint16_t *)(buf+4);
      packet_v3.fsample         = *(uint16_t *)(buf+6);
      packet_v3.smp             = *(uint32_t *)(buf+8);
      /* the packets are quite similar, the data starts at the same location */
      packet.version = packet_v3.version;
      packet.nchans  = packet_v3.nchans;
      packet.nbit    = 16;
      packet.fsample = packet_v3.fsample;
      packet.smp     = packet_v3.smp;
    }
    else {
      fprintf(stderr, "invalid packet version");
      exit(1);
    }

    /* point to the data */
    data->buf = (buf+12);

    /* do some sanity checks */
    if (packet.nchans!=nchans) {
      fprintf(stderr, "jaga2ft: inconsistent number of channels %hu\n", packet.nchans);
      exit(1);
    }
    if (packet.nbit!=nbit) {
      fprintf(stderr, "jaga2ft: inconsistent number of bits %hu\n", packet.nbit);
      exit(1);
    }
    if (packet.fsample!=fsample) {
      fprintf(stderr, "jaga2ft: inconsistent sampling rate %hu\n", packet.fsample);
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

    status = clientrequest(ftSocket, request, &response);
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

  data->buf = NULL; /* this is not allocated on the heap but pointing to somewhere on the stack */
  cleanup_data(&data);
  close(udpsocket);

  if (ftSocket > 0) {
    close_connection(ftSocket);
  } else {
    ft_stop_buffer_server(ftServer);
  }
  return 0;
}
