/*
 * Copyright (C) 2015-2017, Robert Oostenveld
 *
 * Donders Institute for Brain, Cognition and Behaviour,
 * Centre for Cognitive Neuroimaging, Radboud University,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 */

#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <signal.h>

#include "platform.h"
#include "platform_includes.h"
#include "buffer.h"
#include "socketserver.h"
#include "ini.h"
#include "timestamp.h"

#define JAGA_PORT      55000
#define JAGA_NCHAN     16
#define JAGA_BLOCKSIZE 43
#define JAGA_BUFLEN    1388
#define JAGA_NBIT      16
#define JAGA_FSAMPLE   1000

static char usage[] =
"\n" \
"Use as\n" \
"   jaga2ft [ftHost] [ftPort]\n" \
"with the parameters as specified below, or\n" \
"   jaga2ft <config>\n" \
"with the name of a configuration file for detailed setup.\n"
"\n" \
"When ftPort is omitted, it will default to 1972.\n" \
"When ftHost is omitted, it will default to '-'.\n" \
"Using '-' for the buffer hostname (ftHost) starts a local buffer on the given port (ftPort).\n" \
"\n" \
"Example use:\n" \
"   jaga2ft                 # start a local buffer on port 1972\n" \
"   jaga2ft - 1234          # start a local buffer on port 1234\n" \
"   jaga2ft serverpc 1234   # connect to remote buffer running on server PC\n" \
"\n" \
;

typedef struct {
  int        port;
  int        verbose;
  const char *hostname;
  const char *timestamp;
  const char *timeref;

  const char *enable_chan1;
  const char *enable_chan2;
  const char *enable_chan3;
  const char *enable_chan4;
  const char *enable_chan5;
  const char *enable_chan6;
  const char *enable_chan7;
  const char *enable_chan8;
  const char *enable_chan9;
  const char *enable_chan10;
  const char *enable_chan11;
  const char *enable_chan12;
  const char *enable_chan13;
  const char *enable_chan14;
  const char *enable_chan15;
  const char *enable_chan16;

  const char *label_chan1;
  const char *label_chan2;
  const char *label_chan3;
  const char *label_chan4;
  const char *label_chan5;
  const char *label_chan6;
  const char *label_chan7;
  const char *label_chan8;
  const char *label_chan9;
  const char *label_chan10;
  const char *label_chan11;
  const char *label_chan12;
  const char *label_chan13;
  const char *label_chan14;
  const char *label_chan15;
  const char *label_chan16;
  const char *label_chan17;  /* this is for the timestamp */
} configuration;

int keepRunning = 1;

static int iniHandler(void* external, const char* section, const char* name, const char* value) {
  configuration *local = (configuration*)external;

#define MATCH(s, n) strcasecmp(section, s) == 0 && strcasecmp(name, n) == 0
  if (MATCH("General", "port")) {
    local->port = atoi(value);
  } else if (MATCH("General", "verbose")) {
    local->verbose = atoi(value);
  } else if (MATCH("General", "hostname")) {
    local->hostname = strdup(value);
  } else if (MATCH("General", "timestamp")) {
    local->timestamp = strdup(value);
  } else if (MATCH("General", "timeref")) {
    local->timeref = strdup(value);

  } else if (MATCH("ChannelEnable", "chan1")) {
    local->enable_chan1 = strdup(value);
  } else if (MATCH("ChannelEnable", "chan2")) {
    local->enable_chan2 = strdup(value);
  } else if (MATCH("ChannelEnable", "chan3")) {
    local->enable_chan3 = strdup(value);
  } else if (MATCH("ChannelEnable", "chan4")) {
    local->enable_chan4 = strdup(value);
  } else if (MATCH("ChannelEnable", "chan5")) {
    local->enable_chan5 = strdup(value);
  } else if (MATCH("ChannelEnable", "chan6")) {
    local->enable_chan6 = strdup(value);
  } else if (MATCH("ChannelEnable", "chan7")) {
    local->enable_chan7 = strdup(value);
  } else if (MATCH("ChannelEnable", "chan8")) {
    local->enable_chan8 = strdup(value);
  } else if (MATCH("ChannelEnable", "chan9")) {
    local->enable_chan9 = strdup(value);
  } else if (MATCH("ChannelEnable", "chan10")) {
    local->enable_chan10 = strdup(value);
  } else if (MATCH("ChannelEnable", "chan11")) {
    local->enable_chan11 = strdup(value);
  } else if (MATCH("ChannelEnable", "chan12")) {
    local->enable_chan12 = strdup(value);
  } else if (MATCH("ChannelEnable", "chan13")) {
    local->enable_chan13 = strdup(value);
  } else if (MATCH("ChannelEnable", "chan14")) {
    local->enable_chan14 = strdup(value);
  } else if (MATCH("ChannelEnable", "chan15")) {
    local->enable_chan15 = strdup(value);
  } else if (MATCH("ChannelEnable", "chan16")) {
    local->enable_chan16 = strdup(value);

  } else if (MATCH("ChannelLabel", "chan1")) {
    local->label_chan1 = strdup(value);
  } else if (MATCH("ChannelLabel", "chan2")) {
    local->label_chan2 = strdup(value);
  } else if (MATCH("ChannelLabel", "chan3")) {
    local->label_chan3 = strdup(value);
  } else if (MATCH("ChannelLabel", "chan4")) {
    local->label_chan4 = strdup(value);
  } else if (MATCH("ChannelLabel", "chan5")) {
    local->label_chan5 = strdup(value);
  } else if (MATCH("ChannelLabel", "chan6")) {
    local->label_chan6 = strdup(value);
  } else if (MATCH("ChannelLabel", "chan7")) {
    local->label_chan7 = strdup(value);
  } else if (MATCH("ChannelLabel", "chan8")) {
    local->label_chan8 = strdup(value);
  } else if (MATCH("ChannelLabel", "chan9")) {
    local->label_chan9 = strdup(value);
  } else if (MATCH("ChannelLabel", "chan10")) {
    local->label_chan10 = strdup(value);
  } else if (MATCH("ChannelLabel", "chan11")) {
    local->label_chan11 = strdup(value);
  } else if (MATCH("ChannelLabel", "chan12")) {
    local->label_chan12 = strdup(value);
  } else if (MATCH("ChannelLabel", "chan13")) {
    local->label_chan13 = strdup(value);
  } else if (MATCH("ChannelLabel", "chan14")) {
    local->label_chan14 = strdup(value);
  } else if (MATCH("ChannelLabel", "chan15")) {
    local->label_chan15 = strdup(value);
  } else if (MATCH("ChannelLabel", "chan16")) {
    local->label_chan16 = strdup(value);
  } else if (MATCH("ChannelLabel", "chan17")) {
    local->label_chan17 = strdup(value);

  } else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}

void abortHandler(int sig) {
  printf("Ctrl-C pressed -- stopping...\n");
  keepRunning = 0;
}

void diep(char *s) {
  perror(s);
  exit(1);
}

int file_exists(const char *fname) {
  /* see http://stackoverflow.com/questions/230062/whats-the-best-way-to-check-if-a-file-exists-in-c-cross-platform */
  if( access( fname, F_OK ) != -1 )
    return 1;
  else
    return 0;
}

int main(int argc, char *argv[]) {
  struct sockaddr_in si_me, si_other;
  socklen_t slen = sizeof(struct sockaddr_in);
  int udpsocket, n, count = 0, sample = 0, chan = 0, status = 0, labelSize;
  char buf[JAGA_BUFLEN], *labelString;
  host_t host;
  struct timespec tic, toc;

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

  /* these represent the acquisition system properties */
  int nchans         = JAGA_NCHAN;     /* will be updated later on */
  int fsample        = JAGA_FSAMPLE;   /* will be updated later on */
  int blocksize      = JAGA_BLOCKSIZE;

  /* these are used in the communication with the FT buffer and represent statefull information */
  int ftSocket            = -1;
  ft_buffer_server_t *ftServer;
  message_t     *request  = NULL;
  message_t     *response = NULL;
  header_t      *header   = NULL;
  data_t        *data     = NULL;
  ft_chunkdef_t *label    = NULL;

  /* this structure contains the configuration details */
  configuration config;

  /* configure the default settings */
  config.port          = 1972;
  config.verbose       = 1;
  config.hostname      = strdup("-");
  config.timestamp     = strdup("on");
  config.timeref       = strdup("start");

  config.enable_chan1  = strdup("on");
  config.enable_chan2  = strdup("on");
  config.enable_chan3  = strdup("on");
  config.enable_chan4  = strdup("on");
  config.enable_chan5  = strdup("on");
  config.enable_chan6  = strdup("on");
  config.enable_chan7  = strdup("on");
  config.enable_chan8  = strdup("on");
  config.enable_chan9  = strdup("on");
  config.enable_chan10 = strdup("on");
  config.enable_chan11 = strdup("on");
  config.enable_chan12 = strdup("on");
  config.enable_chan13 = strdup("on");
  config.enable_chan14 = strdup("on");
  config.enable_chan15 = strdup("on");
  config.enable_chan16 = strdup("on");

  config.label_chan1  = strdup("1");
  config.label_chan2  = strdup("2");
  config.label_chan3  = strdup("3");
  config.label_chan4  = strdup("4");
  config.label_chan5  = strdup("5");
  config.label_chan6  = strdup("6");
  config.label_chan7  = strdup("7");
  config.label_chan8  = strdup("8");
  config.label_chan9  = strdup("9");
  config.label_chan10 = strdup("10");
  config.label_chan11 = strdup("11");
  config.label_chan12 = strdup("12");
  config.label_chan13 = strdup("13");
  config.label_chan14 = strdup("14");
  config.label_chan15 = strdup("15");
  config.label_chan16 = strdup("16");
  config.label_chan17 = strdup("TimeStamp");

  /* parse the command line arguments */
  if (argc==1) {
    /* use the defaults */
    strcpy(host.name, config.hostname);
    host.port = config.port;
  }
  else if (argc==2 && (strcasecmp(argv[1], "-h")==0 || strcasecmp(argv[1], "-?")==0)) {
    /* show the help */
    printf("%s", usage);
    exit(0);
  }
  else if (argc==2 && (file_exists(argv[1]))) {
    /* the second argument is the configuration file */
    fprintf(stderr, "jaga2ft: loading configuration from '%s'\n", argv[1]);
    if (ini_parse(argv[1], iniHandler, &config) < 0) {
      fprintf(stderr, "Can't load '%s'\n", argv[1]);
      return 1;
    }
    /* get the hostname and port from the configuration file */
    strcpy(host.name, config.hostname);
    host.port = config.port;
  }
  else if (argc==2) {
    /* the second argument is the remote server */
    strcpy(host.name, argv[1]);
  }
  else if (argc==3) {
    /* the second and third argument are the remote server and port */
    strcpy(host.name, argv[1]);
    host.port = atoi(argv[2]);
  }
  else {
    printf("%s", usage);
    exit(0);
  }

  /* count the number of channels that should be sent */
#define ISTRUE(s)  strcasecmp(s, "on")==0
  nchans = 0;
  if (ISTRUE(config.enable_chan1))
    nchans++;
  if (ISTRUE(config.enable_chan2))
    nchans++;
  if (ISTRUE(config.enable_chan3))
    nchans++;
  if (ISTRUE(config.enable_chan4))
    nchans++;
  if (ISTRUE(config.enable_chan5))
    nchans++;
  if (ISTRUE(config.enable_chan6))
    nchans++;
  if (ISTRUE(config.enable_chan7))
    nchans++;
  if (ISTRUE(config.enable_chan8))
    nchans++;
  if (ISTRUE(config.enable_chan9))
    nchans++;
  if (ISTRUE(config.enable_chan10))
    nchans++;
  if (ISTRUE(config.enable_chan11))
    nchans++;
  if (ISTRUE(config.enable_chan12))
    nchans++;
  if (ISTRUE(config.enable_chan13))
    nchans++;
  if (ISTRUE(config.enable_chan14))
    nchans++;
  if (ISTRUE(config.enable_chan15))
    nchans++;
  if (ISTRUE(config.enable_chan16))
    nchans++;
  if (ISTRUE(config.timestamp))
    nchans++;

  fprintf(stderr, "jaga2ft: hostname     =  %s\n", host.name);
  fprintf(stderr, "jaga2ft: port         =  %d\n", host.port);
  fprintf(stderr, "jaga2ft: timestamp    =  %s\n", config.timestamp);
  fprintf(stderr, "jaga2ft: nchans       =  %d\n", nchans);

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
  si_me.sin_port        = htons(JAGA_PORT);
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
  if ((n=recvfrom(udpsocket, buf, JAGA_BUFLEN, 0, &si_other, &slen))==-1)
    diep("recvfrom()");
  if (config.verbose>0)
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
    packet.nbit    = JAGA_NBIT;
    packet.fsample = packet_v3.fsample;
    packet.smp     = packet_v3.smp;
  }
  else {
    fprintf(stderr, "invalid packet version");
    exit(1);
  }

  if (config.verbose>0) {
    printf("jaga2ft: packet.version = %d\n", packet.version);
    printf("jaga2ft: packet.nchans  = %d\n", packet.nchans);
    printf("jaga2ft: packet.nbit    = %d\n", packet.nbit);
    printf("jaga2ft: packet.fsample = %d\n", packet.fsample);
    printf("jaga2ft: packet.smp     = %d\n", packet.smp);
  }

  /* update the defaults */
  fsample = packet.fsample;

  /* define the header */
  header->def->nchans    = nchans;
  header->def->nsamples  = 0;
  header->def->nevents   = 0;
  header->def->fsample   = fsample;
  header->def->data_type = DATATYPE_FLOAT32;
  header->def->bufsize   = 0;

  labelSize = 0; /* count the number of bytes required */
  if (ISTRUE (config.enable_chan1))
    labelSize += strlen (config.label_chan1) + 1;
  if (ISTRUE (config.enable_chan2))
    labelSize += strlen (config.label_chan2) + 1;
  if (ISTRUE (config.enable_chan3))
    labelSize += strlen (config.label_chan3) + 1;
  if (ISTRUE (config.enable_chan4))
    labelSize += strlen (config.label_chan4) + 1;
  if (ISTRUE (config.enable_chan5))
    labelSize += strlen (config.label_chan5) + 1;
  if (ISTRUE (config.enable_chan6))
    labelSize += strlen (config.label_chan6) + 1;
  if (ISTRUE (config.enable_chan7))
    labelSize += strlen (config.label_chan7) + 1;
  if (ISTRUE (config.enable_chan8))
    labelSize += strlen (config.label_chan8) + 1;
  if (ISTRUE (config.enable_chan9))
    labelSize += strlen (config.label_chan9) + 1;
  if (ISTRUE (config.enable_chan10))
    labelSize += strlen (config.label_chan10) + 1;
  if (ISTRUE (config.enable_chan11))
    labelSize += strlen (config.label_chan11) + 1;
  if (ISTRUE (config.enable_chan12))
    labelSize += strlen (config.label_chan12) + 1;
  if (ISTRUE (config.enable_chan13))
    labelSize += strlen (config.label_chan13) + 1;
  if (ISTRUE (config.enable_chan14))
    labelSize += strlen (config.label_chan14) + 1;
  if (ISTRUE (config.enable_chan15))
    labelSize += strlen (config.label_chan15) + 1;
  if (ISTRUE (config.enable_chan16))
    labelSize += strlen (config.label_chan16) + 1;
  if (ISTRUE (config.timestamp))
    labelSize += strlen (config.label_chan17) + 1;

  /* go over all channels for a 2nd time, now copying the strings to the destination */
  labelString = (char *) malloc (labelSize * sizeof(char));
  labelSize   = 0;
  if (ISTRUE (config.enable_chan1)) {
    strcpy (labelString+labelSize, config.label_chan1);
    labelSize += strlen (config.label_chan1) + 1;
  }
  if (ISTRUE (config.enable_chan2)) {
    strcpy (labelString+labelSize, config.label_chan2);
    labelSize += strlen (config.label_chan2) + 1;
  }
  if (ISTRUE (config.enable_chan3)) {
    strcpy (labelString+labelSize, config.label_chan3);
    labelSize += strlen (config.label_chan3) + 1;
  }
  if (ISTRUE (config.enable_chan4)) {
    strcpy (labelString+labelSize, config.label_chan4);
    labelSize += strlen (config.label_chan4) + 1;
  }
  if (ISTRUE (config.enable_chan5)) {
    strcpy (labelString+labelSize, config.label_chan5);
    labelSize += strlen (config.label_chan5) + 1;
  }
  if (ISTRUE (config.enable_chan6)) {
    strcpy (labelString+labelSize, config.label_chan6);
    labelSize += strlen (config.label_chan6) + 1;
  }
  if (ISTRUE (config.enable_chan7)) {
    strcpy (labelString+labelSize, config.label_chan7);
    labelSize += strlen (config.label_chan7) + 1;
  }
  if (ISTRUE (config.enable_chan8)) {
    strcpy (labelString+labelSize, config.label_chan8);
    labelSize += strlen (config.label_chan8) + 1;
  }
  if (ISTRUE (config.enable_chan9)) {
    strcpy (labelString+labelSize, config.label_chan9);
    labelSize += strlen (config.label_chan9) + 1;
  }
  if (ISTRUE (config.enable_chan10)) {
    strcpy (labelString+labelSize, config.label_chan10);
    labelSize += strlen (config.label_chan10) + 1;
  }
  if (ISTRUE (config.enable_chan11)) {
    strcpy (labelString+labelSize, config.label_chan11);
    labelSize += strlen (config.label_chan11) + 1;
  }
  if (ISTRUE (config.enable_chan12)) {
    strcpy (labelString+labelSize, config.label_chan12);
    labelSize += strlen (config.label_chan12) + 1;
  }
  if (ISTRUE (config.enable_chan13)) {
    strcpy (labelString+labelSize, config.label_chan13);
    labelSize += strlen (config.label_chan13) + 1;
  }
  if (ISTRUE (config.enable_chan14)) {
    strcpy (labelString+labelSize, config.label_chan14);
    labelSize += strlen (config.label_chan14) + 1;
  }
  if (ISTRUE (config.enable_chan15)) {
    strcpy (labelString+labelSize, config.label_chan15);
    labelSize += strlen (config.label_chan15) + 1;
  }
  if (ISTRUE (config.enable_chan16)) {
    strcpy (labelString+labelSize, config.label_chan16);
    labelSize += strlen (config.label_chan16) + 1;
  }
  if (ISTRUE (config.timestamp)) {
    strcpy (labelString+labelSize, config.label_chan17);
    labelSize += strlen (config.label_chan17) + 1;
  }

  /* add the channel label chunk to the header */
  label = (ft_chunkdef_t *) malloc (sizeof (ft_chunkdef_t));
  label->type = FT_CHUNK_CHANNEL_NAMES;
  label->size = labelSize;
  header->def->bufsize = append (&header->buf, header->def->bufsize, label, sizeof (ft_chunkdef_t));
  header->def->bufsize = append (&header->buf, header->def->bufsize, labelString, labelSize);
  FREE (label);
  FREE (labelString);

  /* define the constant part of the data and allocate space for the variable part */
  data->def->nchans    = nchans;
  data->def->nsamples  = blocksize;
  data->def->data_type = DATATYPE_FLOAT32;
  data->def->bufsize   = WORDSIZE_FLOAT32*nchans*blocksize;
  data->buf = malloc (data->def->bufsize);

  /* initialization phase, send the header */
  request->def->command = PUT_HDR;
  request->def->bufsize = append(&request->buf, request->def->bufsize, header->def, sizeof(headerdef_t));
  request->def->bufsize = append(&request->buf, request->def->bufsize, header->buf, header->def->bufsize);

  /* this is not needed any more */
  cleanup_header(&header);

  status = clientrequest(ftSocket, request, &response);
  if (config.verbose>1) fprintf(stderr, "jaga2ft: clientrequest returned %d\n", status);
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

  /* determine the reference time for the timestamps */
  if (strcasecmp (config.timeref, "start") == 0) {
      /* since the start of the acquisition */
      get_monotonic_time (&tic, TIMESTAMP_REF_BOOT);
  }
  else if (strcasecmp (config.timeref, "boot") == 0) {
      /* since the start of the day */
      tic.tv_sec = 0;
      tic.tv_nsec = 0;
  }
  else if (strcasecmp (config.timeref, "epoch") == 0) {
      /* since the start of the epoch, i.e. 1-1-1970 */
      tic.tv_sec = 0;
      tic.tv_nsec = 0;
  }
  else {
      fprintf (stderr, "Incorrect specification of timeref, should be 'start', 'day' or 'epoch'\n");
      return 1;
  }

  while (keepRunning) {

    if (config.verbose>2)
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
      packet.nbit    = JAGA_NBIT;
      packet.fsample = packet_v3.fsample;
      packet.smp     = packet_v3.smp;
    }
    else {
      fprintf(stderr, "invalid packet version");
      exit(1);
    }

    /* do some sanity checks */
    if (packet.nchans!=JAGA_NCHAN) {
      fprintf(stderr, "jaga2ft: inconsistent number of channels %hu\n", packet.nchans);
      exit(1);
    }
    if (packet.nbit!=JAGA_NBIT) {
      fprintf(stderr, "jaga2ft: inconsistent number of bits %hu\n", packet.nbit);
      exit(1);
    }
    if (packet.fsample!=fsample) {
      fprintf(stderr, "jaga2ft: inconsistent sampling rate %hu\n", packet.fsample);
      exit(1);
    }

    /* loop over all channels and samples in the block, copy the values from enabled channels */
    sample = 0;
    while (sample<blocksize) {

      UINT16_T *raw = (buf+12);
      chan = 0;
      if (ISTRUE (config.enable_chan1))
	((FLOAT32_T *) (data->buf))[nchans * sample + chan++] = raw[nchans * sample + 0];
      if (ISTRUE (config.enable_chan2))
	((FLOAT32_T *) (data->buf))[nchans * sample + chan++] = raw[nchans * sample + 1];
      if (ISTRUE (config.enable_chan3))
	((FLOAT32_T *) (data->buf))[nchans * sample + chan++] = raw[nchans * sample + 2];
      if (ISTRUE (config.enable_chan4))
	((FLOAT32_T *) (data->buf))[nchans * sample + chan++] = raw[nchans * sample + 3];
      if (ISTRUE (config.enable_chan5))
	((FLOAT32_T *) (data->buf))[nchans * sample + chan++] = raw[nchans * sample + 4];
      if (ISTRUE (config.enable_chan6))
	((FLOAT32_T *) (data->buf))[nchans * sample + chan++] = raw[nchans * sample + 5];
      if (ISTRUE (config.enable_chan7))
	((FLOAT32_T *) (data->buf))[nchans * sample + chan++] = raw[nchans * sample + 6];
      if (ISTRUE (config.enable_chan8))
	((FLOAT32_T *) (data->buf))[nchans * sample + chan++] = raw[nchans * sample + 7];
      if (ISTRUE (config.enable_chan9))
	((FLOAT32_T *) (data->buf))[nchans * sample + chan++] = raw[nchans * sample + 8];
      if (ISTRUE (config.enable_chan10))
	((FLOAT32_T *) (data->buf))[nchans * sample + chan++] = raw[nchans * sample + 9];
      if (ISTRUE (config.enable_chan11))
	((FLOAT32_T *) (data->buf))[nchans * sample + chan++] = raw[nchans * sample + 10];
      if (ISTRUE (config.enable_chan12))
	((FLOAT32_T *) (data->buf))[nchans * sample + chan++] = raw[nchans * sample + 11];
      if (ISTRUE (config.enable_chan13))
	((FLOAT32_T *) (data->buf))[nchans * sample + chan++] = raw[nchans * sample + 12];
      if (ISTRUE (config.enable_chan14))
	((FLOAT32_T *) (data->buf))[nchans * sample + chan++] = raw[nchans * sample + 13];
      if (ISTRUE (config.enable_chan15))
	((FLOAT32_T *) (data->buf))[nchans * sample + chan++] = raw[nchans * sample + 14];
      if (ISTRUE (config.enable_chan16))
	((FLOAT32_T *) (data->buf))[nchans * sample + chan++] = raw[nchans * sample + 15];

      if (ISTRUE (config.timestamp)) {
	if (strcasecmp (config.timeref, "start") == 0)
	  get_monotonic_time (&toc, TIMESTAMP_REF_BOOT);
	else if (strcasecmp (config.timeref, "boot") == 0)
	  get_monotonic_time (&toc, TIMESTAMP_REF_BOOT);
	else if (strcasecmp (config.timeref, "epoch") == 0)
	  get_monotonic_time (&toc, TIMESTAMP_REF_EPOCH);
	((FLOAT32_T *) (data->buf))[nchans * sample + chan++] = get_elapsed_time (&tic, &toc);
      }

      sample++;
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
    if (config.verbose>1) fprintf(stderr, "jaga2ft: clientrequest returned %d\n", status);
    if (status) {
      fprintf(stderr, "jaga2ft: err3\n");
      exit(1);
    }

    if (status) {
      fprintf(stderr, "jaga2ft: err4\n");
      exit(1);
    }

    count += sample;
    printf("jaga2ft: sample count = %i\n", count);

    /* FIXME do someting with the response, i.e. check that it is OK */
    cleanup_message(&request);

    if (response == NULL || response->def == NULL || response->def->command!=PUT_OK) {
      fprintf(stderr, "Error when writing samples.\n");
    }
    cleanup_message(&response);

    /* read the next packet */
    if ((n=recvfrom(udpsocket, buf, JAGA_BUFLEN, 0, &si_other, &slen))==-1)
      diep("recvfrom()");
    if (config.verbose>1)
      printf("jaga2ft: received %d byte packet from %s:%d\n", n, inet_ntoa(si_other.sin_addr), ntohs(si_other.sin_port));

  }				/* while keepRunning */

  cleanup_data(&data);
  close(udpsocket);

  if (ftSocket > 0)
    close_connection(ftSocket);
  else
    ft_stop_buffer_server(ftServer);

  return 0;

}				/* main */
