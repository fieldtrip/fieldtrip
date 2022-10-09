/*
 * Copyright (C) 2022, Robert Oostenveld
 *
 * Donders Institute for Brain, Cognition and Behaviour,
 * Centre for Cognitive Neuroimaging, Radboud University,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 */

#include <stdio.h>
#include <sys/types.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <signal.h>

#include "platform.h"
#include "serial.h"
#include "buffer.h"
#include "socketserver.h"
#include "timestamp.h"

int keepRunning = 1;

/* these sequences comprise the command byte and two CRC bytes */
unsigned char StartAcquisition[] = {0x61, 0x7C, 0x87};
unsigned char StopAcquisition[] = {0x63, 0x5C, 0xC5};

#if defined (PLATFORM_WIN32)
static char usage[] =
"\n" \
    "Use as\n" \
    "   unicorn2ft <device> [ftHost] [ftPort]\n" \
    "\n" \
    "When ftPort is omitted, it will default to 1972.\n" \
    "When ftHost is omitted, it will default to '-'.\n" \
    "Using '-' for the buffer hostname (ftHost) starts a local buffer on the given port (ftPort).\n" \
    "\n" \
    "Example use:\n" \
    "   unicorn2ft COM3:                 # start a local buffer on port 1972\n" \
    "   unicorn2ft COM3: - 1234          # start a local buffer on port 1234\n" \
    "   unicorn2ft COM3: serverpc 1234   # connect to remote buffer running on server PC\n" \
    "\n" \
    ;
#else
static char usage[] =
"\n" \
    "Use as\n" \
    "   unicorn2ft <device> [ftHost] [ftPort]\n" \
    "\n" \
    "When ftPort is omitted, it will default to 1972.\n" \
    "When ftHost is omitted, it will default to '-'.\n" \
    "Using '-' for the buffer hostname (ftHost) starts a local buffer on the given port (ftPort).\n" \
    "\n" \
    "Example use:\n" \
    "   unicorn2ft /dev/tty.UN-20220110                 # start a local buffer on port 1972\n" \
    "   unicorn2ft /dev/tty.UN-20220110 - 1234          # start a local buffer on port 1234\n" \
    "   unicorn2ft /dev/tty.UN-20220110 serverpc 1234   # connect to remote buffer running on server PC\n" \
    "\n" \
    ;
#endif

int serialWriteSlow(SerialPort *SP, int size, void *buffer) {
  int i, retval = 0;
  for (i=0; i<size; i++) {
    retval += serialWrite(SP, 1, buffer+i);
    usleep(100000);
  }
  return retval;
}

void abortHandler(int sig) {
  printf("Ctrl-C pressed -- stopping...\n");
  keepRunning = 0;
  usleep(1000000);
  exit(0);
}

int main(int argc, char *argv[]) {
  int n, c, status = 0;

  char *SerialDevice = NULL;
  SerialPort SP;
  host_t host;

  /* these represent the general acquisition system properties */
  int nchans         = 17;
  float fsample      = 250;

  /* these are used in the communication with the FT buffer and represent statefull information */
  int ftSocket           = -1;
  ft_buffer_server_t *ftServer;
  message_t     *request  = NULL;
  message_t     *response = NULL;
  header_t      *header   = NULL;
  data_t        *data     = NULL;

  if (argc<2) {
    printf("%s", usage);
    exit(0);
  }

  if (argc==2) {
    if (strncmp(argv[1], "/dev", 4)==0 || strncasecmp(argv[1], "COM", 3)==0)
      /* the second argument is the serial port */
      SerialDevice = strdup(argv[1]);
  }
  else  {
    fprintf(stderr, "unicorn2ft: serial device not specified\n");
    return 1;
  }

  if (argc>2)
    strcpy(host.name, argv[2]);
  else {
    strcpy(host.name, "-");
  }

  if (argc>3)
    host.port = atoi(argv[3]);
  else {
    host.port = 1972;
  }

  fprintf(stderr, "unicorn2ft: serial       =  %s\n", SerialDevice);
  fprintf(stderr, "unicorn2ft: hostname     =  %s\n", host.name);
  fprintf(stderr, "unicorn2ft: port         =  %d\n", host.port);
  fprintf(stderr, "unicorn2ft: nchans       =  %d\n", nchans);
  fprintf(stderr, "unicorn2ft: fsample      =  %f\n", fsample);

  /* Spawn tcpserver or connect to remote buffer */
  if (strcmp(host.name, "-") == 0) {
    ftServer = ft_start_buffer_server(host.port, NULL, NULL, NULL);
    if (ftServer==NULL) {
      fprintf(stderr, "unicorn2ft: could not start up a local buffer serving at port %i\n", host.port);
      return 1;
    }
    ftSocket = 0;
    printf("unicorn2ft: streaming to local buffer on port %i\n", host.port);
  }
  else {
    ftSocket = open_connection(host.name, host.port);

    if (ftSocket < 0) {
      fprintf(stderr, "unicorn2ft: could not connect to remote buffer at %s:%i\n", host.name, host.port);
      return 1;
    }
    printf("unicorn2ft: streaming to remote buffer at %s:%i\n", host.name, host.port);
  }

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

  /* define the header */
  header->def->nchans    = nchans;
  header->def->fsample   = fsample;
  header->def->nsamples  = 0;
  header->def->nevents   = 0;
  header->def->data_type = DATATYPE_FLOAT32;
  header->def->bufsize   = 0;

  /* define the constant part of the data and allocate space for the variable part */
  data->def->nchans = nchans;
  data->def->nsamples = 1;
  data->def->data_type = DATATYPE_FLOAT32;
  data->def->bufsize = WORDSIZE_FLOAT32 * nchans * 1;
  data->buf = malloc (data->def->bufsize);

  /* initialization phase, send the header */
  request->def->command = PUT_HDR;
  request->def->bufsize = append (&request->buf, request->def->bufsize, header->def, sizeof (headerdef_t));
  request->def->bufsize = append (&request->buf, request->def->bufsize, header->buf, header->def->bufsize);

  /* this is not needed any more */
  cleanup_header (&header);

  status = clientrequest (ftSocket, request, &response);
  fprintf (stderr, "unicorn2ft: clientrequest returned %d\n", status);

  if (status) {
    fprintf (stderr, "unicorn2ft: could not send request to buffer\n");
    exit (1);
  }

  if (status || response == NULL || response->def == NULL) {
    fprintf (stderr, "unicorn2ft: error in %s on line %d\n", __FILE__, __LINE__);
    exit (1);
  }

  cleanup_message (&request);

  if (response->def->command != PUT_OK) {
    fprintf (stderr, "unicorn2ft: error in 'put header' request.\n");
    exit (1);
  }

  cleanup_message (&response);

  /* open the serial port */
  fprintf (stderr, "unicorn2ft: opening serial port ...\n");
  if (!serialOpenByName (&SP, SerialDevice)) {

    fprintf (stderr, "Could not open serial port %s\n", SerialDevice);
    return 1;
  }

  if (!serialSetParameters (&SP, 115200, 8, 0, 0, 0)) {
    fprintf (stderr, "Could not modify serial port parameters\n");
    return 1;
  }

  fprintf (stderr, "unicorn2ft: opening serial port ... ok\n");

  fprintf (stderr, "unicorn2ft: initializing ... ok\n");

  printf ("Starting to listen - press CTRL-C to quit\n");

  /* register CTRL-C handler */
  signal (SIGINT, abortHandler);

  printf("sizeof(StartAcquisition) = %d\n", sizeof(StartAcquisition));

  /* start streaming data */
  n = serialWrite (&SP, 3, StartAcquisition);
  printf("n = %d\n", n);
    
  /* wait for the response */
  unsigned char buf[3];

  c = 0;
  while (c != 3) {
    usleep (1000);
    n = serialRead (&SP, 1, buf+c);
    c += n;
  }

  while (keepRunning) {
    usleep(1000000);
    printf (".\n");
    
    unsigned long blocksize = 44;
    unsigned char block[blocksize];

  }	/* while keepRunning */

  /* stop streaming data */
  serialWrite (&SP, 3, StopAcquisition);

  cleanup_data (&data);

  if (ftSocket > 0)
    close_connection (ftSocket);
  else
    ft_stop_buffer_server (ftServer);

  return 0;
}	/* main */
