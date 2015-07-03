/*
 * Copyright (C) 2015, Robert Oostenveld
 *
 * Donders Institute for Brain, Cognition and Behaviour,
 * Centre for Cognitive Neuroimaging, Radboud University,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 */

#include <stdio.h>
#include <sys/types.h>
#include <string.h>
#include <stdlib.h>
/* #include <unistd.h> */
#include <math.h>
#include <signal.h>

#include "serial.h"
#include "buffer.h"
#include "socketserver.h"

#define FTHOST          "-"
#define FTPORT          1972

#define BLOCKSIZE       50
#define BUFLEN          33

#define OPENBCI_NCHANS  11 /* 8x EEG, 3x accelerometer */
#define OPENBCI_FSAMPLE 250
#define OPENBCI_CALIB1  (1000000 * 4.5 / 24 / (2^23-1)) /* in uV, for 24x gain */
#define OPENBCI_CALIB2  0.002 / (2^4)                   /* in mG */

static char usage[] =
"\n" \
    "Usage: openbci2ft <device> [ftHost] [ftPort]\n" \
    "\n" \
    "When ftPort is omitted, it will default to 1972.\n" \
    "When ftHost is omitted, it will default to '-'.\n" \
    "Using '-' for the buffer hostname (ftHost) starts a local buffer on the given port (ftPort).\n" \
    "\n" \
    "Example use:\n" \
#ifdef (PLATFORM_WIN32)
    "   openbci2ft COM3:                 # start a local buffer on port 1972\n" \
    "   openbci2ft COM3: - 1234          # start a local buffer on port 1234\n" \
    "   openbci2ft COM3: mentat002 1234  # connect to remote buffer\n" \
#else
    "   openbci2ft /dev/tty.usbserial-DN0094FY                 # start a local buffer on port 1972\n" \
    "   openbci2ft /dev/tty.usbserial-DN0094FY - 1234          # start a local buffer on port 1234\n" \
    "   openbci2ft /dev/tty.usbserial-DN0094FY mentat002 1234  # connect to remote buffer\n" \
#endif
    "\n" \
    ;  

int keepRunning = 1;

void abortHandler(int sig) {
  printf("Ctrl-C pressed -- stopping...\n");
  keepRunning = 0;
}

int main(int argc, char *argv[])
{
  int n, i, c, sample = 0, status = 0, verbose = 0;
  unsigned char buf[BUFLEN], byte;
  SerialPort SP;
  host_t host;

  /* these represent the acquisition system properties */
  int nchans         = OPENBCI_NCHANS;
  int blocksize      = BLOCKSIZE;
  float fsample      = OPENBCI_FSAMPLE;

  /* these are used in the communication with the FT buffer and represent statefull information */
  int ftSocket           = -1;
  ft_buffer_server_t *ftServer;
  message_t    *request  = NULL;
  message_t    *response = NULL;
  header_t     *header   = NULL;
  data_t       *data     = NULL;

  if (argc<2) {
    printf(usage);
    exit(0);
  }

  if (argc>2)
    strcpy(host.name, argv[2]);
  else {
    strcpy(host.name, FTHOST);
  }

  if (argc>3)
    host.port = atoi(argv[3]);
  else {
    host.port = FTPORT;
  }

  fprintf(stderr, "openbci2ft: device       =  %s\n", argv[1]);
  fprintf(stderr, "openbci2ft: hostname     =  %s\n", host.name);
  fprintf(stderr, "openbci2ft: port         =  %d\n", host.port);

  /* Spawn tcpserver or connect to remote buffer */
  if (strcmp(host.name, "-") == 0) {
    ftServer = ft_start_buffer_server(host.port, NULL, NULL, NULL);
    if (ftServer==NULL) {
      fprintf(stderr, "openbci2ft: could not start up a local buffer serving at port %i\n", host.port);
      return 1;
    }
    ftSocket = 0;
    printf("openbci2ft: streaming to local buffer on port %i\n", host.port);
  }
  else {
    ftSocket = open_connection(host.name, host.port);

    if (ftSocket < 0) {
      fprintf(stderr, "openbci2ft: could not connect to remote buffer at %s:%i\n", host.name, host.port);
      return 1;
    }
    printf("openbci2ft: streaming to remote buffer at %s:%i\n", host.name, host.port);
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
  data->def->nchans    = nchans;
  data->def->nsamples  = blocksize;
  data->def->data_type = DATATYPE_FLOAT32;
  data->def->bufsize   = WORDSIZE_FLOAT32*nchans*blocksize;
  data->buf            = malloc(data->def->bufsize);

  /* initialization phase, send the header */
  request->def->command = PUT_HDR;
  request->def->bufsize = append(&request->buf, request->def->bufsize, header->def, sizeof(headerdef_t));
  request->def->bufsize = append(&request->buf, request->def->bufsize, header->buf, header->def->bufsize);

  /* this is not needed any more */
  cleanup_header(&header);

  status = clientrequest(ftSocket, request, &response);
  if (verbose>0)
    fprintf(stderr, "openbci2ft: clientrequest returned %d\n", status);
  if (status) {
    fprintf(stderr, "openbci2ft: could not send request to buffer\n");
    exit(1);
  }

  if (status || response==NULL || response->def == NULL) {
    fprintf(stderr, "openbci2ft: error in %s on line %d\n", __FILE__, __LINE__);
    exit(1);
  }

  cleanup_message(&request);

  if (response->def->command != PUT_OK) {
    fprintf(stderr, "openbci2ft: error in 'put header' request.\n");
    exit(1);
  }

  cleanup_message(&response);

  /* open the serial port */
  fprintf(stderr, "openbci2ft: opening serial port ...\n");
  if (!serialOpenByName(&SP, argv[1])) {
    fprintf(stderr, "Could not open serial port %s\n", argv[1]);
    return 1;
  }

  if (!serialSetParameters(&SP, 115200, 8, 0, 0, 0)) {
    fprintf(stderr, "Could not modify serial port parameters\n");
    return 1;
  }

  fprintf(stderr, "openbci2ft: opening serial port ... ok\n");

  /* 8-bit board will always be initialized upon opening serial port, 32-bit board needs explicit initialization */
  fprintf(stderr, "openbci2ft: initializing ...\n");

  serialWrite(&SP, 1, "v");
  fprintf(stderr, "openbci2ft: press reset on the OpenBCI board if this takes too long\n");
  usleep(1000000);

  /* wait for '$$$' which indicates that the OpenBCI has been initialized */
  c = 0;
  while (c!=3) {
    n = serialRead(&SP, 1, &byte);
    if (n==1) {
      if (byte=='$')
        c++;
      else
        c = 0;
    }
  } /* while waiting for '$$$' */

  fprintf(stderr, "openbci2ft: initializing ... ok\n");

  printf("Starting to listen - press CTRL-C to quit\n");

  /* register CTRL-C handler */
  signal(SIGINT, abortHandler);

  /* start streaming data */
  serialWrite(&SP, 1, "b");

  while (keepRunning) {

    c = 0;
    while (c<blocksize) {
      /* wait for the first byte of the packet */
      buf[0]=0;
      while (buf[0]!=0xA0) {
        n = serialRead(&SP, 1, buf);
      } /* while */

      /* read the remaining 32 bytes of the packet */
      while (n<BUFLEN)
        n += serialRead(&SP, (BUFLEN-n), buf+n);

      if (verbose>1) {
        for (i=0; i<BUFLEN; i++)
          printf("%02x ", buf[i]);
        printf("\n");
      }

      ((FLOAT32_T *)(data->buf))[nchans*c + 0] = OPENBCI_CALIB1 * (buf[ 2]<<24 | buf[ 3]<<16 | buf[ 4]<<8)/255;
      ((FLOAT32_T *)(data->buf))[nchans*c + 1] = OPENBCI_CALIB1 * (buf[ 5]<<24 | buf[ 6]<<16 | buf[ 7]<<8)/255;
      ((FLOAT32_T *)(data->buf))[nchans*c + 2] = OPENBCI_CALIB1 * (buf[ 8]<<24 | buf[ 9]<<16 | buf[10]<<8)/255;
      ((FLOAT32_T *)(data->buf))[nchans*c + 3] = OPENBCI_CALIB1 * (buf[11]<<24 | buf[12]<<16 | buf[13]<<8)/255;
      ((FLOAT32_T *)(data->buf))[nchans*c + 4] = OPENBCI_CALIB1 * (buf[14]<<24 | buf[15]<<16 | buf[16]<<8)/255;
      ((FLOAT32_T *)(data->buf))[nchans*c + 5] = OPENBCI_CALIB1 * (buf[17]<<24 | buf[18]<<16 | buf[19]<<8)/255;
      ((FLOAT32_T *)(data->buf))[nchans*c + 6] = OPENBCI_CALIB1 * (buf[20]<<24 | buf[21]<<16 | buf[22]<<8)/255;
      ((FLOAT32_T *)(data->buf))[nchans*c + 7] = OPENBCI_CALIB1 * (buf[23]<<24 | buf[24]<<16 | buf[25]<<8)/255;

      ((FLOAT32_T *)(data->buf))[nchans*c + 8] = OPENBCI_CALIB2 * (buf[26]<<24 | buf[27]<<16)/32767;
      ((FLOAT32_T *)(data->buf))[nchans*c + 9] = OPENBCI_CALIB2 * (buf[28]<<24 | buf[29]<<16)/32767;
      ((FLOAT32_T *)(data->buf))[nchans*c +10] = OPENBCI_CALIB2 * (buf[28]<<24 | buf[31]<<16)/32767;

      c++;
    } /* while c<blocksize */

    sample += blocksize;
    printf("openbci2ft: sample count = %i\n", sample);

    /*
     * Header
     *   Byte 1: 0xA0
     *   Byte 2: Sample Number
     *
     * EEG Data
     * Note: values are 24-bit signed, MSB first
     *   Bytes 3-5: Data value for EEG channel 1
     *   Bytes 6-8: Data value for EEG channel 2
     *   Bytes 9-11: Data value for EEG channel 3
     *   Bytes 12-14: Data value for EEG channel 4
     *   Bytes 15-17: Data value for EEG channel 5
     *   Bytes 18-20: Data value for EEG channel 6
     *   Bytes 21-23: Data value for EEG channel 6
     *   Bytes 24-26: Data value for EEG channel 8
     *
     * Accelerometer Data
     * Note: values are 16-bit signed, MSB first
     *   Bytes 27-28: Data value for accelerometer channel X
     *   Bytes 29-30: Data value for accelerometer channel Y
     *   Bytes 31-32: Data value for accelerometer channel Z
     *
     * Footer
     *   Byte 33: 0xC0
     */

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
    if (verbose>0)
      fprintf(stderr, "openbci2ft: clientrequest returned %d\n", status);
    if (status) {
      fprintf(stderr, "openbci2ft: error in %s on line %d\n", __FILE__, __LINE__);
      exit(1);
    }

    if (status) {
      fprintf(stderr, "openbci2ft: error in %s on line %d\n", __FILE__, __LINE__);
      exit(1);
    }

    /* FIXME do someting with the response, i.e. check that it is OK */
    cleanup_message(&request);

    if (response == NULL || response->def == NULL || response->def->command!=PUT_OK) {
      fprintf(stderr, "Error when writing samples.\n");
    }
    cleanup_message(&response);

  } /* while keepRunning */

  /* stop streaming data */
  serialWrite(&SP, 1, "s");

  cleanup_data(&data);

  if (ftSocket > 0) {
    close_connection(ftSocket);
  } else {
    ft_stop_buffer_server(ftServer);
  }
  return 0;
} /* main */

