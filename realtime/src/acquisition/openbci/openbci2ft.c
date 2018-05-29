/*
 * Copyright (C) 2015-2017, Robert Oostenveld
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
#include "ini.h"
#include "timestamp.h"

#define OPENBCI_BUFLEN    33
#define OPENBCI_CALIB1    (4.5   / (24 * (2^23-1))) /* in uV, assuming 24x gain */
#define OPENBCI_CALIB2    (0.002 / (2^4))           /* in mG */

int keepRunning = 1;

typedef struct {
  int        blocksize;
  int        port;
  const char *hostname;
  const char *serial;
  const char *reset;
  const char *datalog;
  const char *testsignal;
  const char *sample;
  const char *unwrap;
  const char *timestamp;
  const char *timeref;
  const char *acceleration;
  const char *daisy;
  int        verbose;

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
  const char *label_chan17; /* AccelerationX */
  const char *label_chan18; /* AccelerationY */
  const char *label_chan19; /* AccelerationZ */
  const char *label_chan20; /* this is for the sample    */
  const char *label_chan21; /* this is for the timestamp */

  const char *setting_chan1;
  const char *setting_chan2;
  const char *setting_chan3;
  const char *setting_chan4;
  const char *setting_chan5;
  const char *setting_chan6;
  const char *setting_chan7;
  const char *setting_chan8;
  const char *setting_chan9;
  const char *setting_chan10;
  const char *setting_chan11;
  const char *setting_chan12;
  const char *setting_chan13;
  const char *setting_chan14;
  const char *setting_chan15;
  const char *setting_chan16;

} configuration;

#if defined (PLATFORM_WIN32)
static char usage[] =
"\n" \
    "Use as\n" \
    "   openbci2ft <device> [ftHost] [ftPort]\n" \
    "with the parameters as specified below, or\n" \
    "   openbci2ft <config>\n" \
    "with the name of a configuration file for detailed setup.\n"
    "\n" \
      "When ftPort is omitted, it will default to 1972.\n" \
      "When ftHost is omitted, it will default to '-'.\n" \
      "Using '-' for the buffer hostname (ftHost) starts a local buffer on the given port (ftPort).\n" \
      "\n" \
      "Example use:\n" \
      "   openbci2ft COM3:                 # start a local buffer on port 1972\n" \
      "   openbci2ft COM3: - 1234          # start a local buffer on port 1234\n" \
      "   openbci2ft COM3: serverpc 1234   # connect to remote buffer running on server PC\n" \
      "\n" \
      ;
#else
static char usage[] =
"\n" \
    "Use as\n" \
    "   openbci2ft <device> [ftHost] [ftPort]\n" \
    "with the parameters as specified below, or\n" \
    "   openbci2ft <config>\n" \
    "with the name of a configuration file for detailed setup.\n"
    "\n" \
      "When ftPort is omitted, it will default to 1972.\n" \
      "When ftHost is omitted, it will default to '-'.\n" \
      "Using '-' for the buffer hostname (ftHost) starts a local buffer on the given port (ftPort).\n" \
      "\n" \
      "Example use:\n" \
      "   openbci2ft /dev/tty.usbserial-DN0094FY                 # start a local buffer on port 1972\n" \
      "   openbci2ft /dev/tty.usbserial-DN0094FY - 1234          # start a local buffer on port 1234\n" \
      "   openbci2ft /dev/tty.usbserial-DN0094FY serverpc 1234   # connect to remote buffer running on server PC\n" \
      "\n" \
      ;
#endif

int interpret24bitAsInt32(char *byteArray) {
    int newInt = (
     ((0xFF & byteArray[0]) << 16) |
     ((0xFF & byteArray[1]) << 8) |
     (0xFF & byteArray[2])
    );
    if ((newInt & 0x00800000) > 0) {
      newInt |= 0xFF000000;
    } else {
      newInt &= 0x00FFFFFF;
    }
    return newInt;
}

int interpret16bitAsInt32(char *byteArray) {
    int newInt = (
      ((0xFF & byteArray[0]) << 8) |
       (0xFF & byteArray[1])
      );
    if ((newInt & 0x00008000) > 0) {
          newInt |= 0xFFFF0000;
    } else {
          newInt &= 0x0000FFFF;
    }
    return newInt;
  }

int serialWriteSlow(SerialPort *SP, int size, void *buffer) {
  int i, retval = 0;
  for (i=0; i<size; i++) {
    retval += serialWrite(SP, 1, buffer+i);
    usleep(100000);
  }
  return retval;
}

static int iniHandler(void* external, const char* section, const char* name, const char* value) {
  configuration *local = (configuration*)external;

#define MATCH(s, n) strcasecmp(section, s) == 0 && strcasecmp(name, n) == 0
  if (MATCH("General", "blocksize")) {
    local->blocksize = atoi(value);
  } else if (MATCH("General", "port")) {
    local->port = atoi(value);
  } else if (MATCH("General", "hostname")) {
    local->hostname = strdup(value);
  } else if (MATCH("General", "serial")) {
    local->serial = strdup(value);
  } else if (MATCH("General", "reset")) {
    local->reset = strdup(value);
  } else if (MATCH("General", "datalog")) {
    local->datalog = strdup(value);
  } else if (MATCH("General", "testsignal")) {
    local->testsignal = strdup(value);
  } else if (MATCH("General", "sample")) {
    local->sample = strdup(value);
  } else if (MATCH("General", "unwrap")) {
    local->unwrap = strdup(value);
  } else if (MATCH("General", "timestamp")) {
    local->timestamp = strdup(value);
  } else if (MATCH("General", "timeref")) {
    local->timeref = strdup(value);
  } else if (MATCH("General", "acceleration")) {
    local->acceleration = strdup(value);
  } else if (MATCH("General", "daisy")) {
    local->daisy = strdup(value);
  } else if (MATCH("General", "verbose")) {
    local->verbose = atoi(value);

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
    local->enable_chan14 = strdup(value);
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

    /* channel names for acceleration, sample and timestamp */
  } else if (MATCH("ChannelLabel", "chan17")) {
    local->label_chan17 = strdup(value);
  } else if (MATCH("ChannelLabel", "chan18")) {
    local->label_chan18 = strdup(value);
  } else if (MATCH("ChannelLabel", "chan19")) {
    local->label_chan19 = strdup(value);
  } else if (MATCH("ChannelLabel", "chan20")) {
    local->label_chan20 = strdup(value);
  } else if (MATCH("ChannelLabel", "chan21")) {
    local->label_chan21 = strdup(value);

  } else if (MATCH("ChannelSetting", "chan1")) {
    local->setting_chan1 = strdup(value);
  } else if (MATCH("ChannelSetting", "chan2")) {
    local->setting_chan2 = strdup(value);
  } else if (MATCH("ChannelSetting", "chan3")) {
    local->setting_chan3 = strdup(value);
  } else if (MATCH("ChannelSetting", "chan4")) {
    local->setting_chan4 = strdup(value);
  } else if (MATCH("ChannelSetting", "chan5")) {
    local->setting_chan5 = strdup(value);
  } else if (MATCH("ChannelSetting", "chan6")) {
    local->setting_chan6 = strdup(value);
  } else if (MATCH("ChannelSetting", "chan7")) {
    local->setting_chan7 = strdup(value);
  } else if (MATCH("ChannelSetting", "chan8")) {
    local->setting_chan8 = strdup(value);
  } else if (MATCH("ChannelSetting", "chan9")) {
    local->setting_chan9 = strdup(value);
  } else if (MATCH("ChannelSetting", "chan10")) {
    local->setting_chan10 = strdup(value);
  } else if (MATCH("ChannelSetting", "chan11")) {
    local->setting_chan11 = strdup(value);
  } else if (MATCH("ChannelSetting", "chan12")) {
    local->setting_chan12 = strdup(value);
  } else if (MATCH("ChannelSetting", "chan13")) {
    local->setting_chan13 = strdup(value);
  } else if (MATCH("ChannelSetting", "chan14")) {
    local->setting_chan14 = strdup(value);
  } else if (MATCH("ChannelSetting", "chan15")) {
    local->setting_chan15 = strdup(value);
  } else if (MATCH("ChannelSetting", "chan16")) {
    local->setting_chan16 = strdup(value);

  } else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}

void abortHandler(int sig) {
  printf("Ctrl-C pressed -- stopping...\n");
  keepRunning = 0;
}

int main(int argc, char *argv[]) {
  int n, i, c, prevsample = 0, thissample = 0, sample = 0, count = 0, chan = 0, status = 0, labelSize;
  unsigned char buf[OPENBCI_BUFLEN], byte;
  char *labelString;
  SerialPort SP;
  host_t host;
  struct timespec tic, toc;

  /* these represent the general acquisition system properties */
  int nchans         = -1;
  float fsample      = -1;

  /* these are used in the communication with the FT buffer and represent statefull information */
  int ftSocket           = -1;
  ft_buffer_server_t *ftServer;
  message_t     *request  = NULL;
  message_t     *response = NULL;
  header_t      *header   = NULL;
  data_t        *data     = NULL;
  ft_chunkdef_t *label    = NULL;

  /* this structure contains the configuration details */
  configuration config;

  /* configure the default settings */
  config.blocksize     = 10;
  config.port          = 1972;
  config.hostname      = strdup("-");
  config.serial        = strdup("/dev/tty.usbserial-DN0094FY");
  config.reset         = strdup("on");
  config.datalog       = strdup("off");
  config.testsignal    = strdup("off");
  config.sample        = strdup("on");
  config.unwrap        = strdup("off");
  config.timestamp     = strdup("off");
  config.acceleration  = strdup("on");
  config.timeref       = strdup("start");
  config.daisy         = strdup("off");
  config.verbose       = 0;

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

  config.label_chan1  = strdup("ADC1");
  config.label_chan2  = strdup("ADC2");
  config.label_chan3  = strdup("ADC3");
  config.label_chan4  = strdup("ADC4");
  config.label_chan5  = strdup("ADC5");
  config.label_chan6  = strdup("ADC6");
  config.label_chan7  = strdup("ADC7");
  config.label_chan8  = strdup("ADC8");
  config.label_chan9  = strdup("ADC9");
  config.label_chan10 = strdup("ADC10");
  config.label_chan11 = strdup("ADC11");
  config.label_chan12 = strdup("ADC12");
  config.label_chan13 = strdup("ADC13");
  config.label_chan14 = strdup("ADC14");
  config.label_chan15 = strdup("ADC15");
  config.label_chan16 = strdup("ADC16");
  config.label_chan17 = strdup("AccelerationX");
  config.label_chan18 = strdup("AccelerationY");
  config.label_chan19 = strdup("AccelerationZ");
  config.label_chan20 = strdup("Sample");     /* virtual channel */
  config.label_chan21 = strdup("TimeStamp");  /* virtual channel */

  config.setting_chan1  = strdup("x1060101X");
  config.setting_chan2  = strdup("x2060101X");
  config.setting_chan3  = strdup("x3060101X");
  config.setting_chan4  = strdup("x4060101X");
  config.setting_chan5  = strdup("x5060101X");
  config.setting_chan6  = strdup("x6060101X");
  config.setting_chan7  = strdup("x7060101X");
  config.setting_chan8  = strdup("x8060101X");
  config.setting_chan9  = strdup("xQ060101X");
  config.setting_chan10 = strdup("xW060101X");
  config.setting_chan11 = strdup("xE060101X");
  config.setting_chan12 = strdup("xR060101X");
  config.setting_chan13 = strdup("xT060101X");
  config.setting_chan14 = strdup("xY060101X");
  config.setting_chan15 = strdup("xU060101X");
  config.setting_chan16 = strdup("xI060101X");

  if (argc<2) {
    printf("%s", usage);
    exit(0);
  }

  if (argc==2) {
    if (strncmp(argv[1], "/dev", 4)==0 || strncasecmp(argv[1], "COM", 3)==0)
      /* the second argument is the serial port */
      config.serial = strdup(argv[1]);
    else {
      /* the second argument is the configuration file */
      fprintf(stderr, "openbci2ft: loading configuration from '%s'\n", argv[1]);
      if (ini_parse(argv[1], iniHandler, &config) < 0) {
        fprintf(stderr, "Can't load '%s'\n", argv[1]);
        return 1;
      }
    }
  }

  if (argc>2)
    strcpy(host.name, argv[2]);
  else {
    strcpy(host.name, config.hostname);
  }

  if (argc>3)
    host.port = atoi(argv[3]);
  else {
    host.port = config.port;
  }

#define ISTRUE(s) strcasecmp(s, "on")==0

  if (ISTRUE(config.daisy)) {
    fsample = 125;
    /* have the user specify for all channels whether they are on or off */
  }
  else {
    fsample = 250;
    /* the default setting for all channels is on, switch channels 9-16 off */
    config.enable_chan9  = strdup("off");
    config.enable_chan10 = strdup("off");
    config.enable_chan11 = strdup("off");
    config.enable_chan12 = strdup("off");
    config.enable_chan13 = strdup("off");
    config.enable_chan14 = strdup("off");
    config.enable_chan15 = strdup("off");
    config.enable_chan16 = strdup("off");
  }

  /* count the number of channels that should be sent */
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
  if (ISTRUE(config.acceleration))
    nchans+=3;
  if (ISTRUE(config.sample))
    nchans++;
  if (ISTRUE(config.timestamp))
    nchans++;

  fprintf(stderr, "openbci2ft: serial       =  %s\n", config.serial);
  fprintf(stderr, "openbci2ft: hostname     =  %s\n", host.name);
  fprintf(stderr, "openbci2ft: port         =  %d\n", host.port);
  fprintf(stderr, "openbci2ft: nchans       =  %d\n", nchans);
  fprintf(stderr, "openbci2ft: fsample      =  %d\n", fsample);
  fprintf(stderr, "openbci2ft: blocksize    =  %d\n", config.blocksize);
  fprintf(stderr, "openbci2ft: reset        =  %s\n", config.reset);
  fprintf(stderr, "openbci2ft: datalog      =  %s\n", config.datalog);
  fprintf(stderr, "openbci2ft: sample       =  %s\n", config.sample);
  fprintf(stderr, "openbci2ft: timestamp    =  %s\n", config.timestamp);
  fprintf(stderr, "openbci2ft: acceleration =  %s\n", config.acceleration);
  fprintf(stderr, "openbci2ft: testsignal   =  %s\n", config.testsignal);
  fprintf(stderr, "openbci2ft: daisy        =  %s\n", config.daisy);

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
  if (ISTRUE (config.acceleration))
    labelSize += strlen (config.label_chan17) + 1;
  if (ISTRUE (config.acceleration))
    labelSize += strlen (config.label_chan18) + 1;
  if (ISTRUE (config.acceleration))
    labelSize += strlen (config.label_chan19) + 1;
  if (ISTRUE (config.sample))
    labelSize += strlen (config.label_chan20) + 1;
  if (ISTRUE (config.timestamp))
    labelSize += strlen (config.label_chan21) + 1;

  if (config.verbose > 2)
    fprintf (stderr, "openbci2ft: labelSize = %d\n", labelSize);

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
  if (ISTRUE (config.acceleration)) {
    strcpy (labelString+labelSize, config.label_chan17);
    labelSize += strlen (config.label_chan17) + 1;
  }
  if (ISTRUE (config.acceleration)) {
    strcpy (labelString+labelSize, config.label_chan18);
    labelSize += strlen (config.label_chan18) + 1;
  }
  if (ISTRUE (config.acceleration)) {
    strcpy (labelString+labelSize, config.label_chan19);
    labelSize += strlen (config.label_chan19) + 1;
  }
  if (ISTRUE (config.sample)) {
    strcpy (labelString+labelSize, config.label_chan20);
    labelSize += strlen (config.label_chan20) + 1;
  }
  if (ISTRUE (config.timestamp)) {
    strcpy (labelString+labelSize, config.label_chan21);
    labelSize += strlen (config.label_chan21) + 1;
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
  data->def->nchans = nchans;
  data->def->nsamples = config.blocksize;
  data->def->data_type = DATATYPE_FLOAT32;
  data->def->bufsize = WORDSIZE_FLOAT32 * nchans * config.blocksize;
  data->buf = malloc (data->def->bufsize);

  /* initialization phase, send the header */
  request->def->command = PUT_HDR;
  request->def->bufsize = append (&request->buf, request->def->bufsize, header->def, sizeof (headerdef_t));
  request->def->bufsize = append (&request->buf, request->def->bufsize, header->buf, header->def->bufsize);

  /* this is not needed any more */
  cleanup_header (&header);

  status = clientrequest (ftSocket, request, &response);
  if (config.verbose > 0)
    fprintf (stderr, "openbci2ft: clientrequest returned %d\n", status);
  if (status) {
    fprintf (stderr, "openbci2ft: could not send request to buffer\n");
    exit (1);
  }

  if (status || response == NULL || response->def == NULL) {
    fprintf (stderr, "openbci2ft: error in %s on line %d\n", __FILE__, __LINE__);
    exit (1);
  }

  cleanup_message (&request);

  if (response->def->command != PUT_OK) {
    fprintf (stderr, "openbci2ft: error in 'put header' request.\n");
    exit (1);
  }

  cleanup_message (&response);

  /* open the serial port */
  fprintf (stderr, "openbci2ft: opening serial port ...\n");
  if (!serialOpenByName (&SP, config.serial)) {

    fprintf (stderr, "Could not open serial port %s\n", config.serial);
    return 1;
  }

  if (!serialSetParameters (&SP, 115200, 8, 0, 0, 0)) {
    fprintf (stderr, "Could not modify serial port parameters\n");
    return 1;
  }

  fprintf (stderr, "openbci2ft: opening serial port ... ok\n");

  /* 8-bit board will always be initialized upon opening serial port, 32-bit board needs explicit initialization */

  fprintf (stderr, "openbci2ft: initializing ...\n");
  fprintf (stderr, "openbci2ft: press reset on the OpenBCI board if this takes too long\n");

  if (ISTRUE (config.reset))
    serialWrite (&SP, 1, "v");	/* soft reset, this will return $$$ */
  else
    serialWrite (&SP, 1, "D");	/* query default channel settings, this will also return $$$ */

  /* wait for '$$$' which indicates that the OpenBCI has been initialized */
  c = 0;
  while (c != 3) {
    usleep (1000);
    n = serialRead (&SP, 1, &byte);
    if (n == 1) {
      if (byte == '$')
        c++;
      else
        c = 0;
    }
  }				/* while waiting for '$$$' */

  if (strcasecmp (config.datalog, "14s") == 0)
    serialWrite (&SP, 1, "a");
  else if (strcasecmp (config.datalog, "5min") == 0)
    serialWrite (&SP, 1, "A");
  else if (strcasecmp (config.datalog, "15min") == 0)
    serialWrite (&SP, 1, "S");
  else if (strcasecmp (config.datalog, "30min") == 0)
    serialWrite (&SP, 1, "F");
  else if (strcasecmp (config.datalog, "1hr") == 0)
    serialWrite (&SP, 1, "G");
  else if (strcasecmp (config.datalog, "2hr") == 0)
    serialWrite (&SP, 1, "H");
  else if (strcasecmp (config.datalog, "4hr") == 0)
    serialWrite (&SP, 1, "J");
  else if (strcasecmp (config.datalog, "12hr") == 0)
    serialWrite (&SP, 1, "K");
  else if (strcasecmp (config.datalog, "24hr") == 0)
    serialWrite (&SP, 1, "L");
  else if (strcasecmp (config.datalog, "off") != 0) {
    fprintf (stderr, "Incorrect specification of datalog\n");
    return 1;
  }

  serialWriteSlow (&SP, strlen (config.setting_chan1), config.setting_chan1);
  serialWriteSlow (&SP, strlen (config.setting_chan2), config.setting_chan2);
  serialWriteSlow (&SP, strlen (config.setting_chan3), config.setting_chan3);
  serialWriteSlow (&SP, strlen (config.setting_chan4), config.setting_chan4);
  serialWriteSlow (&SP, strlen (config.setting_chan5), config.setting_chan5);
  serialWriteSlow (&SP, strlen (config.setting_chan6), config.setting_chan6);
  serialWriteSlow (&SP, strlen (config.setting_chan7), config.setting_chan7);
  serialWriteSlow (&SP, strlen (config.setting_chan8), config.setting_chan8);
  if (ISTRUE(config.daisy)) {
    serialWriteSlow (&SP, strlen (config.setting_chan9), config.setting_chan9);
    serialWriteSlow (&SP, strlen (config.setting_chan10), config.setting_chan10);
    serialWriteSlow (&SP, strlen (config.setting_chan11), config.setting_chan11);
    serialWriteSlow (&SP, strlen (config.setting_chan12), config.setting_chan12);
    serialWriteSlow (&SP, strlen (config.setting_chan13), config.setting_chan13);
    serialWriteSlow (&SP, strlen (config.setting_chan14), config.setting_chan14);
    serialWriteSlow (&SP, strlen (config.setting_chan15), config.setting_chan15);
    serialWriteSlow (&SP, strlen (config.setting_chan16), config.setting_chan16);
  }

  if (strcasecmp (config.testsignal, "gnd") == 0)
    serialWrite (&SP, 1, "0");
  else if (strcasecmp (config.testsignal, "dc") == 0)
    serialWrite (&SP, 1, "-");
  else if (strcasecmp (config.testsignal, "1xSlow") == 0)
    serialWrite (&SP, 1, "=");
  else if (strcasecmp (config.testsignal, "1xFast") == 0)
    serialWrite (&SP, 1, "p");
  else if (strcasecmp (config.testsignal, "2xSlow") == 0)
    serialWrite (&SP, 1, "[");
  else if (strcasecmp (config.testsignal, "2xFast") == 0)
    serialWrite (&SP, 1, "]");
  else if (strcasecmp (config.testsignal, "off") != 0) {
    fprintf (stderr, "Incorrect specification of testsignal\n");
    return 1;
  }

  fprintf (stderr, "openbci2ft: initializing ... ok\n");

  printf ("Starting to listen - press CTRL-C to quit\n");

  /* register CTRL-C handler */
  signal (SIGINT, abortHandler);

  /* start streaming data */
  serialWrite (&SP, 1, "b");

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

    sample = 0;
    chan = 0;

    while (sample < config.blocksize) {

      /* wait for the first byte of the following packet */
      buf[0] = 0;
      while (buf[0] != 0xA0) {
        if (serialInputPending (&SP))
          n = serialRead (&SP, 1, buf);
        else
          usleep (1000);
      }			/* while */

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
       *
       */

      /* read the remaining 32 bytes of the packet */
      while (n < OPENBCI_BUFLEN)
        if (serialInputPending (&SP))
          n += serialRead (&SP, (OPENBCI_BUFLEN - n), buf + n);
        else
          usleep (1000);

      if (config.verbose > 1) {
        for (i = 0; i < OPENBCI_BUFLEN; i++)
          printf ("%02x ", buf[i]);
        printf ("\n");
      }

      if (buf[32]!=0xC0)
        continue;

      /* these two booleans are both false if there is no daisy board */
      char lower = (ISTRUE(config.daisy) && (buf[1] % 2 == 1)); /* odd samples are for the lower 8 channels */
      char upper = (ISTRUE(config.daisy) && (buf[1] % 2 == 0)); /* even samples are for the upper 8 channels from the daisy board */

      if (lower || (!lower && !upper)) {
        /* assign the lower 8 channels */
        if (ISTRUE (config.enable_chan1))
          ((FLOAT32_T *) (data->buf))[nchans * sample + (chan++)] =
            OPENBCI_CALIB1 * interpret24bitAsInt32(buf+2);
        if (ISTRUE (config.enable_chan2))
          ((FLOAT32_T *) (data->buf))[nchans * sample + (chan++)] =
            OPENBCI_CALIB1 * interpret24bitAsInt32(buf+5);
        if (ISTRUE (config.enable_chan3))
          ((FLOAT32_T *) (data->buf))[nchans * sample + (chan++)] =
            OPENBCI_CALIB1 * interpret24bitAsInt32(buf+8);
        if (ISTRUE (config.enable_chan4))
          ((FLOAT32_T *) (data->buf))[nchans * sample + (chan++)] =
            OPENBCI_CALIB1 * interpret24bitAsInt32(buf+11);
        if (ISTRUE (config.enable_chan5))
          ((FLOAT32_T *) (data->buf))[nchans * sample + (chan++)] =
            OPENBCI_CALIB1 * interpret24bitAsInt32(buf+14);
        if (ISTRUE (config.enable_chan6))
          ((FLOAT32_T *) (data->buf))[nchans * sample + (chan++)] =
            OPENBCI_CALIB1 * interpret24bitAsInt32(buf+17);
        if (ISTRUE (config.enable_chan7))
          ((FLOAT32_T *) (data->buf))[nchans * sample + (chan++)] =
            OPENBCI_CALIB1 * interpret24bitAsInt32(buf+20);
        if (ISTRUE (config.enable_chan8))
          ((FLOAT32_T *) (data->buf))[nchans * sample + (chan++)] =
            OPENBCI_CALIB1 * interpret24bitAsInt32(buf+23);
      } /* if lower or no daisy board */

      if (upper) {
        /* assign the upper 8 channels of the daisy board */
        if (ISTRUE (config.enable_chan9))
          ((FLOAT32_T *) (data->buf))[nchans * sample + (chan++)] =
            OPENBCI_CALIB1 * interpret24bitAsInt32(buf+2);
        if (ISTRUE (config.enable_chan10))
          ((FLOAT32_T *) (data->buf))[nchans * sample + (chan++)] =
            OPENBCI_CALIB1 * interpret24bitAsInt32(buf+5);
        if (ISTRUE (config.enable_chan11))
          ((FLOAT32_T *) (data->buf))[nchans * sample + (chan++)] =
            OPENBCI_CALIB1 * interpret24bitAsInt32(buf+8);
        if (ISTRUE (config.enable_chan12))
          ((FLOAT32_T *) (data->buf))[nchans * sample + (chan++)] =
            OPENBCI_CALIB1 * interpret24bitAsInt32(buf+11);
        if (ISTRUE (config.enable_chan13))
          ((FLOAT32_T *) (data->buf))[nchans * sample + (chan++)] =
            OPENBCI_CALIB1 * interpret24bitAsInt32(buf+14);
        if (ISTRUE (config.enable_chan14))
          ((FLOAT32_T *) (data->buf))[nchans * sample + (chan++)] =
            OPENBCI_CALIB1 * interpret24bitAsInt32(buf+17);
        if (ISTRUE (config.enable_chan15))
          ((FLOAT32_T *) (data->buf))[nchans * sample + (chan++)] =
            OPENBCI_CALIB1 * interpret24bitAsInt32(buf+20);
        if (ISTRUE (config.enable_chan16))
          ((FLOAT32_T *) (data->buf))[nchans * sample + (chan++)] =
            OPENBCI_CALIB1 * interpret24bitAsInt32(buf+23);
      } /* if upper */

      if (lower || (!lower && !upper)) {
        /* acceleration channels */
        if (ISTRUE (config.acceleration))
          ((FLOAT32_T *) (data->buf))[nchans * sample + (chan++)] =
            OPENBCI_CALIB2 * interpret16bitAsInt32(buf+26);
        if (ISTRUE (config.acceleration))
          ((FLOAT32_T *) (data->buf))[nchans * sample + (chan++)] =
            OPENBCI_CALIB2 * interpret16bitAsInt32(buf+28);
        if (ISTRUE (config.acceleration))
          ((FLOAT32_T *) (data->buf))[nchans * sample + (chan++)] =
            OPENBCI_CALIB2 * interpret16bitAsInt32(buf+30);

        /* virtual channel with sample numbers */
        if (ISTRUE (config.sample)) {
          if (ISTRUE (config.unwrap)) {
            thissample += (buf[1] < prevsample ? 256 + buf[1]-prevsample : buf[1]-prevsample);
            prevsample  = buf[1];
          }
          else {
            thissample = buf[1];
          }
          ((FLOAT32_T *) (data->buf))[nchans * sample + (chan++)] = thissample;
        }

        /* virtual channel with timestamps */
        if (ISTRUE (config.timestamp)) {
          if (strcasecmp (config.timeref, "start") == 0)
            get_monotonic_time (&toc, TIMESTAMP_REF_BOOT);
          else if (strcasecmp (config.timeref, "boot") == 0)
            get_monotonic_time (&toc, TIMESTAMP_REF_BOOT);
          else if (strcasecmp (config.timeref, "epoch") == 0)
            get_monotonic_time (&toc, TIMESTAMP_REF_EPOCH);
          ((FLOAT32_T *) (data->buf))[nchans * sample + (chan++)] =
            get_elapsed_time (&tic, &toc);
        }
      } /* if lower or no daisy board */

      if (upper || (!lower && !upper)) {
        sample++;              /* move to the next column of the data matrix */
        chan = 0;              /* move to the first row of the data matrix */
      }
    }	/* while sample<config.blocksize */

    count += sample;
    printf ("openbci2ft: sample count = %i\n", count);

    /* create the request */
    request = malloc (sizeof (message_t));
    request->def = malloc (sizeof (messagedef_t));
    request->buf = NULL;
    request->def->version = VERSION;
    request->def->bufsize = 0;
    request->def->command = PUT_DAT;

    request->def->bufsize = append (&request->buf, request->def->bufsize, data->def, sizeof (datadef_t));
    request->def->bufsize = append (&request->buf, request->def->bufsize, data->buf, data->def->bufsize);

    status = clientrequest (ftSocket, request, &response);
    if (config.verbose > 0)
      fprintf (stderr, "openbci2ft: clientrequest returned %d\n", status);
    if (status) {
      fprintf (stderr, "openbci2ft: error in %s on line %d\n", __FILE__, __LINE__);
      exit (1);
    }

    if (status) {
      fprintf (stderr, "openbci2ft: error in %s on line %d\n", __FILE__, __LINE__);
      exit (1);
    }

    /* FIXME do someting with the response, i.e. check that it is OK */
    cleanup_message (&request);

    if (response == NULL || response->def == NULL || response->def->command != PUT_OK) {
      fprintf (stderr, "Error when writing samples.\n");
    }
    cleanup_message (&response);

  }	/* while keepRunning */

  /* stop streaming data */
  serialWrite (&SP, 1, "s");

  cleanup_data (&data);

  if (ftSocket > 0)
    close_connection (ftSocket);
  else
    ft_stop_buffer_server (ftServer);

  return 0;
}	/* main */
