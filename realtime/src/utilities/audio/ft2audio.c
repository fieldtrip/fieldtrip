/*
 * This is a ommand-line application to stream signals from a FieldTrip buffer
 * to the sound card. This is based on PortAudio. The FieldTrip buffer
 * should contain one or two channels at an appropriate sampling frequency.
 *
 * (C) 2017, Robert Oostenveld
*/

#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>

#include "portaudio.h"
#include "message.h"
#include "buffer.h"
#include "interface.h"
#include "socketserver.h"
#include "ini.h"

#define TRUE          (1)
#define MIN(x,y) 			(x<y ? x : y)
#define MAX(x,y) 			(x>y ? x : y)
#define MATCH(s, n) 	(strcasecmp(section, s) == 0 && strcasecmp(name, n) == 0)

typedef struct
{
	float *block1;
	float *block2;
	char refresh1;
	char refresh2;
	unsigned int current;
	unsigned int nsamples;
	float fsample;
	unsigned int nchans;
	unsigned int datatype;
	unsigned int sample;
	unsigned int nput;
	unsigned int nget;
	unsigned int skipped;
}
callbackData_t;

typedef struct
{
  int        port;
  int        verbose;
  const char *hostname;
	int				 monitor;
	int				 timeout;
	int				 slip;
	float  	   calibration;
	int        blocksize;
	int        callbacksize;
} configuration_t;

static char usage[] =
"\n" \
"Use as\n" \
"   ft2audio [ftHost] [ftPort]\n" \
"with the parameters as specified below, or\n" \
"   ft2audio <config>\n" \
"with the name of a configuration file for detailed setup.\n"
"\n" \
"Audio output is to the default system output device.\n" \
"When ftPort is omitted, it will default to 1972.\n" \
"When ftHost is omitted, it will default to '-'.\n" \
"Using '-' for the buffer hostname (ftHost) starts a local buffer on the given port (ftPort).\n" \
"\n" \
"Example use:\n" \
"   ft2audio                 # start a local buffer on port 1972\n" \
"   ft2audio - 1234          # start a local buffer on port 1234\n" \
"   ft2audio serverpc 1234   # connect to remote buffer running on server PC\n" \
"\n" \
;

/*******************************************************************/
static int iniHandler(void* external, const char* section, const char* name, const char* value)
{
  configuration_t *local = (configuration_t *)external;

  if (MATCH("General", "port")) {
    local->port = atoi(value);
  } else if (MATCH("General", "verbose")) {
    local->verbose = atoi(value);
	} else if (MATCH("General", "hostname")) {
    local->hostname = strdup(value);

	} else if (MATCH("General", "monitor")) {
		local->monitor = atoi(value);
	} else if (MATCH("General", "timeout")) {
		local->timeout = atoi(value);
	} else if (MATCH("General", "slip")) {
		local->slip = atoi(value);
	} else if (MATCH("General", "calibration")) {
		local->calibration = atof(value);
	} else if (MATCH("General", "blocksize")) {
		local->blocksize = atoi(value);
	} else if (MATCH("General", "callbacksize")) {
		local->callbacksize = atoi(value);

	} else {
    return 0;  /* unknown section/name, error */
  }
  return 1;
}


/*******************************************************************/
static int paWriteCallback(const void *inputBuffer,
		void *outputBuffer,
		unsigned long framesPerBuffer,
		const PaStreamCallbackTimeInfo* timeInfo,
		PaStreamCallbackFlags statusFlags,
		void *ptr)
{
	/* Prevent unused variable warning. */
	(void) inputBuffer;
	/* Cast data passed through stream to our structure. */
	callbackData_t *transfer = ptr;
	float *output = (float *)outputBuffer;

  /* this is to get started */
  if (transfer->current==0 && transfer->refresh1==0)
	  transfer->current = 1; /* switch to the first available block */
	else if (transfer->current==0 && transfer->refresh2==0)
		transfer->current = 2; /* switch to the second available block */

	while (framesPerBuffer>0) {
		framesPerBuffer--;

		if (transfer->current==1) {
			for (unsigned int chan=0; chan<transfer->nchans; chan++)
				output[chan] = transfer->block1[(transfer->sample)*transfer->nchans+chan];
			transfer->sample++;
			transfer->nget++;
		}
		else if (transfer->current==2) {
			for (unsigned int chan=0; chan<transfer->nchans; chan++)
				output[chan] = transfer->block2[(transfer->sample)*transfer->nchans+chan];
			transfer->sample++;
			transfer->nget++;
		}
		else {
			for (unsigned int chan=0; chan<transfer->nchans; chan++)
				output[chan] = 0;
			transfer->skipped++;
		}

		for (unsigned int chan=0; chan<transfer->nchans; chan++)
			output++;

		if (transfer->sample==transfer->nsamples) {
			if (transfer->current==1) {
				/* switch to the start of the other block */
				transfer->refresh1 = 1;
				transfer->current = (transfer->refresh2 ? 0 : 2);
				transfer->sample = 0;
			}
			else if (transfer->current==2) {
				/* switch to the start of the other block */
				transfer->refresh2 = 1;
				transfer->current = (transfer->refresh1 ? 0 : 1);
				transfer->sample = 0;
			}
			else {
				/* both blocks are not available */
				transfer->current = 0;
				transfer->sample = 0;
			}
		} /* if switch to other block */
	} /* while framesPerBuffer */

	return 0;
} /* paWriteCallback */

/*******************************************************************/
int file_exists(const char *fname) {
  /* see http://stackoverflow.com/questions/230062/whats-the-best-way-to-check-if-a-file-exists-in-c-cross-platform */
  if( access( fname, F_OK ) != -1 )
    return 1;
  else
    return 0;
}

/*******************************************************************/
void tic(struct timeval *now)
{
	/* Set a timer. */
  gettimeofday(now, NULL);
  return;
}

/*******************************************************************/
int toc(struct timeval previous)
{
   struct timeval now;
   gettimeofday(&now, NULL);
	 /* Return the elapsed time in milliseconds. */
   return (now.tv_sec - previous.tv_sec)*1000 + (now.tv_usec - previous.tv_usec)/1000;
}

/*******************************************************************/
int main(int argc, char *argv[])
{
	PaStream *stream;
	PaError err;
	callbackData_t transfer;
	struct timeval timeout, monitor;
	int ftSocket, status, retry;
	ft_buffer_server_t *ftServer;
	UINT32_T datatype;
	float fsample;
	unsigned int nchans, nsamples, nevents, begsample, endsample;
	void *bufdata;
  configuration_t config;
	host_t host;

  /* configure the default settings */
	config.verbose       = 1;
	config.hostname      = strdup("-");
  config.port          = 1972;
	config.monitor       = 2000;	/* in milliseconds */
	config.timeout       = 60000;	/* in milliseconds */
	config.slip          = 5000;	/* in milliseconds */
	config.blocksize     = 1000;	/* in milliseconds */
	config.callbacksize  = 100;		/* in milliseconds */
	config.calibration   = 1.0f;

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
		fprintf(stderr, "ft2audio: loading configuration from '%s'\n", argv[1]);
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

	/* Spawn tcpserver or connect to remote buffer */
  if (strcmp(host.name, "-") == 0) {
    ftServer = ft_start_buffer_server(host.port, NULL, NULL, NULL);
    if (ftServer==NULL) {
      fprintf(stderr, "ft2audio: could not start up a local buffer serving at port %i\n", host.port);
      return 1;
    }
    ftSocket = 0;
    printf("ft2audio: streaming from local buffer on port %i\n", host.port);
  }
  else {
		ftServer = 0;
    ftSocket = open_connection(host.name, host.port);

    if (ftSocket < 0) {
      fprintf(stderr, "ft2audio: could not connect to remote buffer at %s:%i\n", host.name, host.port);
      return 1;
    }
    printf("ft2audio: streaming from remote buffer at %s:%i\n", host.name, host.port);
  }

  /* Wait for the first data to arrive, especially important in case a local buffer is started. */
	status = 1;
	tic(&timeout);
	while (status && toc(timeout)<config.timeout) {
		if (config.verbose>0)
			printf("ft2audio: waiting for valid header, timeout in %d seconds\n", (config.timeout-toc(timeout))/1000);
		status = read_header(ftSocket, &datatype, &nchans, &fsample, &nsamples, &nevents);
		if (status)
		  Pa_Sleep(1000);
	}
	if (status) goto error;

	if (config.verbose>0) {
		printf("ft2audio: nchans   = %d\n", nchans);
		printf("ft2audio: fsample  = %.0f\n", fsample);
		printf("ft2audio: nsamples = %d\n", nsamples);
	}

	/* Initialize the structure that is shared with the callback. */
	transfer.refresh1 = 1; /* needs new data */
	transfer.refresh2 = 1; /* needs new data */
	transfer.current  = 0; /* neither one contains data */
	transfer.datatype = datatype;
	transfer.nchans   = nchans;
	transfer.fsample  = fsample;
	transfer.nsamples = fsample*config.blocksize/1000;
	transfer.sample   = 0;
	transfer.nput     = 0; /* how many samples are read from the data stream */
	transfer.nget     = 0; /* how many samples are written to the audio stream */
	transfer.skipped  = 0; /* how often was a sample skipped in the audio stream */
	transfer.block1   = malloc(transfer.nchans*transfer.nsamples*sizeof(float));
	transfer.block2   = malloc(transfer.nchans*transfer.nsamples*sizeof(float));
	DIE_BAD_MALLOC(transfer.block1);
	DIE_BAD_MALLOC(transfer.block2);
	bzero(transfer.block1, transfer.nchans*transfer.nsamples*sizeof(float));
	bzero(transfer.block2, transfer.nchans*transfer.nsamples*sizeof(float));

	/* Initialize to read the data from the buffer. */
	switch (datatype) {
		case DATATYPE_INT16:
			bufdata = malloc(transfer.nchans*transfer.nsamples*WORDSIZE_INT16);
			break;
		case DATATYPE_INT32:
			bufdata = malloc(transfer.nchans*transfer.nsamples*WORDSIZE_INT32);
			break;
		case DATATYPE_INT64:
			bufdata = malloc(transfer.nchans*transfer.nsamples*WORDSIZE_INT64);
			break;
		case DATATYPE_FLOAT32:
			bufdata = malloc(transfer.nchans*transfer.nsamples*WORDSIZE_FLOAT32);
			break;
		case DATATYPE_FLOAT64:
			bufdata = malloc(transfer.nchans*transfer.nsamples*WORDSIZE_FLOAT64);
			break;
		default:
			printf("ft2audio: unsupported datatype = %d\n", datatype);
			goto error;
	}
	DIE_BAD_MALLOC(bufdata);

	if (config.verbose>0) {
		printf( "ft2audio: PortAudio version = 0x%08X\n", Pa_GetVersion());
		printf( "ft2audio: Version text = '%s'\n", Pa_GetVersionInfo()->versionText );
	}

	/* Initialize library before making any other calls. */
	err = Pa_Initialize();
	if (err != paNoError) goto error;

	if (config.verbose>0)
		printf("ft2audio: initialized PortAudio\n");

	/* Open an audio I/O stream. */
	err = Pa_OpenDefaultStream(&stream, 0, nchans, paFloat32, fsample, fsample*config.callbacksize/1000, paWriteCallback, &transfer);
	if (err != paNoError) goto error;

	if (config.verbose>0)
		printf("ft2audio: opened PortAudio stream\n");

	err = Pa_StartStream(stream);
	if (err != paNoError) goto error;

	if (config.verbose>0)
		printf("ft2audio: started PortAudio stream\n");

	/* Jump to the end of the available data, the case of "no data" will be dealt with further down. */
	begsample = nsamples-transfer.nsamples;
	endsample = nsamples-1;

	/* initialize the timers */
	tic(&timeout);
	tic(&monitor);

	while (TRUE) {

		/* monitor the lead or lag between the input and output */
		if (toc(monitor)>config.monitor) {
			tic(&monitor); /* reset the monitor timer */

			status = read_header(ftSocket, &datatype, &nchans, &fsample, &nsamples, &nevents);
			if (status) goto error;

			if (config.verbose>1)
				printf("read_header: datatype = %d, nchans = %d, fsample = %.0f, nsamples = %d, nevents = %d\n", datatype, nchans, fsample, nsamples, nevents);

			/* these should not change once this application is running */
			if (datatype!=transfer.datatype) {
				fprintf(stderr, "Error: unexpected change in datatype\n");
				goto error;
			}
			if (nchans!=transfer.nchans) {
				fprintf(stderr, "Error: unexpected change in number of channels\n");
				goto error;
			}
			if (fsample!=transfer.fsample) {
				fprintf(stderr, "Error: unexpected change in sampling rate\n");
				goto error;
			}

			long data_lead = (long)endsample-(nsamples-1);
			long audio_lag = (long)transfer.nput-transfer.nget-transfer.skipped;
			char jump_to_last = 0, jump_to_next = 0;

			if (config.verbose>0)
				printf("ft2audio: data lead = %ld, audio lag = %ld, skipped = %d, read = %d\n", data_lead, audio_lag, transfer.skipped, transfer.nput);

			if (data_lead > +fsample*config.slip/1000) {
				/* reading data too far in the future, data is coming in too slow */
				printf("ft2audio: data lead too large, reverting to the end\n");
				jump_to_next = 1;
			}
			else if (data_lead < -fsample*config.slip/1000) {
				/* reading data too far in the past, playback is too slow */
				printf("ft2audio: data lag too large, skipping to the end\n");
				jump_to_last = 1;
			}
			else if (audio_lag > +fsample*config.slip/1000) {
				/* writing data slower than it is arriving, playback is too slow */
				printf("ft2audio: audio lag too large, skipping to the end\n");
				jump_to_last = 1;
			}
			else if (audio_lag < -fsample*config.slip/1000) {
				/* writing data faster than it is arriving, data is coming in too slow */
				printf("ft2audio: audio lead too large, reverting to the end\n");
				jump_to_next = 1;
			}

			if (jump_to_last) {
				/* jump to the last available block */
				begsample = nsamples-transfer.nsamples;
				endsample = nsamples-1;
				/* reset the audio stream counters */
				transfer.nput = 0;
				transfer.nget = 0;
				transfer.skipped = 0;
				printf("ft2audio: new endsample = %d\n", endsample);
			}
			else if (jump_to_next) {
				/* jump to the first non-available block */
				begsample = nsamples;
				endsample = nsamples+transfer.nsamples-1;
				/* reset the audio stream counters */
				transfer.nput = 0;
				transfer.nget = 0;
				transfer.skipped = 0;
				printf("ft2audio: new endsample = %d\n", endsample);
			}

		} /* if monitor */

		if (transfer.refresh1) {
			status = wait_data(ftSocket, endsample, 0, 1200.*transfer.nsamples/fsample);
			if (config.verbose>1)
				printf("wait_data: block = 1, endsample = %d, status = %d\n", endsample, status);
			if (status!=0)
			  break;
			status = read_data(ftSocket, begsample, endsample, bufdata);
			if (config.verbose>1)
				printf("read_data: block = 1, endsample = %d, status = %d\n", endsample, status);
			/* status ==0 when new data or !=0 when data not (yet) available */
			if (status==0) {
				tic(&timeout);
				begsample += transfer.nsamples; // in samples
				endsample += transfer.nsamples;

				for (unsigned int sample=0; sample<transfer.nsamples; sample++) {
					for (unsigned int chan=0; chan<transfer.nchans; chan++) {
						unsigned int i = (sample*transfer.nchans)+chan;
						switch (datatype) {
							case DATATYPE_INT16:
								transfer.block1[i] = ((INT16_T *)bufdata)[i];
								break;
							case DATATYPE_INT32:
								transfer.block1[i] = ((INT32_T *)bufdata)[i];
								break;
							case DATATYPE_INT64:
								transfer.block1[i] = ((INT64_T *)bufdata)[i];
								break;
							case DATATYPE_FLOAT32:
								transfer.block1[i] = ((FLOAT32_T *)bufdata)[i];
								break;
							case DATATYPE_FLOAT64:
								transfer.block1[i] = ((FLOAT64_T *)bufdata)[i];
								break;
							default:
								printf("ft2audio: unsupported datatype = %d\n", datatype);
								goto error;
						}
						transfer.block1[i] *= config.calibration;
					}
				}
				transfer.refresh1 = 0;
				transfer.nput += transfer.nsamples;
			}
			else {
				printf("ft2audio: reading block2 failed, timeout in %d seconds\n", (config.timeout-toc(timeout))/1000);
			}
		}

		if (transfer.refresh2) {
			status = wait_data(ftSocket, endsample, 0, 1200.*transfer.nsamples/fsample);
			if (config.verbose>1)
				printf("wait_data: block = 2, endsample = %d, status = %d\n", endsample, status);
			if (status!=0)
			  break;
			status = read_data(ftSocket, begsample, endsample, bufdata);
			if (config.verbose>1)
				printf("read_data: block = 2, endsample = %d, status = %d\n", endsample, status);
			/* status ==0 when new data or !=0 when data not (yet) available */
			if (status==0) {
				tic(&timeout);
				begsample += transfer.nsamples; // in samples
				endsample += transfer.nsamples;

				for (unsigned int sample=0; sample<transfer.nsamples; sample++) {
					for (unsigned int chan=0; chan<transfer.nchans; chan++) {
						unsigned int i = (sample*transfer.nchans)+chan;
						switch (datatype) {
							case DATATYPE_INT16:
								transfer.block2[i] = ((INT16_T *)bufdata)[i];
								break;
							case DATATYPE_INT32:
								transfer.block2[i] = ((INT32_T *)bufdata)[i];
								break;
							case DATATYPE_INT64:
								transfer.block2[i] = ((INT64_T *)bufdata)[i];
								break;
							case DATATYPE_FLOAT32:
								transfer.block2[i] = ((FLOAT32_T *)bufdata)[i];
								break;
							case DATATYPE_FLOAT64:
								transfer.block2[i] = ((FLOAT64_T *)bufdata)[i];
								break;
							default:
								printf("ft2audio: unsupported datatype = %d\n", datatype);
								goto error;
						}
						transfer.block2[i] *= config.calibration;
					}
				}
				transfer.refresh2 = 0;
				transfer.nput += transfer.nsamples;
			}
			else {
				printf("ft2audio: reading block2 failed, timeout in %d seconds\n", (config.timeout-toc(timeout))/1000);
			}
		}

		if (transfer.refresh1==0 && transfer.refresh2==0) {
			/* Both blocks have fresh data, sleep for some time. */
			Pa_Sleep(500.*transfer.nsamples/fsample);
		}

		if (transfer.refresh1!=0 && transfer.refresh2!=0) {
			/* Both blocks need fresh data, sleep for some time. */
			Pa_Sleep(100);
		}

		if (toc(timeout)>config.timeout) {
			/* There is no new data for too long. */
			printf("ft2audio: timeout\n");
			goto error;
		}

	} /* while true */

	err = Pa_StopStream(stream);
	if (err != paNoError) goto error;
	err = Pa_CloseStream(stream);
	if (err != paNoError) goto error;
	Pa_Terminate();

	free(transfer.block1);
	free(transfer.block2);

	if (ftSocket > 0)
    close_connection(ftSocket);
  else
    ft_stop_buffer_server(ftServer);

	printf("ft2audio finished.\n");
	return err;

error:
	Pa_Terminate();
	if (err != paNoError) {
		fprintf(stderr, "An error occured while using the PortAudio stream\n");
		fprintf(stderr, "Error number: %d\n", err);
		fprintf(stderr, "Error message: %s\n", Pa_GetErrorText(err));
	}
	else {
		fprintf(stderr, "An error occured\n");
	}
	return err;

} /* main */
