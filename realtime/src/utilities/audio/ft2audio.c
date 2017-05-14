/*
 * Command-line application to stream signals from a FieldTrip buffer
 * to the sound card. This is based on PortAudio. The FieldTrip buffer
 * should contain one or two channels at an appropriate sampling frequency.
 *
 * (C) 2017, Robert Oostenveld
*/

#include <stdio.h>
#include <math.h>
#include <sys/time.h>

#include "portaudio.h"
#include "message.h"
#include "buffer.h"
#include "interface.h"

#define MONITOR       (2000)  /* in milliseconds */
#define TIMEOUT       (60000) /* in milliseconds */
#define SLIP          (3000)  /* in milliseconds */
#define BLOCK_DIV 		(1)     /* sampling frequency division factor */
#define CALLBACK_DIV  (10)    /* sampling frequency division factor */
#define CALIB         (1.)    /* calibration factor */
#define TRUE          (1)
#define MIN(x,y) (x<y ? x : y)
#define MAX(x,y) (x>y ? x : y)

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


/* This routine will be called by the PortAudio engine when audio is needed.
 ** It may called at interrupt level on some machines so don't do anything
 ** that could mess up the system like calling malloc() or free().
 */
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
int main(void)
{
	PaStream *stream;
	PaError err;
	callbackData_t transfer;
	struct timeval timeout, monitor;
	int server, status, retry;
	UINT32_T datatype;
	float fsample;
	unsigned int nchans, nsamples, nevents, begsample, endsample;
	void *bufdata;

	server = open_connection("localhost", 1972);
	printf("ft2audio: connected to FieldTrip buffer server = %d\n", server);

	status = read_header(server, &datatype, &nchans, &fsample, &nsamples, &nevents);
	if (status) goto error;

	printf("ft2audio: nchans   = %d\n", nchans);
	printf("ft2audio: fsample  = %.0f\n", fsample);
	printf("ft2audio: nsamples = %d\n", nsamples);

	/* Initialize the structure that is shared with the callback. */
	transfer.refresh1 = 1; /* needs new data */
	transfer.refresh2 = 1; /* needs new data */
	transfer.current  = 0; /* neither one contains data */
	transfer.datatype = datatype;
	transfer.nchans   = nchans;
	transfer.fsample  = fsample;
	transfer.nsamples = fsample/BLOCK_DIV;
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

	/* Initialize library before making any other calls. */
	err = Pa_Initialize();
	if (err != paNoError) goto error;

	printf("ft2audio: initialized PortAudio\n");

	/* Open an audio I/O stream. */
	err = Pa_OpenDefaultStream(&stream, 0, nchans, paFloat32, fsample, fsample/CALLBACK_DIV, paWriteCallback, &transfer);
	if (err != paNoError) goto error;

	printf("ft2audio: opened PortAudio stream\n");

	err = Pa_StartStream(stream);
	if (err != paNoError) goto error;

	printf("ft2audio: started PortAudio stream\n");

	/* jump to the end of the available data */
	/* the case of no data will be dealt with further down */
	begsample = nsamples-transfer.nsamples;
	endsample = nsamples-1;

	/* initialize all timers */
	tic(&timeout);
	tic(&monitor);

	while (1) {

		/* monitor the lead or lag between the input and output */
		if (toc(monitor)>MONITOR) {
			tic(&monitor); /* reset the monitor timer */

			status = read_header(server, &datatype, &nchans, &fsample, &nsamples, &nevents);
			if (status) goto error;

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

			printf("ft2audio: data lead = %ld, audio lag = %ld, skipped = %d, read = %d\n", data_lead, audio_lag, transfer.skipped, transfer.nput);

			if (data_lead > +fsample*SLIP/1000) {
				/* reading data too far in the future, data is coming in too slow */
				printf("ft2audio: data lead too large, reverting to the end\n");
				jump_to_next = 1;
			}
			else if (data_lead < -fsample*SLIP/1000) {
				/* reading data too far in the past, playback is too slow */
				printf("ft2audio: data lag too large, skipping to the end\n");
				jump_to_last = 1;
			}
			else if (audio_lag > +fsample*SLIP/1000) {
				/* writing data slower than it is arriving, playback is too slow */
				printf("ft2audio: audio lag too large, skipping to the end\n");
				jump_to_last = 1;
			}
			else if (audio_lag < -fsample*SLIP/1000) {
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

		if (transfer.refresh1==0 && transfer.refresh2==0) {
			/* both blocks have fresh data, sleep for some time */
			Pa_Sleep(500.*transfer.nsamples/fsample);
		}

		if (transfer.refresh1) {
			status = wait_data(server, endsample, 0, 1200.*transfer.nsamples/fsample);
			if (status!=0)
			  break;
			status = read_data(server, begsample, endsample, bufdata);
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
						transfer.block1[i] *= CALIB;
					}
				}
				transfer.refresh1 = 0;
				transfer.nput += transfer.nsamples;
			}
			else {
				printf("ft2audio: reading block1 failed\n");
			}
		}

		if (transfer.refresh2) {
			status = wait_data(server, endsample, 0, 1200.*transfer.nsamples/fsample);
			if (status!=0)
			  break;
			status = read_data(server, begsample, endsample, bufdata);
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
						transfer.block2[i] *= CALIB;
					}
				}
				transfer.refresh2 = 0;
				transfer.nput += transfer.nsamples;
			}
			else {
				printf("ft2audio: reading block2 failed\n");
			}
		}

		if (toc(timeout)>TIMEOUT) {
			printf("ft2audio: timeout\n");
			goto error;
		}

	} /* while */

	err = Pa_StopStream(stream);
	if (err != paNoError) goto error;
	err = Pa_CloseStream(stream);
	if (err != paNoError) goto error;
	Pa_Terminate();

	free(transfer.block1);
	free(transfer.block2);
	close_connection(server);

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
