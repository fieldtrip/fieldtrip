/*
 * Command-line application to stream signals from a FieldTrip buffer
 * to the sound card. This is based on PortAudio. The FieldTrip buffer
 * should contain one channel at the appropriate sampling frequency.
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

#define FS8K   			 	(8000.)
#define FS16K   			(16000.)
#define FS32K   			(32000.)
#define FS48K   			(48000.)
#define FS96K   			(96000.)
#define FS192K   			(192000.)
#define FS384K   			(384000.)

#define FS11K   			(11025.)
#define FS22K   			(22050.)
#define FS44K   			(44100.)
#define FS88K   			(88200.)
#define FS176K   			(176400.)
#define FS358K   			(352800.)

#define MONITOR       (1000)  /* in milliseconds */
#define LASTREAD      (100)
#define TIMEOUT       (15000) /* in milliseconds */
#define SLIP          (3000)  /* in milliseconds */
#define BLOCKSIZE_DIV (1)
#define CALLBACK_DIV  (100)
#define TRUE          (1)
#define CALIB         (1.)
#define MIN(x,y) (x<y ? x : y)
#define MAX(x,y) (x>y ? x : y)

typedef struct
{
	float *block1;
	float *block2;
	char refresh1;
	char refresh2;
	unsigned int current;
	unsigned int sample;
	unsigned int blocksize;
	unsigned int nread;
	unsigned int nwrite;
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
	/* Cast data passed through stream to our structure. */
	callbackData_t *data = ptr;
	float *output = (float*)outputBuffer;
	(void) inputBuffer; /* Prevent unused variable warning. */

  /* this is to get started */
  if (data->current==0 && data->refresh1==0)
	  data->current = 1; /* switch to the first block */
	else if (data->current==0 && data->refresh2==0)
		data->current = 2; /* switch to the second block */

	while (framesPerBuffer>0) {
		framesPerBuffer--;

		if (data->current==1) {
			*output = data->block1[data->sample];
			data->sample++;
			data->nwrite++;
		}
		else if (data->current==2) {
			*output = data->block2[data->sample];
			data->sample++;
			data->nwrite++;
		}
		else {
			*output = 0;
			data->skipped++;
		}
		output++;

		if (data->sample==data->blocksize) {
			/* try to switch to the start of the other block */
			if (data->current==1) {
				data->refresh1 = 1;
				data->current = (data->refresh2 ? 0 : 2);
				data->sample = 0;
			}
			else if (data->current==2) {
				data->refresh2 = 1;
				data->current = (data->refresh1 ? 0 : 1);
				data->sample = 0;
			}
			else {
				data->current = 0;
				data->sample = 0;
			}
		} /* if switch to other block */
	} /* while framesPerBuffer */

	return 0;
} /* paWriteCallback */

/*******************************************************************/
void tic(struct timeval *now) {
	/* set a timer */
  gettimeofday(now, NULL);
  return;
}

/*******************************************************************/
int toc(struct timeval previous) {
	 /* return elapsed time in milliseconds */
   struct timeval now;
   gettimeofday(&now, NULL);
   return (now.tv_sec - previous.tv_sec)*1000 + (now.tv_usec - previous.tv_usec)/1000;
}

/*******************************************************************/
int main(void)
{
	PaStream *stream;
	PaError err;
	callbackData_t data;
	struct timeval timeout, monitor, lastread;
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
	printf("ft2audio: nevents  = %d\n", nevents);

	if (nchans!=1) {
		printf("ft2audio: unsupported nchans = %d\n", nchans);
		goto error;
	}

	if (fsample!=FS8K && fsample!=FS16K && fsample!=FS32K && fsample!=FS48K && fsample!=FS96K && fsample!=FS192K && fsample!=FS384K && fsample!=FS11K && fsample!=FS22K && fsample!=FS44K && fsample!=FS88K && fsample!=FS176K && fsample!=FS358K) {
		printf("ft2audio: unsupported fsample = %.0f\n", fsample);
		goto error;
	}

	/* Initialize our data for use by the callback. */
	data.refresh1 = 1; /* needs new data */
	data.refresh2 = 1; /* needs new data */
	data.current  = 0; /* neither one contains data */
	data.sample   = 0;
	data.nread    = 0; /* how many samples are read from the data stream */
	data.nwrite   = 0; /* how many samples are written to the audio stream */
	data.skipped  = 0; /* how often was a sample skipped in the audio stream */
	data.blocksize = fsample/BLOCKSIZE_DIV;
	data.block1 = malloc(data.blocksize*sizeof(float));
	data.block2 = malloc(data.blocksize*sizeof(float));
	DIE_BAD_MALLOC(data.block1);
	DIE_BAD_MALLOC(data.block2);
	bzero(data.block1, data.blocksize*sizeof(float));
	bzero(data.block2, data.blocksize*sizeof(float));

	switch (datatype) {
		case DATATYPE_INT16:
			bufdata = malloc(nchans*data.blocksize*WORDSIZE_INT16);
			break;
		case DATATYPE_INT32:
			bufdata = malloc(nchans*data.blocksize*WORDSIZE_INT32);
			break;
		case DATATYPE_INT64:
			bufdata = malloc(nchans*data.blocksize*WORDSIZE_INT64);
			break;
		case DATATYPE_FLOAT32:
			bufdata = malloc(nchans*data.blocksize*WORDSIZE_FLOAT32);
			break;
		case DATATYPE_FLOAT64:
			bufdata = malloc(nchans*data.blocksize*WORDSIZE_FLOAT64);
			break;
		default:
			printf("ft2audio: unsupported datatype = %d\n", datatype);
			goto error;
	}

	/* Initialize library before making any other calls. */
	err = Pa_Initialize();
	if(err != paNoError) goto error;

	printf("ft2audio: initialized PortAudio\n");

	/* Open an audio I/O stream. */
	err = Pa_OpenDefaultStream(&stream, 0, 1, paFloat32, fsample, fsample/CALLBACK_DIV, paWriteCallback, &data);
	if(err != paNoError) goto error;

	printf("ft2audio: opened PortAudio stream\n");

	err = Pa_StartStream(stream);
	if(err != paNoError) goto error;

	printf("ft2audio: started PortAudio stream\n");

	/* jump to the end of the available data */
	/* the case of no data will be dealt with further down */
	begsample = nsamples-data.blocksize;
	endsample = nsamples-1;

	/* initialize all timers */
	tic(&timeout);
	tic(&monitor);
	tic(&lastread);

	while (1) {

		/* monitor the lead or lag between the input and output */
		if (toc(monitor)>MONITOR) {
			tic(&monitor); /* reset the monitor timer */

			status = read_header(server, &datatype, &nchans, &fsample, &nsamples, &nevents);
			if (status) goto error;
			tic(&lastread); /* reset the read timer */

			long data_lead = (long)endsample-(nsamples-1);
			long audio_lag = (long)data.nread-data.nwrite-data.skipped;
			char jump_to_last = 0, jump_to_next = 0;

			printf("ft2audio: data lead = %ld, audio lag = %ld, skipped = %d, read = %d\n", data_lead, audio_lag, data.skipped, data.nread);

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
				begsample = nsamples-data.blocksize;
				endsample = nsamples-1;
				/* reset the audio stream counters */
				data.nread = 0;
				data.nwrite = 0;
				data.skipped = 0;
			}
			else if (jump_to_next) {
				/* jump to the first non-available block */
				begsample = nsamples;
				endsample = nsamples+data.blocksize-1;
				/* reset the audio stream counters */
				data.nread = 0;
				data.nwrite = 0;
				data.skipped = 0;
			}

		} /* if monitor */

		if (data.refresh1==0 && data.refresh2==0) {
			/* both blocks have fresh data, sleep for some time */
			Pa_Sleep(1000.*data.blocksize/fsample);
		}

		if (data.refresh1 && toc(lastread)>LASTREAD) {
			status = read_data(server, begsample, endsample, bufdata);
			tic(&lastread); /* reset the read timer */
			/* status ==0 when new data or !=0 when data not (yet) available */
			if (status==0) {
				retry = 0;
				begsample += data.blocksize;
				endsample += data.blocksize;
				for (unsigned int i=0; i<data.blocksize; i++) {
					switch (datatype) {
						case DATATYPE_INT16:
							data.block1[i] = CALIB*((INT16_T *)bufdata)[i];
							break;
						case DATATYPE_INT32:
							data.block1[i] = CALIB*((INT32_T *)bufdata)[i];
							break;
						case DATATYPE_INT64:
							data.block1[i] = CALIB*((INT64_T *)bufdata)[i];
							break;
						case DATATYPE_FLOAT32:
							data.block1[i] = CALIB*((FLOAT32_T *)bufdata)[i];
							break;
						case DATATYPE_FLOAT64:
							data.block1[i] = CALIB*((FLOAT64_T *)bufdata)[i];
							break;
						default:
							printf("ft2audio: unsupported datatype = %d\n", datatype);
							goto error;
					}
				}
				data.refresh1 = 0;
				data.nread += data.blocksize;
				// printf("ft2audio: refreshed block 1\n");
			}
			else {
				retry++;
			}
		}

		if (data.refresh2 && toc(lastread)>LASTREAD) {
			status = read_data(server, begsample, endsample, bufdata);
			tic(&lastread); /* reset the read timer */
			/* status ==0 when new data or !=0 when data not (yet) available */
			if (status==0) {
				retry = 0;
				begsample += data.blocksize;
				endsample += data.blocksize;
				for (unsigned int i=0; i<data.blocksize; i++) {
					switch (datatype) {
						case DATATYPE_INT16:
							data.block2[i] = CALIB*((INT16_T *)bufdata)[i];
							break;
						case DATATYPE_INT32:
							data.block2[i] = CALIB*((INT32_T *)bufdata)[i];
							break;
						case DATATYPE_INT64:
							data.block2[i] = CALIB*((INT64_T *)bufdata)[i];
							break;
						case DATATYPE_FLOAT32:
							data.block2[i] = CALIB*((FLOAT32_T *)bufdata)[i];
							break;
						case DATATYPE_FLOAT64:
							data.block2[i] = CALIB*((FLOAT64_T *)bufdata)[i];
							break;
						default:
							printf("ft2audio: unsupported datatype = %d\n", datatype);
							goto error;
					}
				}
				data.refresh2 = 0;
				data.nread += data.blocksize;
				// printf("ft2audio: refreshed block 2\n");
			}
			else {
				retry++;
			}
		}

		if (!retry) {
			/* reading was fine, reset the timer */
			tic(&timeout);
		}
		else {
			if (retry && toc(timeout)>TIMEOUT) {
				printf("ft2audio: timeout\n");
				goto error;
			}
			else {
				/* wait a bit for the next attempt */
				Pa_Sleep(100);
			}
		}

	} /* while */

	err = Pa_StopStream(stream);
	if(err != paNoError) goto error;
	err = Pa_CloseStream(stream);
	if(err != paNoError) goto error;
	Pa_Terminate();

	free(data.block1);
	free(data.block2);
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
		fprintf(stderr, "An unknown error occured\n");
	}
	return err;

} /* main */
