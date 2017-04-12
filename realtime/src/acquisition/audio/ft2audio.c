/*
 * Acquisition tool to stream audio signals from a FieldTrip buffer,
 * to the sound card.This is  based on PortAudio.
 * 	
 * (C) 2017, Robert Oostenveld
*/

#include <stdio.h>
#include <math.h>

#include "portaudio.h"
#include "message.h"
#include "buffer.h"
#include "interface.h"

#define NUM_SECONDS   (2)
#define TRUE          (1)
#define CALIB         (1000.)
#define MAX_RETRY     (500)
#define SAMPLE_RATE   44100
#define BLOCK_SIZE     (SAMPLE_RATE/2)
#define MIN(x,y) (x<y ? x : y)
#define MAX(x,y) (x>y ? x : y)

typedef struct
{
	float *block1;
	float *block2;
	int refresh1;
	int refresh2;
	int current;
	int sample;
}
userData_t;

/* This routine will be called by the PortAudio engine when audio is needed.
 ** It may called at interrupt level on some machines so don't do anything
 ** that could mess up the system like calling malloc() or free().
 */
static int paWriteCallback( const void *inputBuffer,
		void *outputBuffer,
		unsigned long framesPerBuffer,
		const PaStreamCallbackTimeInfo* timeInfo,
		PaStreamCallbackFlags statusFlags,
		void *userData )
{
	/* Cast data passed through stream to our structure. */
	userData_t *userdata = (userData_t*)userData;
	float *output = (float*)outputBuffer;
	(void) inputBuffer; /* Prevent unused variable warning. */
	unsigned int i, dropped = 0;

	for( i=0; i<framesPerBuffer; i++ ) {
		if (userdata->current==1) {
			*output = userdata->block1[userdata->sample];
		}
		else if (userdata->current==2) {
			*output = userdata->block2[userdata->sample];
		}
		else {
			*output = 0;
			dropped++;
		}
		output++;
		userdata->sample++;
		if (userdata->sample==BLOCK_SIZE) {
			/* switch to the start of the other block */
			if (userdata->current==1) {
				userdata->refresh1 = 1;
				userdata->current = (userdata->refresh2 ? 0 : 2);
				userdata->sample = 0;
			}
			else if (userdata->current==2) {
				userdata->refresh2 = 1;
				userdata->current = (userdata->refresh1 ? 0 : 1);
				userdata->sample = 0;
			}
		}
	}
	if (dropped)
		printf("ft2audio: dropped = %d\n", dropped);

	return 0;
}

/*******************************************************************/
static userData_t data;

int main(void)
{
	PaStream *stream;
	PaError err;
	int i, server = 0, status, retry = MAX_RETRY+1;
	float phase = 0;
	/* header information */
	UINT32_T datatype;
	float fsample;
	int nchans, nsamples, nevents;
	/* data from the fieldtrip buffer */
	int begsample, endsample;
	void *rawdata;

	/* Initialize library before making any other calls. */
	err = Pa_Initialize();
	if( err != paNoError ) goto error;

	/* Open an audio I/O stream. */
	err = Pa_OpenDefaultStream( &stream, 0, 1, paFloat32, SAMPLE_RATE, SAMPLE_RATE/100, paWriteCallback, &data );
	if( err != paNoError ) goto error;

	err = Pa_StartStream( stream );
	if( err != paNoError ) goto error;

	server = open_connection("localhost", 1972);
	printf("ft2audio: server = %d\n", server);

	/* Initialize our data for use by callback. */
	data.block1 = malloc(BLOCK_SIZE*sizeof(float));
	data.block2 = malloc(BLOCK_SIZE*sizeof(float));
	DIE_BAD_MALLOC(data.block1);
	DIE_BAD_MALLOC(data.block2);
	bzero(data.block1, BLOCK_SIZE*sizeof(float));
	bzero(data.block2, BLOCK_SIZE*sizeof(float));
	data.refresh1 = 1; /* needs new data */
	data.refresh2 = 1; /* needs new data */
	data.current  = 0; /* neither one contains data */
	data.sample   = 0;

	status = read_header(server, &datatype, &nchans, &fsample, &nsamples, &nevents);
	if (status) goto error;

	printf("ft2audio: nchans   = %d\n", nchans);
	printf("ft2audio: fsample  = %f\n", fsample);
	printf("ft2audio: nsamples = %d\n", nsamples);
	printf("ft2audio: nevents  = %d\n", nevents);

	if (nchans!=1) {
		printf("ft2audio: unsupported nchans = %d\n", nchans);
		goto error;
	}
	if (fsample!=SAMPLE_RATE) {
		printf("ft2audio: unsupported fsample = %f\n", fsample);
		goto error;
	}

	/* junp to the end of the available data */
	begsample = nsamples-BLOCK_SIZE;
	endsample = nsamples-1;

	switch (datatype) {
		case DATATYPE_INT16:
			rawdata = malloc(nchans*BLOCK_SIZE*WORDSIZE_INT16);
			break;
		case DATATYPE_INT32:
			rawdata = malloc(nchans*BLOCK_SIZE*WORDSIZE_INT32);
			break;
		case DATATYPE_FLOAT32:
			rawdata = malloc(nchans*BLOCK_SIZE*WORDSIZE_FLOAT32);
			break;
		default:
			printf("ft2audio: unsupported datatype = %d\n", datatype);
			goto error;
	}

	while (1) {
		if (data.refresh1) {
			retry++;
			status = read_data(server, begsample, endsample, rawdata);
			/* status ==0 when new data or !=0 when data not yet ready */
			if (status==0) {
				printf("ft2audio: refresh1\n");
				begsample += BLOCK_SIZE;
				endsample += BLOCK_SIZE;
				for (i=0; i<BLOCK_SIZE; i++) {
					switch (datatype) {
						case DATATYPE_INT16:
							data.block1[i] = CALIB*((INT16_T *)rawdata)[i];
							break;
						case DATATYPE_INT32:
							data.block1[i] = CALIB*((INT32_T *)rawdata)[i];
							break;
						case DATATYPE_FLOAT32:
							data.block1[i] = CALIB*((FLOAT32_T *)rawdata)[i];
							break;
						default:
							printf("ft2audio: unsupported datatype = %d\n", datatype);
							goto error;
					}
				}
				data.refresh1 = 0;
				retry = 0;
				if (data.current==0)
					data.current = 1; /* this is needed to get started */
			}
		}

		if (data.refresh2) {
			retry++;
			status = read_data(server, begsample, endsample, rawdata);
			/* status ==0 when new data or !=0 when data not yet ready */
			if (status==0) {
				printf("ft2audio: refresh2\n");
				begsample += BLOCK_SIZE;
				endsample += BLOCK_SIZE;
				for (i=0; i<BLOCK_SIZE; i++) {
					switch (datatype) {
						case DATATYPE_INT16:
							data.block2[i] = CALIB*((INT16_T *)rawdata)[i];
							break;
						case DATATYPE_INT32:
							data.block2[i] = CALIB*((INT32_T *)rawdata)[i];
							break;
						case DATATYPE_FLOAT32:
							data.block2[i] = CALIB*((FLOAT32_T *)rawdata)[i];
							break;
					}
				}
				data.refresh2 = 0;
				retry = 0;
			}
		}

		if (retry>MAX_RETRY) {
			status = read_header(server, &datatype, &nchans, &fsample, &nsamples, &nevents);
			if (status) goto error;
			printf("ft2audio: nsamples = %d, waiting for %d\n", nsamples, endsample);
			retry = 0;
		}

		if (!data.refresh1 && !data.refresh2) {
			/* sleep for a short time in miliseconds */
			Pa_Sleep(1000*(BLOCK_SIZE/10)/SAMPLE_RATE);
		}

	} /* while */

	err = Pa_StopStream( stream );
	if( err != paNoError ) goto error;
	err = Pa_CloseStream( stream );
	if( err != paNoError ) goto error;
	Pa_Terminate();

	free(data.block1);
	free(data.block2);
	close_connection(server);

	printf("ft2audio finished.\n");
	return err;
error:
	Pa_Terminate();
	fprintf( stderr, "An error occured while using the portaudio stream\n" );
	fprintf( stderr, "Error number: %d\n", err );
	fprintf( stderr, "Error message: %s\n", Pa_GetErrorText( err ) );
	return err;
}
