/*
 * Simple application for playing back an online experiment through the FieldTrip buffer.
 * This is the opposite functionality of the 'recording' application.
 *
 * (C) 2010 Stefan Klanke
 */

#include <signal.h>
#include <pthread.h>
#include <string.h>
#if !defined(WIN32) || defined(COMPILER_MINGW)
#include <sys/time.h>
#endif

#include "buffer.h"
#include "rdadefs.h"

#define MAXLINE 256
#define MAX_PRINT_CHN  300

typedef struct {
	int numSamples;	/* number of samples (or 0 for a write-events-operation)  */
	int numEvents;  /* number of events  (or 0 for a write-samples-operation) */
	long offset;    /* offset into events file and eventBuffer memory blob    */
	int size;       /* number of bytes to transmit for events or samples      */
	double time;    /* time when this needs to be sent, relative to PUT_HDR   */
} WriteOperation;

char directory[MAXLINE];
int ftSocket = -1;
int numWriteOps, allocedWriteOps;
INT64_T totalSamples, totalEvents;
INT64_T sizeSamples, sizeEvents;
int curSamplesFile = 0, numSamplesFiles = 1;
unsigned int bytesPerSample;

WriteOperation *writeOps = NULL;
FILE *fSamples;
headerdef_t *header = NULL;
char *eventBuffer = NULL;
UINT32_T headerSize;

static char usage[] = "Usage: playback <directory> [hostname=localhost [port=1972]]\n";

double getCurrentTime() {
#if defined(WIN32) && !defined(COMPILER_MINGW)
	return timeGetTime() * 0.001;
#else
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + tv.tv_usec*1e-6;
#endif
}

int readHeader(const char *directory) {
	char filename[MAXLINE];
	FILE *f;
	long size;

	snprintf(filename, MAXLINE, "%s/header", directory);
	f = fopen(filename, "rb");
	if (f == NULL) {
		fprintf(stderr, "Can not read file %s\n", filename);
		return -1;
	}

	fseek(f, 0, SEEK_END);
	size = ftell(f);
	fseek(f, 0, SEEK_SET);

	if (header != NULL) free(header);
	header = (headerdef_t *) malloc(size);

	if (header == NULL) {
		fprintf(stderr, "Cannot allocate %li bytes for reading the header.\n", size);
		exit(1);
	}

	fread(header, 1, size, f);
	fclose(f);

	headerSize = size;

	return headerSize;
}


int readTiming(const char *directory, double speedup) {
	FILE *f;
	char filename[MAXLINE];
	long offEvts = 0;

	snprintf(filename, MAXLINE, "%s/timing", directory);
	f = fopen(filename, "r");
	if (f == NULL) {
		fprintf(stderr, "Can not read file %s\n", filename);
		return -1;
	}

	allocedWriteOps = 1000;
	numWriteOps = 0;
	totalSamples = totalEvents = 0;
	if (writeOps != NULL) free(writeOps);
	writeOps = (WriteOperation *) malloc(allocedWriteOps * sizeof(WriteOperation));

	if (writeOps == NULL) {
		fprintf(stderr, "Out of memory\n");
		exit(1);
	}

	while (!feof(f)) {
		char type;
		double time;
		int num, r, i;
		WriteOperation *wop;
		eventdef_t *evdef;

		r = fscanf(f, "%c %i %lf\n", &type, &num, &time);

		if (r!=3) continue;

		if (numWriteOps == allocedWriteOps) {
			void *newmem = realloc(writeOps, (allocedWriteOps + 1000) * sizeof(WriteOperation));
			if (newmem == NULL) {
				fprintf(stderr, "Out of memory\n");
				exit(1);
			}
			writeOps = (WriteOperation *) newmem;
			allocedWriteOps += 1000;
		}

		wop = writeOps + numWriteOps;
		wop->time = time / speedup;
		switch(type) {
			case 'S':
				wop->numSamples = num;
				wop->numEvents = 0;
				wop->size = num * bytesPerSample;
				totalSamples += num;
				break;
			case 'E':
				wop->numSamples = 0;
				wop->numEvents = num;
				wop->offset = offEvts;
				wop->size = 0;
				for (i=0;i<num;i++) {
					int siz;

					evdef = (eventdef_t *) (eventBuffer+offEvts);
					siz = sizeof(eventdef_t) + evdef->bufsize;

					totalEvents++;
					printf("%li. event:  %i bytes @ %li\n", (long) totalEvents, siz, offEvts);

					wop->size += siz;
					offEvts += siz;
				}
				break;
			default:
				printf("Invalid timing definition\n");
				continue;
		}
		numWriteOps++;
	}
	fclose(f);
	return numWriteOps;
}

INT64_T openSamplesFile(const char *directory, int counter) {
	char filename[MAXLINE];
	long size;

	if (counter < 1) {
		snprintf(filename, MAXLINE, "%s/samples", directory);
	} else {
		snprintf(filename, MAXLINE, "%s/samples%i", directory, counter);
	}
	fSamples = fopen(filename, "rb");
	if (fSamples == NULL) {
		if (counter < numSamplesFiles) {
			fprintf(stderr, "Can not read file %s\n", filename);
		}
		return -1;
	}

	fseek(fSamples, 0, SEEK_END);
	size = ftell(fSamples);
	fseek(fSamples, 0, SEEK_SET);
	return size;
}

int readAllEvents(const char *directory) {
	char filename[MAXLINE];
	FILE *f;

	snprintf(filename, MAXLINE, "%s/events", directory);
	f = fopen(filename, "rb");
	if (f == NULL) {
		fprintf(stderr, "Can not read file %s\n", filename);
		return -1;
	}

	fseek(f, 0, SEEK_END);
	sizeEvents = ftell(f);
	fseek(f, 0, SEEK_SET);
	if (sizeEvents > 0) {
		eventBuffer = (char *) malloc(sizeEvents);
		if (eventBuffer == NULL) {
			fprintf(stderr, "Cannot allocate %li bytes for reading events\n", (long) sizeEvents);
		}

		fread(eventBuffer, 1, sizeEvents, f);
	}
	fclose(f);
	return sizeEvents;
}

void run() {
	datadef_t *ddef;
	int maxSize = 0;
	int r;
	int nextSampleOp = numWriteOps;
	messagedef_t reqdef;
	message_t request, *response;
	double T0, t;
	int op;

	/* search for first sample operation */
	for (op = 0; op < numWriteOps; op++) {
		if (writeOps[op].numSamples > 0) {
			nextSampleOp = op;
			break;
		}
	}

	for (op = 0; op < numWriteOps; op++) {
		if (writeOps[op].numSamples > 0 && writeOps[op].size > maxSize) {
			maxSize = writeOps[op].size;
		}
	}

	if (maxSize > 0) {
		ddef = (datadef_t *) malloc(sizeof(datadef_t) + maxSize);
		if (ddef == NULL) {
			printf("Cannot allocate temporary buffer for writing samples\n");
			exit(1);
		}
		ddef->nchans    = header->nchans;
		ddef->data_type = header->data_type;
		ddef->bufsize   = writeOps[nextSampleOp].numSamples * bytesPerSample;
		ddef->nsamples  = writeOps[nextSampleOp].numSamples;
		fread(ddef+1, bytesPerSample, writeOps[nextSampleOp].numSamples, fSamples);
	}

	/* write out header */
	request.def = &reqdef;
	reqdef.version = VERSION;
	reqdef.command = PUT_HDR;
	reqdef.bufsize = headerSize;
	request.buf = header;

	T0 = getCurrentTime();
	printf("Writing header...\n");
	r = clientrequest(ftSocket, &request, &response);

	if (r!=0 || response->def == NULL || response->def->command != PUT_OK) {
		fprintf(stderr, "Error in FieldTrip request\n");
	}
	if (response->buf != NULL) free(response->buf);
	free(response->def);
	free(response);

	op=0;

	for (op=0;op<numWriteOps;op++) {
		t = getCurrentTime() - T0;
		if (writeOps[op].time > t) {
			usleep(1.0e6*(writeOps[op].time -  t));
			t = getCurrentTime() - T0;
		}

		if (writeOps[op].numSamples > 0) {
			reqdef.command = PUT_DAT;
			reqdef.bufsize = sizeof(datadef_t) + writeOps[op].size;
			request.buf = ddef;
			printf("%.3f: Writing %i sample(s)\n", t, writeOps[op].numSamples);
		} else {
			reqdef.command = PUT_EVT;
			reqdef.bufsize = writeOps[op].size;
			request.buf = eventBuffer + writeOps[op].offset;
			printf("%.3f: Writing %i event(s)\n", t, writeOps[op].numEvents);
		}
		r = clientrequest(ftSocket, &request, &response);

		if (r!=0 || response->def == NULL || response->def->command != PUT_OK) {
			fprintf(stderr, "Error in FieldTrip request\n");
		}
		if (response->buf != NULL) free(response->buf);
		free(response->def);
		free(response);

		if (writeOps[op].numSamples > 0) {
			/* pre-load next bunch of samples */
			for (nextSampleOp = op+1; nextSampleOp < numWriteOps; nextSampleOp++) {
				if (writeOps[nextSampleOp].numSamples > 0) break;
			}
			if (nextSampleOp < numWriteOps) {
				INT64_T numSamplesRead;
				char *sampleBuffer;

				ddef->nchans    = header->nchans;
				ddef->data_type = header->data_type;
				ddef->bufsize   = writeOps[nextSampleOp].numSamples * bytesPerSample;
				ddef->nsamples  = writeOps[nextSampleOp].numSamples;

				sampleBuffer = (char *) (ddef+1);

				numSamplesRead = fread(sampleBuffer, bytesPerSample, writeOps[nextSampleOp].numSamples, fSamples);
				if (numSamplesRead < writeOps[nextSampleOp].numSamples) {
					if (curSamplesFile == numSamplesFiles) {
						fprintf(stderr, "Error reading samples!\n");
						break;
					} else {
						int remain = writeOps[nextSampleOp].numSamples - numSamplesRead;

						fclose(fSamples);
						if (openSamplesFile(directory, ++curSamplesFile) < 0) break;

						sampleBuffer += bytesPerSample * numSamplesRead;
						numSamplesRead = fread(sampleBuffer, bytesPerSample, remain, fSamples);
						if (numSamplesRead < remain) {
							fprintf(stderr, "Error reading samples!\n");
							break;
						}
					}
				}
			}
		}
	}
	printf("Done!\n");
}

int main(int argc, char **argv) {
	char hostname[MAXLINE] = "localhost";
	int port, nops;
	double speed = 1.0;

	#if defined(WIN32) && !defined(COMPILER_MINGW)
	timeBeginPeriod(1);
	#endif

	if (argc<2) {
		fputs(usage, stderr);
		exit(1);
	}

	strncpy(directory, argv[1], MAXLINE);
	if (argc>2) {
		strncpy(hostname, argv[2], MAXLINE);
	}
	if (argc>3) {
		port = atoi(argv[3]);
	} else {
		port = 1972;
	}

	if (argc>4) {
		speed = strtod(argv[4], NULL);
		if (speed <= 0.0) {
			fprintf(stderr, "4th argument, if given, must be positive speedup-factor\n");
			exit(1);
		}
	}

	if (readHeader(directory) < 0) {
		exit(1);
	}
	bytesPerSample = header->nchans * wordsize_from_type(header->data_type);

	sizeSamples = openSamplesFile(directory, 0);
	if (sizeSamples < 0) {
		exit(1);
	}
	fclose(fSamples);
	while (1) {
		INT64_T addSize = openSamplesFile(directory, numSamplesFiles);
		if (addSize < 0) break;
		numSamplesFiles++;
		sizeSamples += addSize;
		fclose(fSamples);
	}
	printf("Total size of samples: %li MB\n", (long) (sizeSamples >> 20));
	/* re-open first samples file */
	if (openSamplesFile(directory, 0) < 0) exit(1);

	if (readAllEvents(directory) < 0) exit(1);

	nops = readTiming(directory, speed);
	if (nops < 0) {
		exit(1);
	}
	if (totalSamples == 0 && totalEvents == 0) {
		printf("No samples or events defined\n");
		exit(0);
	} else {
		INT64_T siz = totalSamples * (INT64_T) bytesPerSample;

		printf("Total samples: %li  events: %li\n", (long) totalSamples, (long) totalEvents);

		if (siz > sizeSamples) {
			fputs("Error: 'samples' file(s) too small for given 'timing' definition\n", stderr);
			exit(1);
		}
		if (siz < sizeSamples) {
			printf("Warning: 'samples' file contains %li bytes, but 'timing' definition specifies %li bytes\n", (long) sizeSamples, (long) siz);
		}
	}

	printf("Trying to connect to %s:%i...\n", hostname, port);
	ftSocket = open_connection(hostname, port);
	if (ftSocket < 0) {
		exit(1);
	}

	run();

	close_connection(ftSocket);

	fclose(fSamples);
	free(writeOps);
	free(header);
	if (eventBuffer != NULL) free(eventBuffer);

	#if defined(WIN32) && !defined(COMPILER_MINGW)
	timeEndPeriod(1);
	#endif

	return 0;
}
