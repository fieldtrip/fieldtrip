/*
 * Simple program that fakes a FieldTrip buffer server, and just writes incoming data
 * to NIFTI-files (in case it's MRI data).
 * 
 * Copyright (C) 2010, Stefan Klanke
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "buffer.h"
#include "socketserver.h"
#include <signal.h>
#include "nifti1.h"

#ifdef WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#include <sys/time.h>
#endif


ft_buffer_server_t *S;
volatile int keepRunning = 1;
char baseDirectory[512];
int setCounter = 0;
int sampleCounter = 0;
int hasHeader = 0;

char currentHeader[352];  /* 348 for Nifti-header itself, 4 bytes padding */

void printCorrectNifti() {
	struct nifti_1_header *H = (struct nifti_1_header *) currentHeader;
	
	printf("Dimensions: %i x %i x %i\n", H->dim[1], H->dim[2], H->dim[3]);
	/* apply some corrections to the header, ie. offset=352, magic="n+1", ... */
	H->magic[0] = 'n';
	H->magic[1] = '+';
	H->magic[2] = '1';
	H->magic[3] = 0;
	H->vox_offset = 352.0;
}

void writeProtocol(const ft_chunk_t *chunk) {
	char filename[800];
	FILE *fp;
	
	sprintf(filename, "%s/%03i_siemens_prot.txt", baseDirectory, setCounter);
	fp = fopen(filename, "wb");
	fwrite(chunk->data, 1, chunk->def.size, fp);
	fclose(fp);
}

void writeOneFile(void *data, int size) {
	char filename[800];
	FILE *fp;
	
	sprintf(filename, "%s/%03i_%05i.nii", baseDirectory, setCounter, sampleCounter);
	fp = fopen(filename, "wb");
	fwrite(currentHeader, 1, 352, fp);
	fwrite(data, 1, size, fp);
	fclose(fp);
}

int my_request_handler(const message_t *request, message_t **response, void *user_data) {
	message_t *resp;
	const ft_chunk_t *chunk;
	short retVal;
	
	switch(request->def->command) {
		case PUT_HDR:
			printf("Header received, bufsize = %i\n", request->def->bufsize);
			chunk = find_chunk(request->buf, sizeof(headerdef_t), request->def->bufsize, FT_CHUNK_NIFTI1);
			if (chunk == NULL || chunk->def.size != 348) {
				printf("No NIFTI-1 header present - ignoring\n");
				hasHeader = 0;
				retVal = PUT_ERR;
			} else {
				hasHeader = 1;
				++setCounter;
				sampleCounter = 0;
				memcpy(currentHeader, chunk->data, 348);
				retVal = PUT_OK;
				printCorrectNifti();
				
				chunk = find_chunk(request->buf, sizeof(headerdef_t), request->def->bufsize, FT_CHUNK_SIEMENS_AP);
				if (chunk != NULL) {
					writeProtocol(chunk);
				}
			}
			break;
		case PUT_DAT:
			if (hasHeader) {
				int i, size;
				char *smp;
				
				const datadef_t *ddef = (const datadef_t *) request->buf;
				smp = (char *) (ddef+1); /* points at first byte of first sample/scan */
				size = ddef->nchans * wordsize_from_type(ddef->data_type);
				
				printf("%i scan(s) received\n", ddef->nsamples);
				
				/* ddef->nsamples should be 1, but who knowns what the future brings */
				for (i=0;i<ddef->nsamples;i++) {
					++sampleCounter;
					writeOneFile(smp + i*size, size);
				}
				retVal = PUT_OK;
			} else {
				retVal = PUT_ERR;
			}
			break;
		case PUT_EVT:
			printf("Received event, ignoring...\n");
			retVal = PUT_OK;
			break;
		default:
			retVal = GET_ERR;
			/* all other requests yield an error at the other side */
			return -1;
	}
	/* compose a response so writing end does not get suspicious */
	resp = (message_t *) malloc(sizeof(message_t));
	resp->def = (messagedef_t *) malloc(sizeof(messagedef_t));
	resp->def->version = 1;
	resp->def->command = retVal;
	resp->def->bufsize = 0;
	resp->buf = 0;
	*response = resp;
	return 0;
}

void abortHandler(int sig) {
	keepRunning = 0;
}

int main(int argc, char *argv[]) {
	int port;

	/* padding between nifti header and data */
	currentHeader[348] = 0;
	currentHeader[349] = 0;
	currentHeader[350] = 0;
	currentHeader[351] = 0;

	/* USAGE: 
		save_as_nifti [directory=.] [port=1972]\n");
	*/

	if (argc>1) {
		strncpy(baseDirectory, argv[1], sizeof(baseDirectory));
	} else {
		strcpy(baseDirectory, ".");
	}

	if (argc>2) {
		port = atoi(argv[2]);
	} else {
		port = 1972;
	}
	
	S = ft_start_buffer_server(port, NULL, my_request_handler, NULL);
	if (S==NULL) return 1;
	
	printf("Listening on port %i - press Ctrl-C to stop...\n", port);

	signal(SIGINT, abortHandler);
	while (keepRunning) {
		#ifdef WIN32
		Sleep(1000);
		#else
		sleep(1);
		#endif
	}
	printf("Ctrl-C pressed -- stopping operation...\n");
	ft_stop_buffer_server(S);
	printf("Done.\n");

	return 0;
}

