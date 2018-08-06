/*
 * (C) 2010 Stefan Klanke
 */
#include <signal.h>
#include <string.h>

#include "buffer.h"
#include "socketserver.h"
#include "rdadefs.h"
#include "platform.h"

#define MAXLINE 		256
#define MAX_PRINT_CHN 	300
#define WAIT_INFINITELY	{while (1) {sleep(1);}}
#define ABORT			{exit(1);}

const UINT8_T _rda_guid[16]={
		0x8E,0x45,0x58,0x43,0x96,0xC9,0x86,0x4C,0xAF,0x4A,0x98,0xBB,0xF6,0xC9,0x14,0x50
};

int keepRunning = 1;
int numChannels, ftSocket = -1, goodToSend = 0;
int samplesWritten;

static char usage[] = 
"\n" \
    "Usage: rda2ft rdaHost rdaPort [ftHost ftPort]\n" \
    "\n" \
    "Calling 'rda2ft' with only the first two arguments starts a local buffer on port 1972.\n" \
    "Using '-' for the buffer hostname (ftHost) starts a local buffer on the given port (ftPort).\n" \
    "The port number for the 32 bit float RDA interface is 51244\n" \
    "\n" \
    "Example use:\n" \
    "   rda2ft localhost 51244 - 1972\n" \
    ;

#ifndef PLATFORM_LINUX

size_t strnlen(const char *s, size_t max) {
		size_t n = 0;
		while (n<max && *s++) {
				++n;
		}
		return n;
}

#endif

void abortHandler(int sig) {
		printf("Ctrl-C pressed -- stopping...\n");
		keepRunning = 0;
}

void handleStartPacket(int ftSocket, int size, void *buf) {
		rda_msg_start_t *header = (rda_msg_start_t *) buf;
		double *dRes = (double *) ((char *) buf + sizeof(rda_msg_start_t));
		char *nameStart;
		int sizeNames;

		printf("\nRDA start packet (%i bytes) -- %i channels, Ts=%.1fus\n\n", size, header->nChannels, header->dSamplingInterval);
		numChannels = header->nChannels;

		goodToSend = 0;

		if (numChannels > 0) {
				int i;
				char *name = (char *) buf + sizeof(rda_msg_start_t) + numChannels*sizeof(double);

				nameStart = name;

				printf("Channel  Resolution  Name\n");

				for (i=0;i<numChannels;i++) {
						int maxLen = size - (name - (char *) buf);
						int len = strnlen(name, maxLen);

						if (len>=maxLen) {
								fprintf(stderr, "Invalid START packet received (unterminated channel names) -- ignoring.\n");
								return;
						}
						if (i<MAX_PRINT_CHN) {
								printf("%5i    %10.3f  '%s'\n", i, dRes[i], name);
						}
						name+=len+1;
				}

				sizeNames = name - nameStart;

				if (numChannels > MAX_PRINT_CHN) {
						printf("Suppressed output of further %i channels\n", numChannels - MAX_PRINT_CHN);
				}
		} else {
				fprintf(stderr, "Invalid START packet received -- channel number = %i\n", numChannels);
				return;
		}

		if (ftSocket != -1) {
				message_t request, *response;
				messagedef_t msgdef;
				char *msgbuf;
				headerdef_t *hdr;
				ft_chunk_t *chunk;
				int res;
				int sizeChunks = 2*sizeof(ft_chunkdef_t) + sizeNames + numChannels*sizeof(double);
				int sizeAll    = sizeof(headerdef_t) + sizeChunks;

				msgbuf = (char *) malloc(sizeAll);
				hdr = (headerdef_t *) msgbuf;
				if (hdr == NULL) {
						fprintf(stderr, "Out of memory -- not sending header\n");
						return;
				}

				request.def = &msgdef;
				request.buf = msgbuf;

				msgdef.version = VERSION;
				msgdef.command = PUT_HDR;
				msgdef.bufsize = sizeAll;

				hdr->nchans = numChannels;
				hdr->nsamples  = 0;
				hdr->nevents   = 0;
				hdr->fsample   = 1.0e6/header->dSamplingInterval;
				hdr->data_type = DATATYPE_FLOAT32; 
				hdr->bufsize   = sizeChunks;

				chunk = (ft_chunk_t *) (msgbuf + sizeof(headerdef_t));
				chunk->def.type = FT_CHUNK_CHANNEL_NAMES;
				chunk->def.size = sizeNames;
				memcpy(chunk->data, nameStart, sizeNames);

				chunk = (ft_chunk_t *) (msgbuf + sizeof(headerdef_t) + sizeof(ft_chunkdef_t) + sizeNames);
				chunk->def.type = FT_CHUNK_RESOLUTIONS;
				chunk->def.size = numChannels*sizeof(double);
				memcpy(chunk->data, (char *) buf + sizeof(rda_msg_start_t), numChannels*sizeof(double));

				res = clientrequest(ftSocket, &request, &response);

				if (res == 0) {
						if (response->def->command == PUT_OK) {
								goodToSend = 1;
								samplesWritten = 0;
						}
						free(response->def);
						if (response->buf) free(response->buf);
						free(response);
				}
		}
}

void handleDataPacket(int ftSocket, int size, void *buf) {
		rda_msg_data_t *data = (rda_msg_data_t *) buf;
		int isInt, sizeData;

		message_t request, *response;
		messagedef_t reqdef;
		request.def = &reqdef;

		isInt = (data->hdr.nType == RDA_INT_MSG) ? 1 : 0;
		sizeData = (isInt ? sizeof(INT16_T) : sizeof(float)) * data->nPoints * numChannels;

		printf("RDA data block %4i (%i bytes):  %4i samples, %2i markers\n", data->nBlock, data->hdr.nSize, data->nPoints, data->nMarkers);	

		if (data->nMarkers > 0) {
				char *evbuf;

				/* offset of markers into RDA packet */
				int offset = sizeof(rda_msg_data_t) + sizeData;
				/* maximum length of all the "type" strings combined (actually a bit more because of trailing zeros) */
				int maxSizeTypes = size - offset - sizeof(rda_marker_t)*data->nMarkers;
				/* offset + size of events in evbuf */
				int evsiz = 0;

				evbuf = (char *) malloc(data->nMarkers*(sizeof(eventdef_t) + sizeof(INT32_T)) + maxSizeTypes);
				if (evbuf == NULL) {
						fprintf(stderr, "Out of memory\n");
						return;
				}

				while (offset + sizeof(rda_marker_t) <= size) {
						eventdef_t *evdef;
						rda_marker_t *marker = (rda_marker_t *) ((char *) buf + offset);
						char *markerType     = (char *) buf + offset + sizeof(rda_marker_t);
						int maxLen           = size - (markerType - (char *) buf);
						int typeLen          = strnlen(markerType, maxLen);
						char *markerValue    = (char *) buf + offset + sizeof(rda_marker_t) + typeLen + 1;
						int valueLen         = strnlen(markerValue, maxLen);

						/* The type and the value are both strings that are represented in each marker like this "Stimulus\0S  1\0"
						 * which codes for a Marker type "Stimulus" and the value "S  1". 
						 */ 

						printf("Marker: Pos=%i  Length=%i  Size=%i  Channel=%i  maxLen=%i  markerType=%s  typeLen=%d  markerValue=%s  valueLen=%d\n",
										marker->nPosition, marker->nPoints, marker->nSize, marker->nChannel, maxLen, markerType, typeLen, markerValue, valueLen);

						offset += marker->nSize;

						evdef = (eventdef_t *) evbuf + evsiz;
						evdef->type_type = DATATYPE_CHAR;
						evdef->type_numel = typeLen;
						evdef->value_type = DATATYPE_CHAR;
						evdef->value_numel = valueLen;
						evdef->sample = samplesWritten + marker->nPosition;
						evdef->offset = 0;
						evdef->duration = marker->nPoints;
						evdef->bufsize = typeLen + valueLen;

						memcpy(evbuf + evsiz + sizeof(eventdef_t)        , markerType, typeLen);
						memcpy(evbuf + evsiz + sizeof(eventdef_t)+typeLen, markerValue, valueLen);

						evsiz += evdef->bufsize + sizeof(eventdef_t);
				}

				reqdef.version = VERSION;
				reqdef.command = PUT_EVT;
				reqdef.bufsize = evsiz;
				request.buf = evbuf;
				if (clientrequest(ftSocket, &request, &response)) {
						fprintf(stderr, "Error when writing events to buffer.\n");
				} else {
						if (response->def->command != PUT_OK) {
								fprintf(stderr, "Buffer server returned an error (writing events).\n");
						}
						free(response->def);
						if (response->buf) free(response->buf);
						free(response);
				}
				free(evbuf);
		}	

		if (goodToSend && data->nPoints > 0) {
				/* Write samples to FieldTrip buffer */
				datadef_t *ddef = (datadef_t *) malloc(sizeof(datadef_t) + sizeData);
				if (ddef == NULL) {
						fprintf(stderr, "Out of memory!\n");
						return;
				}
				ddef->nsamples = data->nPoints;
				ddef->nchans = numChannels;
				ddef->data_type = isInt ? DATATYPE_INT16 : DATATYPE_FLOAT32;
				ddef->bufsize = sizeData;
				/* ddef+1 points at first byte after datadef, data+1 points at first byte after RDA data */
				memcpy(ddef + 1, data + 1, sizeData);

				reqdef.version = VERSION;
				reqdef.command = PUT_DAT;
				reqdef.bufsize = sizeof(datadef_t) + sizeData;
				request.buf = ddef;
				if (clientrequest(ftSocket, &request, &response)) {
						fprintf(stderr, "Error when writing samples to buffer.\n");
				} else {
						if (response->def->command != PUT_OK) {
								fprintf(stderr, "Buffer server returned an error (writing samples).\n");
						} else {
								samplesWritten += data->nPoints;
						}
						free(response->def);
						if (response->buf) free(response->buf);
						free(response);
				}
				free(ddef);
		}	
}

int main(int argc, char **argv) {
		host_t ftHost, rdaHost;
		int rdaSocket;
		ft_buffer_server_t *S;

		if (sizeof(rda_msg_hdr_t)!=24) {
				fprintf(stderr, "Compiled-in datatypes do not match RDA protocol\n");
				exit(1);
		}


		if (argc!=5 && argc!=3) {
				fputs(usage, stderr);
				exit(1);
		}

		strncpy(rdaHost.name, argv[1], sizeof(rdaHost.name));
		rdaHost.name[sizeof(rdaHost.name)-1]=0;
		rdaHost.port = atoi(argv[2]);

		if (argc==5) {
				strncpy(ftHost.name, argv[3], sizeof(ftHost.name));
				ftHost.name[sizeof(ftHost.name)-1]=0;
				ftHost.port = atoi(argv[4]);
		} else {
				strcpy(ftHost.name, "-");
				ftHost.port = 1972;
		}

		/* this isn't a FieldTrip connection, but we can use the same function */
		rdaSocket = open_connection(rdaHost.name, rdaHost.port);
		if (rdaSocket < 0) {
				fprintf(stderr, "Could not connect to RDA server at %s:%i\n", rdaHost.name, rdaHost.port);
				return 1;
		}

		/* Spawn tcpserver or connect to remote buffer */
		if (strcmp(ftHost.name, "-") == 0) {
				S = ft_start_buffer_server(ftHost.port, NULL, NULL, NULL);
				if (S==NULL) {
						fprintf(stderr, "Could not start up a FieldTrip buffer serving at port %i\n", ftHost.port);
						return 1;
				}
				ftSocket = 0; 
				printf("Streaming from %s:%i to local buffer on port %i\n", rdaHost.name, rdaHost.port, ftHost.port);
		} else {
				ftSocket = open_connection(ftHost.name, ftHost.port);

				if (ftSocket < 0) {
						fprintf(stderr, "Could not connect to FieldTrip buffer at %s:%i\n", ftHost.name, ftHost.port);
						return 1;
				}
				printf("Streaming from %s:%i to remote buffer at %s:%i\n", rdaHost.name, rdaHost.port, ftHost.name, ftHost.port);
		}		

		/* register CTRL-C handler */
		signal(SIGINT, abortHandler);

		printf("Starting to listen - press CTRL-C to quit\n");

		while (keepRunning) {
				rda_msg_hdr_t header;
				char *buf;
				int n,s;

				n = bufread(rdaSocket, &header, sizeof(header));

				if (n!=sizeof(header)) {
						fprintf(stderr, "Error while reading packet header from the RDA server\n");
						break;
				}

				if (memcmp(_rda_guid, header.guid, 16)!=0) {
						int i;
						fprintf(stderr, "Incorrect GUID in received packet:\n");
						for (i=0;i<16;i++) {
								fprintf(stderr, "0x%02x   should be 0x%02x\n", header.guid[i], _rda_guid[i]);
						}
						break;
				}

				buf = (char *) malloc(header.nSize);
				if (buf == NULL) {
						fprintf(stderr, "Out of memory\n");
						break;
				}

				memcpy(buf, &header, sizeof(header));
				s = header.nSize - sizeof(header);
				n = bufread(rdaSocket, buf + sizeof(header), s);
				if (s!=n) {
						fprintf(stderr, "Error while reading packet remainder from the RDA server\n");
						free(buf);
						break;
				}

				switch(header.nType) {
						case RDA_START_MSG:
								handleStartPacket(ftSocket, header.nSize, buf);
								break;
						case RDA_INT_MSG:
						case RDA_FLOAT_MSG:			
								handleDataPacket(ftSocket, header.nSize, buf);
								break;
						case RDA_STOP_MSG:
								printf("\nRemote Data Acquisition stopped\n\n");
								break;
				}

				free(buf);
		}

		close_connection(rdaSocket);

		if (ftSocket > 0) {
				close_connection(ftSocket);
		} else {
				ft_stop_buffer_server(S);
		}

		return 0;
}
