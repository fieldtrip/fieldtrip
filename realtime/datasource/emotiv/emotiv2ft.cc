#include <iostream>
#include <fstream>
#include <conio.h>
#include <sstream>
#include <windows.h>
#include <map>
#include <signal.h>
#include <FtBuffer.h>
#include <socketserver.h>

#include "edk.h"
#include "edkErrorCode.h"

#define TOTAL_CHANNELS 22

EE_DataChannel_t targetChannelList[TOTAL_CHANNELS] = {
		ED_COUNTER,
		ED_AF3, ED_F7, ED_F3, ED_FC5, ED_T7, 
		ED_P7, ED_O1, ED_O2, ED_P8, ED_T8, 
		ED_FC6, ED_F4, ED_F8, ED_AF4, ED_GYROX, ED_GYROY, ED_TIMESTAMP, 
		ED_FUNC_ID, ED_FUNC_VALUE, ED_MARKER, ED_SYNC_SIGNAL
};

const char *labels[TOTAL_CHANNELS] = {
	"COUNTER",
	"AF3",
	"F7",
	"F3",
	"FC5",
	"T7",
	"P7",
	"O1",
	"O2",
	"P8",
	"T8", 
	"FC6",
	"F4",
	"F8",
	"AF4",
	"GYROX", 
	"GYROY",
	"TIMESTAMP",
	"FUNC_ID", 
	"FUNC_VALUE", 
	"MARKER", 
	"SYNC_SIGNAL"
};


volatile bool keepRunning = true;

bool writeHeader(int ftSocket, float fSample) {
	FtBufferRequest req;
	FtBufferResponse resp;
	char *chunk_data;
	int N=0,P;
	
	for (int n=0;n<TOTAL_CHANNELS;n++) {
		N+=strlen(labels[n])+1;
	}
	chunk_data = new char[N];
	
	P=0;
	for (int n=0;n<TOTAL_CHANNELS;n++) {
		int Ln = strlen(labels[n])+1;
		memcpy(chunk_data + P, labels[n], Ln);
		P+=Ln;
	}
	

	req.prepPutHeader(TOTAL_CHANNELS, DATATYPE_FLOAT64, fSample);
	req.prepPutHeaderAddChunk(FT_CHUNK_CHANNEL_NAMES, N, chunk_data);
	
	delete[] chunk_data;
	
	int err = clientrequest(ftSocket, req.out(), resp.in());
	if (err || !resp.checkPut()) {
		fprintf(stderr, "Could not write header to FieldTrip buffer\n.");
		return false;
	}
	return true;
}


void abortHandler(int sig) {
	printf("Ctrl-C pressed -- stopping...\n");
	keepRunning = false;
}

					  
int main(int argc, char** argv) {

	EmoEngineEventHandle eEvent;
	EmoStateHandle eState;
	unsigned int userID					= 0;
	//const unsigned short composerPort	= 1726;
	float secs							= 1;
	unsigned int samplingRate = 0;
	bool readyToCollect	= false;
	ft_buffer_server_t *server = NULL;
	int port, ftSocket;
	int sampleCounter = 0;
	char hostname[256];
	FtBufferResponse resp;
	FtSampleBlock sampleBlock(DATATYPE_FLOAT64);
	
	if (argc>1) {
		strncpy(hostname, argv[1], sizeof(hostname));
	} else {
		strcpy(hostname, "localhost");
	}
	
	if (argc>2) {
		port = atoi(argv[2]);
	} else {
		port = 1972;
	}	
	
	if (strcmp(hostname, "-")) {
		ftSocket = open_connection(hostname, port);
		if (ftSocket < 0) {
			fprintf(stderr, "Could not connect to FieldTrip buffer on %s:%i\n", hostname, port);
			return 1;
		}
	} else {
		server = ft_start_buffer_server(port, NULL, NULL, NULL);
		if (server == NULL) {
			fprintf(stderr, "Could not spawn FieldTrip buffer on port %i\n", port);
			return 1;
		}
		ftSocket = 0;
	}
	
	eEvent = EE_EmoEngineEventCreate();
	eState = EE_EmoStateCreate();
		
	if (EE_EngineConnect() != EDK_OK) {
		fprintf(stderr, "Emotiv Engine start up failed.");
		exit(1);
	}
	// TODO: provide alternative connection
	// ...	if (EE_EngineRemoteConnect(input.c_str(), composerPort) != EDK_OK) {
	
	DataHandle hData = EE_DataCreate();
	EE_DataSetBufferSizeInSec(secs);
	
	/* register CTRL-C handler */
	signal(SIGINT, abortHandler);
	printf("Starting to transfer data - press CTRL-C to quit\n");

	while (keepRunning) {
		int state = EE_EngineGetNextEvent(eEvent);

		if (state == EDK_OK) {

			EE_Event_t eventType = EE_EmoEngineEventGetType(eEvent);
			EE_EmoEngineEventGetUserId(eEvent, &userID);

			// Log the EmoState if it has been updated
			if (eventType == EE_UserAdded) {
				if (EE_DataGetSamplingRate(userID, &samplingRate) != EDK_OK) {
					fprintf(stderr, "Cannot retrieve sampling rate\n");
					break;
				}
				
				if (!writeHeader(ftSocket, (float) samplingRate)) return 1;
				
				EE_DataAcquisitionEnable(userID,true);
				readyToCollect = true;
				printf("EDK: User added, will now start acquisition.\n");
			}
		}

		if (readyToCollect) {
			EE_DataUpdateHandle(0, hData);

			unsigned int nSamplesTaken=0;
			EE_DataGetNumberOfSample(hData,&nSamplesTaken);
	
			if (nSamplesTaken != 0) {
				double* data = new double[nSamplesTaken];
				double* dest = (double *) sampleBlock.getMatrix(TOTAL_CHANNELS, nSamplesTaken); 
				for (int i=0;i<TOTAL_CHANNELS;i++) {
					EE_DataGet(hData, targetChannelList[i], data, nSamplesTaken);
						
					for (unsigned int j=0;j<nSamplesTaken;j++) {
						dest[i + j*TOTAL_CHANNELS] = data[j];
					}
				}
				delete[] data;
				int err = clientrequest(ftSocket, sampleBlock.asRequest(), resp.in());
				if (err || !resp.checkPut()) {
					fprintf(stderr, "Could not write samples to FieldTrip buffer\n.");
				}
				sampleCounter += nSamplesTaken;
				printf("Wrote %2i samples (%i total)\n", nSamplesTaken, sampleCounter);
			}
		}
		Sleep(10);
	}
	EE_DataFree(hData);
	if (ftSocket > 0) close_connection(ftSocket);
	if (server != NULL) ft_stop_buffer_server(server);

	EE_EngineDisconnect();
	EE_EmoStateFree(eState);
	EE_EmoEngineEventFree(eEvent);

	return 0;
}


