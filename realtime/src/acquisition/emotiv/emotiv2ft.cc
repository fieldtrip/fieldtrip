#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

#ifdef PLATFORM_WIN32
    #include <windows.h>
    void dosleep(unsigned milliseconds) {
        Sleep(milliseconds);
    }
#else
    #include <unistd.h>
    void dosleep(unsigned milliseconds) {
        usleep(milliseconds * 1000); // takes microseconds
    }
#endif

#include <OnlineDataManager.h>
#include <ConsoleInput.h>
#include <StringServer.h>

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
/*
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
*/

EmoEngineEventHandle eEvent;
EmoStateHandle eState;	
unsigned int userID	= 0;

int port, ctrlPort;
char hostname[256];
StringServer ctrlServ;
ConsoleInput conIn;

void acquisition(const char *configFile, unsigned int fSample) {
	int sampleCounter = 0;
	
	OnlineDataManager<double, double> ODM(0, TOTAL_CHANNELS, (float) fSample);
	
	if (ODM.configureFromFile(configFile) != 0) {
		fprintf(stderr, "Configuration %s file is invalid\n", configFile);
		return;
	} else {
		printf("Streaming %i out of %i channels\n", ODM.getSignalConfiguration().getStreamingSelection().getSize(), TOTAL_CHANNELS);
	}
	if (!strcmp(hostname, "-")) {
		if (!ODM.useOwnServer(port)) {
			fprintf(stderr, "Could not spawn buffer server on port %d.\n",port);
			return;
		}
	} else {
		if (!ODM.connectToServer(hostname, port)) {
			fprintf(stderr, "Could not connect to buffer server at %s:%d.\n",hostname, port);
			return;
		}
	}
	
	ODM.enableStreaming();

	DataHandle hData = EE_DataCreate();
	EE_DataSetBufferSizeInSec(1.0);
	
	printf("Starting to transfer data - press [Escape] to quit\n");

	while (1) {
		if (conIn.checkKey() && conIn.getKey()==27) break;	
		ctrlServ.checkRequests(ODM);
	
		EE_DataUpdateHandle(0, hData);

		unsigned int nSamplesTaken=0;
		EE_DataGetNumberOfSample(hData,&nSamplesTaken);

		if (nSamplesTaken != 0) {
			double* data = new double[nSamplesTaken];
			double* dest = ODM.provideBlock(nSamplesTaken); 
			for (int i=0;i<TOTAL_CHANNELS;i++) {
				EE_DataGet(hData, targetChannelList[i], data, nSamplesTaken);
					
				for (unsigned int j=0;j<nSamplesTaken;j++) {
					dest[i + j*TOTAL_CHANNELS] = data[j];
				}
			}
			delete[] data;
			
			if (!ODM.handleBlock()) break;
			sampleCounter += nSamplesTaken;
			printf("Wrote %2i samples (%i total)\n", nSamplesTaken, sampleCounter);
		}
		dosleep(10);
	}
	EE_DataFree(hData);
}

					  
int main(int argc, char** argv) {
	//const unsigned short composerPort	= 1726;
	unsigned int samplingRate = 0;
	
	if (argc<2) {
		fprintf(stderr, "Usage:  emotiv2ft <config-file> [hostname=localhost [port=1972 [ctrlPort=8000]]]\n");
		fprintf(stderr, "Passing a minus (-) for the hostname tells this application to spawn its own buffer server\n");
		return 1;
	}
	
	if (argc>2) {
		strncpy(hostname, argv[2], sizeof(hostname));
	} else {
		strcpy(hostname, "localhost");
	}
	
	if (argc>3) {
		port = atoi(argv[3]);
	} else {
		port = 1972;
	}	

	if (argc>4) {
		ctrlPort = atoi(argv[4]);
	} else {
		ctrlPort = 8000;
	}	
	
	if (!ctrlServ.startListening(ctrlPort)) {
		fprintf(stderr, "Cannot listen on port %d for configuration commands\n", ctrlPort);
		return 1;
	}
	
	eEvent = EE_EmoEngineEventCreate();
	eState = EE_EmoStateCreate();
		
	if (EE_EngineConnect() != EDK_OK) {
		fprintf(stderr, "Emotiv Engine start up failed.");
		exit(1);
	}
	// TODO: provide alternative connection
	// ...	if (EE_EngineRemoteConnect(input.c_str(), composerPort) != EDK_OK) {
	
	// Loop here until acquisition can start and sampling rate is known
		
	printf("Waiting for the device to become ready - press [Escape] to quit\n");
	while (1) {
		int state = EE_EngineGetNextEvent(eEvent);

		if (state == EDK_OK) {
			EE_Event_t eventType = EE_EmoEngineEventGetType(eEvent);
			EE_EmoEngineEventGetUserId(eEvent, &userID);

			// Log the EmoState if it has been updated
			if (eventType == EE_UserAdded) {
				if (EE_DataGetSamplingRate(userID, &samplingRate) != EDK_OK) {
					fprintf(stderr, "Cannot retrieve sampling rate\n");
					break;
				} else {
					EE_DataAcquisitionEnable(userID,true);
					printf("EDK: User added, will now start acquisition.\n");
				}
				break;
			}
		}
		
		if (conIn.checkKey() && conIn.getKey()==27) break;		
	}
	
	if (samplingRate!=0) {
		acquisition(argv[1], samplingRate);
	}
	
	EE_EngineDisconnect();
	EE_EmoStateFree(eState);
	EE_EmoEngineEventFree(eEvent);

	return 0;
}


