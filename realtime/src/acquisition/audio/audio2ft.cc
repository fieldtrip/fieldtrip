/*
 * Acquisition tool to stream audio signals to a FieldTrip buffer,
 * based on PortAudio for grabbing the data from the sound card.
 * 	
 * (C) 2010 S. Klanke
*/

#include <portaudio.h>
#include <stdio.h>
#include <OnlineDataManager.h>
#include <ConsoleInput.h>
#include <StringServer.h>

#define FSAMPLE 44100.0

PaStream *stream = NULL;
int numInputs = 0;
OnlineDataManager<float, float> *ODM = 0;
long numTotal = 0;


static int audioCallback(const void *input, void *output, unsigned long frameCount, 
                  const PaStreamCallbackTimeInfo *timeInfo, PaStreamCallbackFlags statusFlags,
                  void *userData) {
				  
	
	if (ODM==0) return paAbort;
	
	float *block = ODM->provideBlock(frameCount);
	if (block == 0) {
		fprintf(stderr, "Out of memory\n");
		return paAbort;
	}
	memcpy(block, input, frameCount * numInputs * sizeof(float));
		
	if (!ODM->handleBlock()) {
		fprintf(stderr, "Error in handling this data block - stopping\n");
		return paAbort;
	}
	numTotal += frameCount;
	printf("%4lu / %8lu\r", frameCount, numTotal);
	return paContinue;
}


int listDevices() {
	int numInputDevs = 0;
	
	int numDevices = Pa_GetDeviceCount();
	if (numDevices < 0) {
		fprintf(stderr, "Pa_GetDeviceCount() returned 0x%4x\n", numDevices );
		return numDevices;
	}
	
	if (numDevices == 0) {
		printf("No devices found by PortAudio.\n");
	}

	for (int i=0; i<numDevices; i++ ) {
		const PaDeviceInfo *deviceInfo;
		const PaHostApiInfo *apiInfo;
		
		deviceInfo = Pa_GetDeviceInfo(i);
		
		if (deviceInfo->maxInputChannels == 0) continue;
		
		apiInfo = Pa_GetHostApiInfo(deviceInfo->hostApi);
		
		printf("%2i: %s\n    %i inputs, %.0f Hz, %.3fs latency, driver = %s\n\n", 
			i, deviceInfo->name, 
			deviceInfo->maxInputChannels,
			deviceInfo->defaultSampleRate,
			deviceInfo->defaultLowInputLatency,
			apiInfo->name);
		numInputDevs++;
    }
    return numInputDevs;
}

int openDevice(int nr) {
	const PaDeviceInfo *deviceInfo = Pa_GetDeviceInfo(nr);
	if (deviceInfo == NULL) return 0;
	
	numInputs = deviceInfo->maxInputChannels;
	if (numInputs == 0) return 0;
	
	PaStreamParameters parIn;
	parIn.device = nr;
	parIn.channelCount = numInputs;
	parIn.sampleFormat = paFloat32;
	parIn.suggestedLatency = 0.0;
	parIn.hostApiSpecificStreamInfo = NULL;

	PaError err = Pa_OpenStream(&stream, &parIn, NULL, FSAMPLE, 0, paNoFlag, audioCallback, NULL);
	if (err != paNoError) return 0;
	
	const PaStreamInfo *info = Pa_GetStreamInfo(stream);
	
	printf("Opened autio stream with %f Hz, %fs IO-latency\n", info->sampleRate, info->inputLatency);
	return numInputs;
}


int main(int argc, char *argv[]) {
	PaError err;
	ConsoleInput conIn;
	char hostname[256] = "localhost";
	int port = 1972;
   int device;
	
	err = Pa_Initialize();
	if( err != paNoError ) {
		fprintf(stderr, "Pa_Initialize() returned error %i\n", err);
		return 0;
	}
	
	if (argc<2) {
		listDevices();
		printf("\nUsage: audio2ft deviceNr [hostname [port]]\n");
		printf("The 2nd parameter (device number) must be one of the numbers listed above.\n");
		goto cleanup;
	}
	
	device = atoi(argv[1]);
	if (!openDevice(device)) {
		goto cleanup;
	}
	
	if (argc>2) {
		strncpy(hostname, argv[2], sizeof(hostname));
	}
	if (argc>3) {
		port = atoi(argv[3]);
	}
	
	ODM = new OnlineDataManager<float, float>(0, numInputs, FSAMPLE, FSAMPLE);
	
	if (0==strcmp(hostname, "-")) {
		if (!ODM->useOwnServer(port)) {
			fprintf(stderr, "Could not spawn buffer server on port %i.\n", port);
			goto cleanup;
		}
	} else {
		if (!ODM->connectToServer(hostname, port)) {
			fprintf(stderr, "Could not connect to buffer server on %s:%i.\n", hostname,port);
			goto cleanup;
		}
	}
	
	if (ODM->configureFromFile("config.txt") != 0) {
		fprintf(stderr, "Configuration file config.txt not found or invalid\n");
		goto cleanup;
	}	

	ODM->enableStreaming();
	
	err = Pa_StartStream(stream);
	if (err != paNoError) {
		fprintf(stderr, "Could not start stream\n");
		goto cleanup;
	}
	
	printf("Starting - press ESC to quit\n");
	while (1) {
		if (conIn.checkKey()) {
			int c = conIn.getKey();
			if (c==27) break; // quit
		}
		conIn.milliSleep(100);
	}
	
cleanup:
	if (stream!=NULL) {
		err = Pa_StopStream(stream);   
		err = Pa_CloseStream(stream);
	}
	err = Pa_Terminate();
	delete ODM;
	return 0;
}

