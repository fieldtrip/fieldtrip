/*
 * Copyright (C) 2010, Stefan Klanke
 * Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
 * Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 */

#include <stdio.h>
#include <math.h>

#include <PixelDataGrabber.h>

int main(int argc, char *argv[]) {
	PixelDataGrabber pdg;
	char hostname[256];
	int port;
	char directory[256];
	
	timeBeginPeriod(1);
	
	if (argc>=2) {
		strncpy(hostname, argv[1], 256);
	} else {
		strncpy(hostname, "localhost", 256);
	}
	
	if (argc>=3) {
		port = atoi(argv[2]);
	} else {
		port = 1972;
	}
	
	if (argc>=4) {
		strncpy(directory, argv[3], 256);
	} else {
		strncpy(directory, "E:\\IMAGE", 256);
	}
	
	printf("Trying to connect to fieldtrip buffer on %s:%d\n", hostname, port);
	
	if (pdg.connectToFieldTrip(hostname, port)) {
		printf("OK\n");
	} else {
		printf("Failure!\n");
		exit(1);
	}
	
	printf("Trying to open directory %s for monitoring...\n", directory);
	if (pdg.monitorDirectory(directory)) {
		printf("OK\n");
	} else {
		printf("Failure!\n");
		exit(1);
	}
	
	printf("Starting to monitor... - press CTRL-C to cancel\n");
	
	while (pdg.run(200) >= 0) {}
	
	return 0;
}
