/** NeuroSky ThinkGear acqusition tool to stream data to a FieldTrip buffer,
	and write data to one or multiple GDF files (if you ever reach the size limit...).
	(C) 2010 S. Klanke
*/
#include <serial.h>

#include <OnlineDataManager.h>
#include <ConsoleInput.h>
#include <StringServer.h>

/**
ThinkCap is a 7 channel dry-sensor EEG device. It is basically using
similar communication protocol with the rest of the other NeuroSky
devices.  They differ only in the number of channels available and
sampling rate, i.e. Thinkcap: 7ch/250Hz; Mindset: 1ch/512Hz.
*/

#define NUMCHANS    7
#define FSAMPLE     250.0
#define MAXBLOCK    20

unsigned char packet[170];


int tryReadPacket(SerialPort *SP) {
	unsigned char byte;
	unsigned char checksum;
	int len;
	int nr, syncBytes = 0;
   
sync:
	while (syncBytes < 2) {
		nr = serialRead(SP, 1, &byte);
		if (nr<0) return nr;
		if (nr==0) continue;
		if (byte == 0xAA) {
			syncBytes++;
		} else {
			syncBytes = 0;
		}
	}
   
	do {
		nr = serialRead(SP, 1, &byte);
		if (nr < 0) return -1;
		if (nr==0) continue;
		if (byte > 0xAA) {
			syncBytes = 0;
			goto sync;
		}
	} while (!nr || byte==0xAA);

	len = byte;
	checksum = 0;

	for (int i=0;i<len;i++) {
		do {
			nr = serialRead(SP, 1, packet + i);
			if (nr==0) continue;
			if (nr < 0) return -1;
		} while (!nr);
		checksum += packet[i];
	}
	checksum = ~checksum;
	do {
		nr = serialRead(SP, 1, &byte);
		if (nr < 0) return -1;
	} while (!nr);
	if (byte != checksum) {
		syncBytes = 0;
		goto sync;
	}
	return len;
}

int main(int argc, char *argv[]) {
	ConsoleInput conIn;
	StringServer ctrlServ;
	SerialPort SP;
	char hostname[256];
	int port, ctrlPort;
	int counter = 0;
	OnlineDataManager<short, float> ODM(0, NUMCHANS, FSAMPLE);

	if (argc<3) {
		printf("Usage: thinkgear2ft <device> <config-file> [hostname=localhost [port=1972 [[ctrlPort=8000]]]\n");
		return 0;
	}

	if (argc>3) {
		strncpy(hostname, argv[3], sizeof(hostname));
	} else {
		strcpy(hostname, "localhost");
	}

	if (argc>4) {
		port = atoi(argv[4]);
	} else {
		port = 1972;
	}

	if (argc>5) {
		ctrlPort = atoi(argv[5]);
	} else {
		ctrlPort = 8000;
	}

	if (!serialOpenByName(&SP, argv[1])) {
		fprintf(stderr, "Could not open serial port %s\n", argv[1]);
		return 1;
	}

	// last parameter is timeout in 1/10 of a second
	if (!serialSetParameters(&SP, 57600, 8, 0, 1, 1)) {
		fprintf(stderr, "Could not modify serial port parameters\n");
		return 1;
	}

	if (!strcmp(hostname, "-")) {
		if (!ODM.useOwnServer(port)) {
			fprintf(stderr, "Could not spawn buffer server on port %d.\n",port);
			return 0;
		}
	} else {
		if (!ODM.connectToServer(hostname, port)) {
			fprintf(stderr, "Could not connect to buffer server at %s:%d.\n",hostname, port);
			return 0;
		}
	}
	if (ODM.configureFromFile(argv[2]) != 0) {
		fprintf(stderr, "Configuration file is invalid\n");
		return 0;
	}

	ctrlServ.startListening(ctrlPort);
	ODM.enableStreaming();

	printf("Starting - press ESC to quit\n");
	while (1) {
		short samples[MAXBLOCK*NUMCHANS];
		int nSamples = 0;

		if (conIn.checkKey()) {
			int c = conIn.getKey();
			if (c==27) break; // quit
		}

		ctrlServ.checkRequests(ODM);

		while (nSamples < MAXBLOCK && serialInputPending(&SP)>=6+2*NUMCHANS) {
			int len = tryReadPacket(&SP);

			if (len != 2*NUMCHANS+2 || packet[0] != 0xB0 || packet[1] != 2*NUMCHANS) {
				printf("Unrecognized packet: %2d %02X %02X\n", len, packet[0], packet[1]);
				continue;
			}
			short *dest = samples + nSamples * NUMCHANS;

			for (int i=0;i<NUMCHANS;i++) {
				short high = packet[2+i*2];
				short low  = packet[3+i*2];
				if (high & 0x10 && low==3) {
					low=2;
				}
				high &= 0x0F;
				short val = (high << 8) | low;
				dest[i] = val;
			}
			nSamples++;
		} 
		
		if (nSamples == 0) {
			conIn.milliSleep(10);
		} else {
			short *block = ODM.provideBlock(nSamples);
			memcpy(block, samples, nSamples*NUMCHANS*sizeof(short));
			if (!ODM.handleBlock()) {
				fprintf(stderr, "Error in handling this data block - stopping\n");
				break;
			}
			counter+=nSamples;
			printf("Samples: %i\r", counter);
		}
	}

	ODM.disableStreaming();

	serialClose(&SP);
	return 0;
}
