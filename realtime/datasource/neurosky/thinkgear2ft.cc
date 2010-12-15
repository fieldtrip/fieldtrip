/** NeuroSky ThinkGear acqusition tool to stream data to a FieldTrip buffer,
	and write data to one or multiple GDF files (if you ever reach the size limit...).
	(C) 2010 S. Klanke
*/
#include <serial.h>

#include <OnlineDataManager.h>
#include <ConsoleInput.h>
#include <StringServer.h>

#define NUMCHANS    7
#define FSAMPLE     128.0
#define MAXBLOCK    20

/*
	// printf("Looking for sync bytes 0xA5 0x5A (%c%c)\n", 0xA5, 0x5A);
	
	// read bytes until we get 0xA5,0x5A,...
	for (int iter = 0;iter < 200; iter++) {
		unsigned char byte;
		int nr;
		
		nr = serialRead(&SP, 1, &byte);
		if (nr<0) {
			fprintf(stderr, "Error when reading from serial port - exiting\n");
			return 1;
		} 
		if (nr==0) {
			printf(".");
			usleep(10000);	// sleep for 10 ms if no byte received
			continue;
		}
		
		//printf("%02X %c\n", byte, byte);
		
		if (byte == 0xA5) {
			serialBuffer[0] = byte;
			leftOverBytes = 1;
		} else if (leftOverBytes == 1 && byte == 0x5A) {
			serialBuffer[1] = byte;
			leftOverBytes = 2;
			break; // success !
		} else {
			leftOverBytes = 0;
		}
	}
	if (leftOverBytes != 2) {
		fprintf(stderr, "Could not read synchronisation bytes from ModularEEG\n");
		goto cleanup;
	}
	
	printf("Got synchronization bytes - starting acquisition\n");
	
	while (keepRunning) {
		int numRead, numTotal, numSamples, numPending, maxReadNow;
		
		if (ConIn.checkKey()) {
			int c = ConIn.getKey();
			if (c==27) break; // quit
		}
		
		numPending = serialInputPending(&SP);
		if (numPending < 0) {
			fprintf(stderr, "Error when reading from serial port - exiting\n");
			break;
		}
		if (numPending == 0) {
			usleep(10000);
			continue;
		}
		
		maxReadNow = sizeof(serialBuffer) - leftOverBytes;
		if (numPending > maxReadNow) {
			numPending = maxReadNow;
		}
		
		numRead = serialRead(&SP, numPending, serialBuffer + leftOverBytes);
		if (numRead != numPending) {
		    fprintf(stderr, "Error when reading from serial port - exiting\n");
			break;
		}	
		eventChain.clear();
		
		numTotal   = leftOverBytes + numRead;
		numSamples = numTotal / PACKET_LEN;
		
		if (numSamples > FSAMPLE) {
			fprintf(stderr, "Received too much data from serial port - exiting.\n");
			break;
		}
		
		if (numSamples == 0) {
			leftOverBytes += numRead;
			continue;
		}
		
		// first decode into switchData + sampleData arrays
		for (int j=0;j<numSamples;j++) {
			int soff = j*PACKET_LEN;
			int doff = j*NUM_HW_CHAN;
			
			if (serialBuffer[soff] != 0xA5 || serialBuffer[soff+1] != 0x5A) {
				fprintf(stderr, "ModularEEG out of sync in sample %i - exiting.\n", sampleCounter + j);
				keepRunning = 0;
				break;
			}
				
			for (int i=0;i<NUM_HW_CHAN;i++) {
				sampleData[i+doff] = serialBuffer[soff+4+2*i]*256 + serialBuffer[soff+5+2*i];
			}
		
			switchData[j] = serialBuffer[soff+16];
			if (switchData[j] != switchState) {
				switchState = switchData[j];
				if (switchState!=0) {
					eventChain.add(sampleCounter + j, "Switch", switchState);
				}
			}
		}
	
		sampleCounter += numSamples;
		leftOverBytes = numTotal - numSamples * PACKET_LEN;
		
		// copy left-over bytes to the beginning of the serialBuffer
		if (leftOverBytes > 0) {
			memcpy(serialBuffer, serialBuffer + numSamples*PACKET_LEN, leftOverBytes);
		}
	}


	
}
*/



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
      if (byte == 0xAA) {
         syncBytes++;
      } else {
         syncBytes = 0;
      }
   }
   
   do {
      nr = serialRead(SP, 1, &byte);
      if (nr < 0) return -1;
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
   int port;
   int counter = 0;
   OnlineDataManager<short, float> ODM(0, NUMCHANS, FSAMPLE, GDF_INT16, DATATYPE_FLOAT32);
	   
   
	if (argc<3) {
		printf("Usage: thinkgear2ft <device> <config-file> [hostname=localhost [port=1972]]\n");
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
	
	if (!serialOpenByName(&SP, argv[1])) {
		fprintf(stderr, "Could not open serial port %s\n", argv[1]);
		return 1;
	}
	
	// last parameter is timeout in 1/10 of a second
	if (!serialSetParameters(&SP, 57600, 8, 0, 1, 0)) {
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
	
	ctrlServ.startListening(8000);
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
      
      do {
         int len = tryReadPacket(&SP);
      
         if (len != 2*NUMCHANS+2 || packet[0] != 0xB0 || packet[1] != 2*NUMCHANS) {
            printf("Unrecognized packet: %2d %02X %02X\n", len, packet[0], packet[1]);
            continue;
         }
         short *dest = samples + nSamples * NUMCHANS;
         
         for (int i=0;i<NUMCHANS;i++) {
            short high = packet[2+i*2];
            short low  = packet[3+i*2];
            if (high & 0x10) {
               low--;
            }
            high &= 0x8F;
            short val = (high << 8) | low;
            dest[i] = val;
         }
         nSamples++;
      } while (nSamples < MAXBLOCK && serialInputPending(&SP)>=6+2*NUMCHANS);
      
      short *block = ODM.provideBlock(nSamples);
      memcpy(block, samples, nSamples*NUMCHANS*sizeof(short));
      if (!ODM.handleBlock()) {
         fprintf(stderr, "Error in handling this data block - stopping\n");
         break;
      }
      counter+=nSamples;
      printf("Samples: %i\r", counter);
	}
	
	ODM.disableStreaming();

  	serialClose(&SP);
	return 0;
}