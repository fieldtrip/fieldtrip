/*
   ModularEEG acqusition tool to stream data to a FieldTrip buffer,
   and write data to one or multiple GDF files (if you ever reach the size limit...).

   (C) 2010 S. Klanke
 */

#include <stdio.h>
#include <pthread.h>

#include <socketserver.h>
#include <serial.h>

#include <FtBuffer.h>
#include <LocalPipe.h>
#include <GdfWriter.h>
#include <ConsoleInput.h>

#define NUM_HW_CHAN 6
#define FSAMPLE     256
#define PACKET_LEN  17

int64_t maxFileSize = 1024*1024*1024; // 1GB

// lock-free FIFO for saving thread
int16_t *rbInt;
int rbIntSize, rbIntChans;
int rbIntWritePos;
int rbIntReadPos;

// pipe or socketpair for inter-thread communication
LocalPipe threadpipe;
// GDF writing object

GDF_Writer *gdfWriter = NULL;

// Name of file to write to
char baseFilename[1024];
char curFilename[1024];
int fileCounter = 0;

unsigned char serialBuffer[1024];
int leftOverBytes = 0;

/** Read data from serial port until 0xA5 0x5A shows up.
  Those two bytes will be written to the beginning of
  the global array "serialBuffer". Also, leftOverBytes
  will be set to 2 if this functions returns successfully.
  Returns 2 on success, 0 or 1 on error.
 */
int readSyncBytes(SerialPort *SP) {
		// printf("Looking for sync bytes 0xA5 0x5A (%c%c)\n", 0xA5, 0x5A);

		// read bytes until we get 0xA5,0x5A,...
		for (int iter = 0;iter < 200; iter++) {
				unsigned char byte;
				int nr;

				nr = serialRead(SP, 1, &byte);
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
		return leftOverBytes;
}

void *savingThreadFunction(void *arg) {
		int writePtr,n;
		int64_t fileSize = 256*(1+rbIntChans);

		while (1) {
				int newSamplesA, newSamplesB;
				int16_t *rbIntPtr;
				int64_t newSize;

				n = threadpipe.read(sizeof(int), &writePtr);
				if (n!=sizeof(int)) {
						// this should never happen for blocking sockets/pipes
						fprintf(stderr, "Unexpected error in pipe communication\n");
						break;
				}
				if (writePtr < 0) {
						printf("\nSaving thread received %i - exiting...\n", writePtr);
						break;
				}

				rbIntPtr = rbInt + rbIntReadPos*rbIntChans;
				if (writePtr > rbIntReadPos) {
						newSamplesA = writePtr - rbIntReadPos;
						newSamplesB = 0;
				} else {
						newSamplesA = rbIntSize - rbIntReadPos;
						newSamplesB = writePtr;
				}

				newSize = fileSize + (newSamplesA+newSamplesB) * rbIntChans * sizeof(int32_t);
				if (newSize > maxFileSize) {
						gdfWriter->close();

						fileCounter++;

						snprintf(curFilename, sizeof(curFilename), "%s_%i.gdf", baseFilename, fileCounter);
						if (!gdfWriter->createAndWriteHeader(curFilename)) {
								fprintf(stderr, "Error: could not create GDF file %s\n", curFilename);
								break;
						}
						newSize = 256*(1+rbIntChans) + (newSamplesA+newSamplesB) * rbIntChans * sizeof(int16_t);
				}

				gdfWriter->addSamples(newSamplesA, rbIntPtr);
				// mark as empty
				for (int i=0;i<newSamplesA;i++) {
						rbIntPtr[i*rbIntChans] = -1;
				}
				if (newSamplesB > 0) {
						gdfWriter->addSamples(newSamplesB, rbInt);
						for (int i=0;i<newSamplesB;i++) {
								rbInt[i*rbIntChans] = -1;
						}
				}
				rbIntReadPos = writePtr;
				fileSize = newSize;
		}
		return NULL;
}


int main(int argc, char *argv[]) {
		int port, ftSocket, err;
		int sampleCounter = 0;
		char hostname[256];
		FtBufferRequest req;
		FtBufferResponse resp;
		FtEventList   eventChain;
		FtSampleBlock sampleBlock(DATATYPE_FLOAT32);
		int16_t switchState = 0;
		ft_buffer_server_t *ftServer = NULL;
		pthread_t savingThread;
		ConsoleInput ConIn;
		SerialPort SP;
		int keepRunning = 1;
		int numTimeouts = 0;

		int16_t sampleData[NUM_HW_CHAN * FSAMPLE]; // holds up to 1 seconds of data
		int16_t switchData[FSAMPLE];       // again, 1 second of data (switches)

		if (argc<3) {
				printf("Usage: modeeg2ft_2chn <device> <gdf-file> [hostname=localhost [port=1972]]\n\n");
				printf("Remarks\n");
				printf(" 1) <device> must be your serial port, e.g, COM3: on Windows, or /dev/ttyS0 on Linux\n");
				printf(" 2) putting a minus (-) for the <gdf-file> parameter disables saving to disk.\n");
				printf(" 3) putting a minus (-) for the hostname parameter creates a new buffer server locally.\n");
				return 0;
		}

		strcpy(baseFilename, argv[2]);
		if (strcmp(baseFilename, "-")) {
				gdfWriter = new GDF_Writer(3, FSAMPLE, GDF_INT16);
				gdfWriter->setLabel(0, "Switches");
				gdfWriter->setLabel(1, "ModEeg1");
				gdfWriter->setLabel(2, "ModEeg2");
				gdfWriter->setPhysicalLimits(1, -256.0, 255.5);
				gdfWriter->setPhysicalLimits(2, -256.0, 255.5);
				gdfWriter->setPhysDimCode(1, GDF_MICRO + GDF_VOLT);
				gdfWriter->setPhysDimCode(2, GDF_MICRO + GDF_VOLT);


				char *lastDotPos = strrchr(baseFilename, '.');
				// cut off .gdf suffix, if given
				if (lastDotPos != NULL) {
						if (!strcasecmp(lastDotPos, ".gdf")) {
								*lastDotPos = 0;
						}
				}
				strcpy(curFilename, baseFilename);
				strcat(curFilename, ".gdf");
				if (!gdfWriter->createAndWriteHeader(curFilename)) {
						fprintf(stderr, "Could not open GDF file %s for writing\n", curFilename);
						return 1;
				}

				// allocate and initialise internal ring buffer for saving to GDF
				rbIntSize = 2*FSAMPLE;	// enough space for 2 seconds of data
				rbIntChans = 3;
				rbInt = new int16_t[rbIntChans*rbIntSize];
				rbIntWritePos = 0;
				rbIntReadPos = 0;
				for (int i=0;i<rbIntSize;i++) rbInt[i*rbIntChans]=-1;
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

		if (strcmp(hostname, "-")) {
				ftSocket = open_connection(hostname, port);
				if (ftSocket < 0) {
						fprintf(stderr, "Could not connect to FieldTrip buffer on %s:%i\n", hostname, port);
						return 1;
				}
		} else {
				ftServer = ft_start_buffer_server(port, NULL, NULL, NULL);
				if (ftServer == NULL) {
						fprintf(stderr, "Could not spawn TCP server on port %i\n", port);
						return 1;
				}
				else {
						fprintf(stderr, "Spawned TCP server on port %i\n", port);
				}
				ftSocket = 0;
		}

		// write header to FieldTrip buffer
		req.prepPutHeader(2, DATATYPE_FLOAT32, (float) FSAMPLE);
		err = clientrequest(ftSocket, req.out(), resp.in());
		if (err || !resp.checkPut()) {
				fprintf(stderr, "Could not write header to FieldTrip buffer\n.");
				return 1;
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

		if (gdfWriter != NULL) {
				if (pthread_create(&savingThread, NULL, savingThreadFunction, NULL)) {
						fprintf(stderr, "Could not spawn GDF saving thread.\n");
						return 1;
				}
		}

		if (readSyncBytes(&SP) != 2) {
				fprintf(stderr, "Could not read synchronisation bytes from ModularEEG\n");
				goto cleanup;
		}

		printf("Got synchronization bytes - starting acquisition\n");
		printf("\nPress <Esc> to quit\n\n");

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
						if (++numTimeouts > 250) {
								// write a timeout event

								eventChain.clear();
								eventChain.add(sampleCounter, "TIMEOUT", 0);

								clientrequest(ftSocket, eventChain.asRequest(), resp.in());
								// silently ignore errors here

								fprintf(stderr, "Timeout -- re-opening serial port\n");
								serialClose(&SP);

								if (!serialOpenByName(&SP, argv[1])) {
										fprintf(stderr, "Could not open serial port %s\n", argv[1]);
										break;
								}
								if (!serialSetParameters(&SP, 57600, 8, 0, 1, 0)) {
										fprintf(stderr, "Could not modify serial port parameters\n");
										break;
								}
								if (readSyncBytes(&SP) != 2) {
										fprintf(stderr, "Got synchronization bytes - re-starting acquisition\n");
										continue;
								} else {
										fprintf(stderr, "Could not read synchronization bytes - exiting.\n");
										break;
								}
								numTimeouts = 0;
						} else {
								usleep(20000);
						}
						continue;
				}
				numTimeouts = 0; // we read something, so reset timeout counter

				maxReadNow = sizeof(serialBuffer) - leftOverBytes;
				if (numPending > maxReadNow) {
						numPending = maxReadNow;
				}

				numRead = serialRead(&SP, numPending, serialBuffer + leftOverBytes);
				if (numRead != numPending) {
						fprintf(stderr, "Error when reading from serial port - exiting\n");
						break;
				}

				// remove any events in our list
				eventChain.clear();

				numTotal = leftOverBytes + numRead;
				numSamples = numTotal / PACKET_LEN;

				if (numSamples > FSAMPLE) {
						fprintf(stderr, "Received too much data from serial port - exiting.\n");
						break;
				}

				if (numSamples == 0) {
						leftOverBytes += numRead;
						continue;
				}

				printf("%i sample(s) , %i bytes\n", numSamples, numTotal);

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
						// if the switch byte changes its state, add this as an event
						if (switchData[j] != switchState) {
								switchState = switchData[j];
								if (switchState!=0) {
										eventChain.add(sampleCounter + j, "Switch", switchState);
								}
						}
				}

				// if saving is enabled, append to GDF file
				if (gdfWriter != NULL) {
						for (int j=0;j<numSamples;j++) {
								int doff = rbIntWritePos * rbIntChans;

								if (rbInt[doff] != -1) {
										fprintf(stderr, "Error: saving thread does not keep up with load\n");
										break;
								}

								rbInt[doff+0] = switchData[j];
								rbInt[doff+1] = sampleData[0 + j*NUM_HW_CHAN];
								rbInt[doff+2] = sampleData[1 + j*NUM_HW_CHAN];

								if (++rbIntWritePos == rbIntSize) rbIntWritePos = 0;
						}
						threadpipe.write(sizeof(int), static_cast<void *>(&rbIntWritePos));
				}

				// streaming stuff
				if (ftSocket != -1) {
						float *dest = (float *) sampleBlock.getMatrix(2, numSamples);
						if (dest==NULL) {
								fprintf(stderr, "Out of memory\n");
								break;
						}
						for (int j=0;j<numSamples;j++) {
								dest[0 + 2*j] = (float) (sampleData[0 + j*NUM_HW_CHAN] - 512) * 0.5;
								dest[1 + 2*j] = (float) (sampleData[1 + j*NUM_HW_CHAN] - 512) * 0.5;
						}
						err = clientrequest(ftSocket, sampleBlock.asRequest(), resp.in());
						if (err || !resp.checkPut()) {
								fprintf(stderr, "Could not write samples to FieldTrip buffer\n.");
						}

						if (eventChain.count() > 0) {
								err = clientrequest(ftSocket, eventChain.asRequest(), resp.in());
								if (err || !resp.checkPut()) {
										fprintf(stderr, "Could not write events to FieldTrip buffer\n.");
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

cleanup:
		if (gdfWriter != NULL) {
				int quitValue = -1;
				threadpipe.write(sizeof(int), &quitValue);
		}

		serialClose(&SP);
		if (ftSocket > 0) close_connection(ftSocket);
		if (ftServer != NULL) ft_stop_buffer_server(ftServer);
		if (gdfWriter != NULL) {
				gdfWriter->close();
				delete gdfWriter;
		}

		return 0;
}
