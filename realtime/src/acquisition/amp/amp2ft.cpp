/* Copyright (C) 2013 Federico Raimondo
 * Applied Artificial Intelligence Lab
 * Computer Sciences Department
 * University of Buenos Aires, Argentina
 *
 * This file is part of Amp2ft
 *
 * Amp2ft is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Amp2ft is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Amp2ft.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <AmpServerClient.h>
#include <stdio.h>
#include <unistd.h>
#include <ConsoleInput.h>
#include <StringServer.h>
#include <errno.h>



struct AmpServerClientConfig config;
StringServer ctrlServ;
char * configFile;
char * gdfFile;
ConsoleInput conIn;
bool haveGDF;

char* getParam(const char * needle, char* haystack[], int count) {
	int i = 0;
	for (i = 0; i < count; i ++) {
		if (strcmp(needle, haystack[i]) == 0) {
			if (i < count -1) {
				return haystack[i+1];
			}
		}
	}
	return 0;
}


int isParam(const char * needle, char* haystack[], int count) {
	int i = 0;
	for (i = 0; i < count; i ++) {
		if (strcmp(needle, haystack[i]) == 0) {
			return 1;
		}
	}
	return 0;
}


void acquisition(AmpServerClient& client, int numHwChan, int fSample) {
	int sampleCounter = 0;
	DPRINTF("Starting ODM with sample frequency %d, nchannels %d\n", fSample, numHwChan);
	OnlineDataManager<int, float> ODM(0, numHwChan, (float) fSample);

	if (ODM.configureFromFile(configFile) != 0) {
		fprintf(stderr, "Configuration %s file is invalid\n", configFile);
		return;
	} else {
		printf("Streaming %i out of %i channels\n", ODM.getSignalConfiguration().getStreamingSelection().getSize(), numHwChan);
	}
	DPRINTF("Connecting to buffer server at %s:%d.\n",config.hostname, config.port);

	if (!strcmp(config.hostname, "-")) {
		if (!ODM.useOwnServer(config.port)) {
			fprintf(stderr, "Could not spawn buffer server on port %d.\n",config.port);
			return;
		}
	} else {
		if (!ODM.connectToServer(config.hostname, config.port)) {
			fprintf(stderr, "Could not connect to buffer server at %s:%d.\n",config.hostname, config.port);
			return;
		}
	}
	DPRINTF("Connected!\n");

	/*
	 * As I cannot use 24 bits integers, the maximum is established
	 * as the maximum voltage for 32 bits numbers, keeping the relation
	 * 1 unidad = 0.488281308 micro volts
	 * For the minimum, the relation is 0.488 micro volts.
	 */

	ODM.setPhysicalLimits(-4.096, 4.096); //Page 87 of the manual
	ODM.setPhysicalDimCode(GDF_MICRO + GDF_VOLT);
	ODM.setDigitalLimits(-8388608, 8388607);
	ODM.setSlopeAndOffset(0.488281308, 0);
	if (haveGDF) ODM.setFilename(gdfFile);


	DPRINTF("Enabling streaming\n");
	if (!ODM.enableStreaming()) return;
	DPRINTF("Enabled\n");

	if (haveGDF) {
		printf("\nPress <Esc> to quit, <S> to enable saving, <D> to disable saving\n\n");
	} else {
		printf("\nPress <Esc> to quit\n\n");
	}

	client.start();
	while (1) {
		if (conIn.checkKey()) {
			int c = conIn.getKey();
			if (c==27) break; // quit

			if (haveGDF) {
				if (c=='s' || c=='S') {
					if (!ODM.enableSaving()) {
						fprintf(stderr, "Cannot enable saving to GDF\n");
						break;
					}
				}

				if (c=='d' || c=='D') {
					ODM.disableSaving();
				}
			}
		}
		ctrlServ.checkRequests(ODM);
		int readed;
		char * ptr;
		int subsample = client.getSubsample();
		if ((readed = client.checkNewData()) > 0) {
			int topass = readed/subsample;
			int passed = 0;
			if ((readed+sampleCounter) % subsample < sampleCounter % subsample) {
				topass++;
			}
			if (topass > 0) {
				ptr = (char*)ODM.provideBlock(topass);
				if (ptr == NULL) {
					fprintf(stderr, "Out of memory\n");
					break;
				}
			}
			passed = client.readNewData((int32_t*)ptr, topass, ODM.getEventList());
			if (passed > 0) {
				if (!ODM.handleBlock()) break;
			}

			if (passed != topass) {DPRINTF("WTF\n");}

			sampleCounter += passed;
		}
		usleep(10);
	}

	 client.stop();

}

char programname[2048];

void help(void) {
	printf("This is %s\n\n", programname);
	printf("Usage: %s <config-file> [OPTIONS]\n", programname);
	printf("\nOptional parameters:\n");
	printf("\t-gdf <gdf-file>\tEnable saving to <gdf-file> (default: disabled)\n");
	printf("\t-host <ip>\tUse fieldtrip buffer in <ip> (default: create a new buffer in localhost)\n");
	printf("\t-port <port>\tUse <port> to connect to fieldtrip buffer (default: 1972)\n");
	printf("\t-amph <ip>\tSet AmpServer ip address (default: localhost)\n");
	printf("\t-ampcp <port>\tSet AmpServer command port (default: 9877)\n");
	printf("\t-ampsp <port>\tSet AmpServer stream port (default: 9877)\n");
	printf("\t-sf <x>\tSet sampling frequency to x (default: 1000)\n");
	printf("\t-h\t Print this help\n");
}


int main(int argc, char *argv[]) {
		check_datatypes();

		strncpy((char*)&programname, argv[0], 2048);

		if (isParam("-h", argv, argc)) {
			help();
			return 0;
		}

		if (argc<2) {
			help();
			return 0;
		}

		haveGDF = false;

		configFile =strdup(argv[1]);
		if (access(configFile, R_OK)) {
			fprintf(stderr, "Cannot open configuration file %s\n", configFile);
			return 1;
		}

		if (isParam("-gdf", argv, argc)) {
			gdfFile = getParam("-gdf", argv, argc);
			haveGDF = true;
			if (access(gdfFile, W_OK)) {
				int errsv = errno;
				if (errsv != ENOENT) {
					fprintf(stderr, "Cannot open GDF file %s (%s)\n", gdfFile, strerror(errno));
					return 1;
				}
			}
		}

		if (isParam("-host", argv, argc)) {
			config.hostname = strdup(getParam("-host", argv, argc));
		} else {
			config.hostname = strdup("-");
		}

		if (isParam("-port", argv, argc)) {
			config.port = atoi(getParam("-port", argv, argc));
		} else {
			config.port = 1972;
		}

		if (isParam("-amph", argv, argc)) {
			config.amphostname = strdup(getParam("-amph", argv, argc));
		} else {
			config.amphostname = strdup("localhost");
		}

		if (isParam("-ampcp", argv, argc)) {
			config.ampcommandport = atoi(getParam("-ampcp", argv, argc));
		} else {
			config.ampcommandport = 9877;
		}

		if (isParam("-ampsp", argv, argc)) {
			config.ampstreamport = atoi(getParam("-ampsp", argv, argc));
		} else {
			config.ampstreamport = 9879;
		}

		if (isParam("-sf", argv, argc)) {
			config.sfreq = atoi(getParam("-sf", argv, argc));
		} else {
			config.sfreq = 1000;
		}

		printf("Starting AmpServerClient with connection to %s on ports %d (command) and %d (stream)\n", \
			config.amphostname, config.ampcommandport, config.ampstreamport);

		AmpServerClient client(&config);
		client.connectClient();

		if (!ctrlServ.startListening(8000)) {
			fprintf(stderr, "Cannot listen on port %d for configuration commands\n", 8000);
			return 1;
		}

		acquisition(client, client.getNumChannels(), client.getSamplingFreq());


		if (argc > 2) free(gdfFile);
		free(config.hostname);
		free(config.amphostname);
		free(configFile);

        return 0;
}
