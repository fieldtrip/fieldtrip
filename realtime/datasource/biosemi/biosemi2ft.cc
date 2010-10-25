/** Biosemi acqusition tool to stream (downsampled) data to a FieldTrip buffer,
	and write full-rate data to one or multiple GDF files.
	(C) 2010 S. Klanke
*/
#include <BioSemiClient.h>
#include <SignalConfiguration.h>
#include <stdio.h>
//#include <signal.h>
#include <FtBuffer.h>
#include <socketserver.h>
#include <pthread.h>
#include <LocalPipe.h>
#include <GdfWriter.h>
#include <MultiChannelFilter.h>
#include <ConsoleInput.h>

/* TODOs: 
	- add proper clean-up (instead of return 1; in main routine)
	- add TCP command server mechanism to START and STOP acquisition etc.
	- add documentation
*/

double chebyCoefsB[32][5] = {
  {  1.00000000,  0.00000000,  0.00000000,  0.00000000,  0.00000000},
  {  0.10289966, -0.12065116,  0.18416503, -0.12065116,  0.10289966},
  {  0.08816900, -0.20929048,  0.27825091, -0.20929048,  0.08816900},
  {  0.08471123, -0.25175920,  0.34686785, -0.25175920,  0.08471123},
  {  0.08421877, -0.27799101,  0.39317078, -0.27799101,  0.08421877},
  {  0.08467437, -0.29611308,  0.42573086, -0.29611308,  0.08467437},
  {  0.08543489, -0.30946972,  0.44966776, -0.30946972,  0.08543489},
  {  0.08626647, -0.31975126,  0.46793316, -0.31975126,  0.08626647},
  {  0.08707953, -0.32792119,  0.48229840, -0.32792119,  0.08707953},
  {  0.08784086, -0.33457419,  0.49387754, -0.33457419,  0.08784086},
  {  0.08854053, -0.34009907,  0.50340192, -0.34009907,  0.08854053},
  {  0.08917840, -0.34476144,  0.51136970, -0.34476144,  0.08917840},
  {  0.08975834, -0.34874923,  0.51813121, -0.34874923,  0.08975834},
  {  0.09028555, -0.35219930,  0.52393959, -0.35219930,  0.09028555},
  {  0.09076551, -0.35521371,  0.52898216, -0.35521371,  0.09076551},
  {  0.09120337, -0.35787021,  0.53340038, -0.35787021,  0.09120337},
  {  0.09160386, -0.36022903,  0.53730302, -0.36022903,  0.09160386},
  {  0.09197114, -0.36233760,  0.54077505, -0.36233760,  0.09197114},
  {  0.09230890, -0.36423375,  0.54388381, -0.36423375,  0.09230890},
  {  0.09262037, -0.36594808,  0.54668333, -0.36594808,  0.09262037},
  {  0.09290834, -0.36750555,  0.54921746, -0.36750555,  0.09290834},
  {  0.09317528, -0.36892674,  0.55152214, -0.36892674,  0.09317528},
  {  0.09342333, -0.37022882,  0.55362712, -0.37022882,  0.09342333},
  {  0.09365435, -0.37142615,  0.55555723, -0.37142615,  0.09365435},
  {  0.09387001, -0.37253088,  0.55733337, -0.37253088,  0.09387001},
  {  0.09407174, -0.37355337,  0.55897321, -0.37355337,  0.09407174},
  {  0.09426083, -0.37450246,  0.56049185, -0.37450246,  0.09426083},
  {  0.09443840, -0.37538579,  0.56190222, -0.37538579,  0.09443840},
  {  0.09460545, -0.37620996,  0.56321550, -0.37620996,  0.09460545},
  {  0.09476289, -0.37698073,  0.56444135, -0.37698073,  0.09476289},
  {  0.09491149, -0.37770312,  0.56558823, -0.37770312,  0.09491149},
  {  0.09505199, -0.37838155,  0.56666351, -0.37838155,  0.09505199}
};
double chebyCoefsA[32][5] = {

  {  1.00000000,  0.00000000,  0.00000000,  0.00000000,  0.00000000},
  {  1.00000000, -2.18401672,  2.06596957, -0.89798564,  0.16469484},
  {  1.00000000, -2.79930569,  3.08123529, -1.54853970,  0.30261806},
  {  1.00000000, -3.10256029,  3.69606227, -1.99009736,  0.40936729},
  {  1.00000000, -3.28342935,  4.10065719, -2.30170922,  0.49010768},
  {  1.00000000, -3.40361345,  4.38563025, -2.53148857,  0.55232521},
  {  1.00000000, -3.48927737,  4.59676762, -2.70730851,  0.60141637},
  {  1.00000000, -3.55342822,  4.75931407, -2.84593312,  0.64101084},
  {  1.00000000, -3.60326600,  4.88824356, -2.95792561,  0.67356313},
  {  1.00000000, -3.64309984,  4.99297419, -3.05023292,  0.70076944},
  {  1.00000000, -3.67566681,  5.07971690, -3.12759630,  0.72383104},
  {  1.00000000, -3.70278895,  5.15272970, -3.19335614,  0.74361903},
  {  1.00000000, -3.72572620,  5.21502760, -3.24993079,  0.76077882},
  {  1.00000000, -3.74537762,  5.26880449, -3.29911278,  0.77579802},
  {  1.00000000, -3.76240196,  5.31569445, -3.34225820,  0.78905146},
  {  1.00000000, -3.77729293,  5.35693933, -3.38041146,  0.80083177},
  {  1.00000000, -3.79042784,  5.39349965, -3.41438976,  0.81137062},
  {  1.00000000, -3.80209998,  5.42612993, -3.44484153,  0.82085372},
  {  1.00000000, -3.81254078,  5.45543102, -3.47228778,  0.82943165},
  {  1.00000000, -3.82193531,  5.48188735, -3.49715187,  0.83722774},
  {  1.00000000, -3.83043330,  5.50589376, -3.51978134,  0.84434394},
  {  1.00000000, -3.83815724,  5.52777528, -3.54046410,  0.85086527},
  {  1.00000000, -3.84520827,  5.54780194, -3.55944071,  0.85686318},
  {  1.00000000, -3.85167064,  5.56619987, -3.57691371,  0.86239812},
  {  1.00000000, -3.85761510,  5.58315993, -3.59305483,  0.86752163},
  {  1.00000000, -3.86310152,  5.59884426, -3.60801063,  0.87227786},
  {  1.00000000, -3.86818086,  5.61339153, -3.62190694,  0.87670487},
  {  1.00000000, -3.87289681,  5.62692102, -3.63485242,  0.88083565},
  {  1.00000000, -3.87728701,  5.63953589, -3.64694135,  0.88469895},
  {  1.00000000, -3.88138408,  5.65132584, -3.65825597,  0.88831989},
  {  1.00000000, -3.88521644,  5.66236918, -3.66886833,  0.89172057},
  {  1.00000000, -3.88880893,  5.67273465, -3.67884179,  0.89492046}
};

//volatile bool keepRunning = true;

// scale streamed data to microvolts
static const double fixedGain = 1.0/8192.0; 

int64_t maxFileSize = 1024*1024*1024; // 1GB for testing
// lock-free FIFO for saving thread
int32_t *rbInt;
int rbIntSize, rbIntChans;
int rbIntWritePos;
int rbIntReadPos;
// pipe or socketpair for inter-thread communication
LocalPipe locPipe;
// GDF writing object
GDF_Writer *gdfWriter = NULL;
// Name of file to write to
char baseFilename[1024];
char curFilename[1024];
int fileCounter = 0;


bool writeHeader(int ftSocket, float fSample, const ChannelSelection &streamSel) {
	FtBufferRequest req;
	FtBufferResponse resp;
	char *chunk_data;
	int N=0,P;

	for (int n=0;n<streamSel.getSize();n++) {
		int Ln = strlen(streamSel.getLabel(n))+1;
		N+=Ln;
	}
	chunk_data = new char[N];
	
	P=0;
	for (int n=0;n<streamSel.getSize();n++) {
		int Ln = strlen(streamSel.getLabel(n))+1;
		memcpy(chunk_data + P, streamSel.getLabel(n), Ln);
		P+=Ln;
	}

	req.prepPutHeader(streamSel.getSize(), DATATYPE_FLOAT32, fSample);
	req.prepPutHeaderAddChunk(FT_CHUNK_CHANNEL_NAMES, N, chunk_data);
	
	delete[] chunk_data;
	
	int err = clientrequest(ftSocket, req.out(), resp.in());
	if (err || !resp.checkPut()) {
		fprintf(stderr, "Could not write header to FieldTrip buffer\n.");
		return false;
	}
	return true;
}

void *savingThreadFunction(void *arg) {
	int writePtr,n;
	int64_t fileSize = 256*(1+rbIntChans);

	while (1) {
		int newSamplesA, newSamplesB;
		int32_t *rbIntPtr;
		int64_t newSize;
		
		n = locPipe.read(sizeof(int), &writePtr);
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
			newSize = 256*(1+rbIntChans) + (newSamplesA+newSamplesB) * rbIntChans * sizeof(int32_t);
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

/*
void abortHandler(int sig) {
	printf("Ctrl-C pressed -- stopping...\n");
	keepRunning = false;
}
*/


int main(int argc, char *argv[]) {
	int port, ftSocket;
	int sampleCounter = 0;
	char hostname[256];
	BioSemiClient BS;
	FtBufferResponse resp;
	FtEventList   eventChain;
	FtSampleBlock sampleBlock(DATATYPE_FLOAT32);
	SignalConfiguration signalConf;
	int triggerState = 0;
	ft_buffer_server_t *ftServer = NULL;
	pthread_t savingThread;
	float *auxVec;
	int skipSamples = 0;
	MultiChannelFilter<float> *lpFilter = NULL;
	ConsoleInput ConIn;
	int acqState = 0; // 0 = stopped, 1 = running, 2 = paused
	static char *stateDescr[3] = {"INACTIVE", "RUNNING ", "PAUSED  "};
	int lenBaseFilename, sessionCount=0;
	int nsBattery = 0, nsCMS = 0;
	int cmsInRange = 0;
	
	if (argc<3) {
		printf("Usage: biosemi2ft <config-file> <gdf-file> [hostname=localhost [port=1972]]\n");
		return 0;
	}
	
	int err = signalConf.parseFile(argv[1]);
	if (err == -1) {
		fprintf(stderr, "Could not read configuration file %s\n", argv[1]);
		return 1;
	}
	if (err > 0) {
		fprintf(stderr, "Encountered %i errors in configuration file - aborting\n", err);
		return 1;
	}
	
	// enable for debugging if you like
	if (0) {
		const ChannelSelection& streamSel = signalConf.getStreamingSelection();
		const ChannelSelection& saveSel = signalConf.getSavingSelection();
		
		printf("Selected for streaming:\n");
		for (int i=0;i<streamSel.getSize();i++) {
			printf("%3i: %s\n", streamSel.getIndex(i)+1, streamSel.getLabel(i));
		}
		printf("Selected for saving:\n");
		for (int i=0;i<saveSel.getSize();i++) {
			printf("%3i: %s\n", saveSel.getIndex(i)+1, saveSel.getLabel(i));
		}
	}
	
	strcpy(baseFilename, argv[2]);
	char *lastDotPos = strrchr(baseFilename, '.');
	// cut off .gdf suffix, if given
	if (lastDotPos != NULL) {
		if (!strcasecmp(lastDotPos, ".gdf")) {
			*lastDotPos = 0;
		}
	}
	lenBaseFilename = strlen(baseFilename);
	
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
		ftSocket = 0;
	}
	
	if (!BS.openDevice()) return 1;
	
	double T0 = BS.getCurrentTime();
	
	// the number of HW channels (excluding sync + status)
	int numChanHW = BS.getNumChannels() + BS.getNumChanAIB();
	
	auxVec = new float[numChanHW];
	
	printf("Speed mode.........: %i\n", BS.getSpeedMode());
	printf("Number of channels.: %i\n", BS.getNumChannels());
	printf("Number of AIB chns.: %i\n", BS.getNumChanAIB());
	printf("Sampling frequency.: %i Hz\n", BS.getSamplingFreq());
	
	if (numChanHW < signalConf.getMaxSavingChannel() || numChanHW < signalConf.getMaxStreamingChannel() ) {
		fprintf(stderr, "Configuration file specifies channel numbers beyond hardware limits.\n");
		return 1;
	}
	
	const ChannelSelection& streamSel = signalConf.getStreamingSelection();
	const ChannelSelection& saveSel = signalConf.getSavingSelection();
	
	printf("Selected %i channels for streaming, %i channels for saving to GDF.\n", streamSel.getSize(), saveSel.getSize());
	
	if (streamSel.getSize() > 0) {
		if (signalConf.getBandwidth() <= 0) {
			lpFilter = new MultiChannelFilter<float>(streamSel.getSize(), 4);
			int coefInd = signalConf.getDownsampling() - 1;
			if (coefInd > 31) {
				fprintf(stderr, "Warning, no filter coefficients for more than 32x downsampling\n");
				coefInd = 31;
			}
			lpFilter->setCoefficients(chebyCoefsB[coefInd], chebyCoefsA[coefInd]);
		} else {
			lpFilter = new MultiChannelFilter<float>(streamSel.getSize(), signalConf.getOrder());
			lpFilter->setButterLP(signalConf.getBandwidth() / (0.5*BS.getSamplingFreq()));
		}
	}
	
	if (saveSel.getSize() > 0) {
		gdfWriter = new GDF_Writer(1+saveSel.getSize(), BS.getSamplingFreq(), GDF_INT32);
		gdfWriter->setLabel(0, "STATUS");
		for (int i=0;i<saveSel.getSize();i++) {
			//int idx = saveSel.getLabel(i);
			gdfWriter->setLabel(1+i, saveSel.getLabel(i));
			// min:  (-2^31)/8192.0  (int32 -> microVolt)
			// max: (2^31-1)/8192.0  (int32 -> microVolt)
			gdfWriter->setPhysicalLimits(1+i, -262144.0, 262143.99987792969);
			gdfWriter->setPhysDimCode(1+i, GDF_MICRO + GDF_VOLT);
		}
	}
	
	rbIntSize = 5*BS.getSamplingFreq();	// enough space for 5 seconds of data
	rbIntChans = 1+saveSel.getSize();
	rbInt = new int32_t[rbIntChans*rbIntSize];
	rbIntWritePos = 0;
	rbIntReadPos = 0;
	for (int i=0;i<rbIntSize;i++) rbInt[i*rbIntChans]=-1;
	
	if (pthread_create(&savingThread, NULL, savingThreadFunction, NULL)) {
		fprintf(stderr, "Could not spawn GDF saving thread.\n");
		return 1;
	}
		
	/* register CTRL-C handler */
	//signal(SIGINT, abortHandler);
	//printf("Starting to listen - press CTRL-C to quit\n");
	
	printf("\nPress <S> to start/stop acquisition, <P> to pause/unpause, <Esc> to quit\n\n");
	
	//while (keepRunning) {
	while (1) {
		BioSemiBlock block;
		bool newBlock;
		int nStream = (ftSocket >= 0) ? streamSel.getSize() : 0;
		int nSave = signalConf.getSavingSelection().getSize();
		int deci = signalConf.getDownsampling();
		
		if (ConIn.checkKey()) {
			int c = ConIn.getKey();
			if (c==27) break; // quit
			
			if (c=='s' || c=='S') {
				if (acqState == 0) {
					// start acquisition
					++sessionCount;
					fileCounter = 0;
					sprintf(baseFilename + lenBaseFilename, "_S%i", sessionCount);
					strcpy(curFilename, baseFilename);
					strcat(curFilename, ".gdf");
					if (!gdfWriter->createAndWriteHeader(curFilename)) {
						fprintf(stderr, "Could not open GDF file %s for writing\n", curFilename);
						break;
					}
					if (!writeHeader(ftSocket, BS.getSamplingFreq() / signalConf.getDownsampling(), streamSel)) {
						break;
					}
					nsBattery = nsCMS = 0;
					acqState = 1;
				} else {
					// stop acquisition: first wait for writing thread to finish
					while (rbIntReadPos != rbIntWritePos) {
						BS.msleep(1);
					}
					// then close the file
					gdfWriter->close();
					acqState = 0;
				}
			} else if (c=='p' || c=='P') {
				if (acqState == 1) {
					acqState = 2;
				} else if (acqState == 2) {
					acqState = 1;
				}
			}
		}
		
		newBlock = BS.checkNewBlock(block);
		if (!newBlock) {
			BS.msleep(1);
			continue;
		}
		
		double batt = BS.getValue(block.startIndex + 279);
		batt=(batt/2097152)-175.5;
		if (batt < 0) batt = 0;
		if (batt > 100) batt = 100;
		
		printf("%s T=%8.3f Ptr=%8i, samples=%3i, batt=%5.1f%%\r", stateDescr[acqState], BS.getCurrentTime() - T0, block.startIndex, block.numSamples, batt);
		if (block.numSamples != block.numInSync) {
			fprintf(stderr, "USB device out of sync (%i / %i) -- exiting\n", block.numInSync, block.numSamples);
			break;
		}
		
		eventChain.clear();
		if (nsBattery >= 0) {
			if (nsBattery<block.numSamples) {
				int dSamC = (sampleCounter + deci - 1)/deci;

				// we're due a report of the battery level
				eventChain.add(dSamC, "BATTERY", (float) batt);
				//printf("\n-!- Battery level = %4.1f\n", batt);
				// next report in getBatteryRefresh() seconds
				nsBattery = signalConf.getBatteryRefresh() * BS.getSamplingFreq();
				// if this is zero, user doesn't want events,
				// so prevent further output by setting it to -1
				if (nsBattery == 0) nsBattery = -1; 
			} else {
				nsBattery -= block.numSamples;
			}
		}
		
		for (int j=0;j<block.numSamples;j++) {
			int value  = BS.getValue(block.startIndex + 1 + j*block.stride);
			int cmsBit = (value & 0x10000000) ? 0:1;
			int dSamC = (sampleCounter + j + deci - 1)/deci;
			
			if (nsCMS == 0 || cmsBit != cmsInRange) {
				eventChain.add(dSamC, "CMS_IN_RANGE", cmsBit);
				//printf("\n-!- CMS in range: %i\n", cmsBit);
				cmsInRange = cmsBit;
				nsCMS = signalConf.getStatusRefresh() * BS.getSamplingFreq();
			} else {
				--nsCMS;
			}
			
			value = (value & 0x00FFFF00) >> 8;
			if (value && value!=triggerState) {
				if (!signalConf.useSplittedTrigger())  {
					eventChain.add(dSamC, "TRIGGER", value);
					printf("\n-!- Trigger at sample %i => %i\n", dSamC, value);
				} else {
					if ((value & 0xFF00) != (triggerState & 0xFF00)) {
						eventChain.add(dSamC, signalConf.getHighTriggerName(), value >> 8);
						printf("\n-!- %s at sample %i => %i\n", signalConf.getHighTriggerName(), dSamC, value >> 8);
					}
					if ((value & 0x00FF) != (triggerState & 0x00FF)) {
						eventChain.add(dSamC, signalConf.getLowTriggerName(), value & 0x00FF);
						printf("\n-!- %s at sample %i => %i\n", signalConf.getLowTriggerName(), dSamC, value & 0x00FF);
					}
				}
			}
			triggerState = value;
		}
		
		if (acqState != 1) continue;
		
		if (nSave > 0) {
			for (int j=0;j<block.numSamples;j++) {
				int doff = rbIntChans*rbIntWritePos;
				int soff = block.startIndex + j*block.stride;
				
				if (rbInt[doff] != -1) {
					fprintf(stderr, "Error: saving thread does not keep up with load\n");
					break;
				}
				
				rbInt[doff++] = BS.getValue(soff + 1);
				for (int i=0;i<nSave;i++) {
					int ihw = 2 + streamSel.getIndex(i);
					rbInt[doff+i] = BS.getValue(soff + ihw);
				}
				if (++rbIntWritePos == rbIntSize) rbIntWritePos = 0;
			}
			locPipe.write(sizeof(int), static_cast<void *>(&rbIntWritePos)); 
		}
		
		if (nStream > 0) {
			int deci = signalConf.getDownsampling();
			int numThisTime = (block.numSamples - skipSamples + deci - 1)/deci;
		
			float *dest = (float *) sampleBlock.getMatrix(nStream, numThisTime);
			if (dest==NULL) {
				printf("Out of memory\n");
				break;
			}
			for (int j=0;j<block.numSamples;j++) {
				int soff = block.startIndex + 2 + j*block.stride;
				for (int i=0;i<nStream;i++) {
					int ihw = streamSel.getIndex(i);
					auxVec[i] = fixedGain * BS.getValue(soff + ihw);
				}
				
				if (skipSamples == 0) {
					lpFilter->process(dest, auxVec);
					dest += nStream;
				} else {
					lpFilter->process(auxVec);
				}
				
				if (--skipSamples < 0) skipSamples = deci-1;
			}
			int err = clientrequest(ftSocket, sampleBlock.asRequest(), resp.in());
			if (err || !resp.checkPut()) {
				fprintf(stderr, "Could not write samples to FieldTrip buffer\n.");
			}
		}

		if (ftSocket != -1 && eventChain.count() > 0) {
			int err = clientrequest(ftSocket, eventChain.asRequest(), resp.in());
			if (err || !resp.checkPut()) {
				fprintf(stderr, "Could not write events to FieldTrip buffer\n.");
			}
		}
		
		sampleCounter += block.numSamples;
	}
	
	skipSamples = -1;
	locPipe.write(sizeof(int), &skipSamples);
	
	BS.closeDevice();
	if (ftSocket > 0) close_connection(ftSocket);
	if (ftServer != NULL) ft_stop_buffer_server(ftServer);
	if (gdfWriter) {
		gdfWriter->close();
		delete gdfWriter;
	}
	delete[] auxVec;
	
	return 0;
}
