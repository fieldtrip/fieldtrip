#include <GdfWriter.h>
#include <SignalConfiguration.h>
#include <stdio.h>
#include <FtBuffer.h>
#include <socketserver.h>
#include <pthread.h>
#include <LocalPipe.h>
#include <MultiChannelFilter.h>
#include <StringServer.h>

bool convertToInt(const std::string& in, int& value) {
	const char *start = in.c_str();
	char *end;
	long v = strtol(start, &end, 10);
	if (start == end) return false;
	if (*end!=0) return false;
	value = (int) v;
	return true;
} 

bool convertToDouble(const std::string& in, double& value) {
	const char *start = in.c_str();
	char *end;
	double v = strtod(start, &end);
	if (start == end) return false;
	if (*end!=0) return false;
	value = v;
	return true;
} 


template <typename To, typename Ts> 
class OnlineDataManager : public StringRequestHandler {
	public:
	
	OnlineDataManager(int nStatus, int nCont, float fSample, GDF_Type gdfType, int ftType) : sampleBlock(ftType) {
		this->nStatus = nStatus;
		this->nCont   = nCont;
		this->gdfType = gdfType;
		this->ftType  = ftType;
		this->fSample = fSample;
		
		pBlock = 0;
		allocSizeBlock = 0;
		auxVec = new Ts[nCont];
		
		ftSocket = -1;
		ftServer = 0;
		
		sampleCounter = 0;
		skipSamples = 0;
		
		gdfWriter = 0;
		lpFilter = 0;
		
		rbWritePos = rbReadPos = 0;
		rbChans = rbSize = 0;
		rbData = 0;
		rbMaxSize = (int) (5*fSample); // 5 seconds of data for saving ring buffer
		
		savingEnabled = false;
		isValid = false;
		
		maxFileSize = 1024*1024*1024; // 1GB
		fileCounter = 0;
	}
	
	virtual ~OnlineDataManager() {
		// stop saving thread
		int64_t minusOne = -1;
		locPipe.write(sizeof(int64_t), &minusOne);
		// TODO: some more of this
	
		// stop buffer server, if spawned
		if (ftServer) ft_stop_buffer_server(ftServer);
		// FtConnection is cleaned up automatically, if needed
		
		// clean up GDF writer, if needed
		if (gdfWriter) {
			gdfWriter->close();
			delete gdfWriter;
		}
		// clean up variables
		delete lpFilter;
		delete[] auxVec;
		delete[] pBlock;
		delete[] rbData;
	}
	
	virtual std::string handleStringRequest(const std::string& request) {
		static const std::string ok("OK\n");
		static const std::string unknown("ERROR: UNKNOWN COMMAND\n");
		static const std::string malform("ERROR: MALFORMED COMMAND\n");
		static const std::string stopFirst("ERROR: STOP FIRST\n");
		int target = 0; // 1 = STREAM, 2= SAVE
		unsigned int pos = 0;
		
		std::string token1 = StringServer::getNextToken(request, pos);
		if (token1.empty()) return ok; // empty command -> nothing to do
		
		if (token1.compare("STREAM") == 0) {
			target = 1;
		} else if (token1.compare("SAVE") == 0) {
			target = 2;
		} else {
			return unknown;
		}
		
		std::string token2 = StringServer::getNextToken(request, pos);
		
		if (token2.compare("START") == 0) {
			if (!StringServer::getNextToken(request, pos).empty()) return malform;
			if (target == 1) {
				if (streamingEnabled) return ok; // silently ignore that it's already running
			
			} else {
				if (savingEnabled) return ok; // silently ignore that it's already running
			
			}
			return ok;
		} else if (token2.compare("STOP") == 0) {
			if (!StringServer::getNextToken(request, pos).empty()) return malform;
			if (target == 1) {
				if (!streamingEnabled) return ok; // silently ignore that it's already stopped
			
			} else {
				if (!savingEnabled) return ok; // silently ignore that it's already stopped
			
			}
			return ok;
		} else if (token2.compare("SELECT") == 0) {
			if ((target==1 && streamingEnabled) || (target==2 && savingEnabled)) {
				return stopFirst;
			}
			ChannelSelection cs;
			if (!cs.parseString(request.size() - pos, request.data() + pos)) return malform;

			printf("New selection:\n");
			for (int i=0;i<cs.getSize();i++) {
				printf("%i : '%s'\n", cs.index[i], cs.label[i].c_str());
			}
			return ok;

		} else if (target == 2 && token2.compare("FILE") == 0) {
			if (savingEnabled) return stopFirst;
			
			std::string filename = StringServer::getNextToken(request, pos);
			if (!StringServer::getNextToken(request, pos).empty()) return malform;
			
			printf("new filename: %s\n", filename.c_str());
			return ok;
		} else if (target == 1 && token2.compare("FILTER") == 0) {
			if (streamingEnabled) return stopFirst;
			
			double bandwidth;
			int order;
			int factor;
						
			if (convertToDouble(StringServer::getNextToken(request, pos), bandwidth) 
					&& convertToInt(StringServer::getNextToken(request, pos), order)
					&& convertToInt(StringServer::getNextToken(request, pos), factor)
					&& StringServer::getNextToken(request, pos).empty()) {
							
				printf("stream filter %f , %i , %i\n", bandwidth, order, factor);
				return ok;
			} else {
				return malform;
			}
		} 
		
		return unknown;
	}
	
	bool connectToServer(char *address) {
		if (ftSocket != -1) return false;
		
		if (!ftConnection.connect(address)) return false;
		ftSocket = ftConnection.getSocket();
		return true;
	}
	
	bool useOwnServer(int port) {
		if (ftSocket != -1) return false;
		
		ftServer = ft_start_buffer_server(port, NULL, NULL, NULL);
		if (ftServer == 0) return false;
		
		ftSocket = 0; // => dma
		return true;
	}
	
	int configureFromFile(const char *filename) {
		int numErr = signalConf.parseFile(filename);
		if (numErr == 0) {
			configureSaving();
			configureStreaming();
			isValid = true;
		} else {
			isValid = false;
		}
		return numErr;
	}
	
	To *provideBlock(int N) {
		int needed = N * (nStatus + nCont);
		if (needed > allocSizeBlock) {
			delete[] pBlock;
			pBlock = new To[needed];
			allocSizeBlock = needed;
		}
		nThisBlock = N;
		
		eventList.clear();
		return pBlock;
	}
	
	bool handleBlock() {
		if (!isRunning) return false;
		
		if (ftSocket != -1) {
			if (!handleStreaming()) return false;
		}
		
		if (savingEnabled) {
			return handleSaving();
		}
		return true;
	}
	
	bool enableSaving(const char *filename) {
		strcpy(baseFilename, filename);
		char *lastDotPos = strrchr(baseFilename, '.');
		// cut off .gdf suffix, if given
		if (lastDotPos != NULL) {
			if (!strcasecmp(lastDotPos, ".gdf")) {
				*lastDotPos = 0;
			}
		}
		
		fileCounter = 0;
		strcpy(curFilename, baseFilename);
		strcat(curFilename, ".gdf");
		if (!gdfWriter->createAndWriteHeader(curFilename)) {
			fprintf(stderr, "Could not open GDF file %s for writing\n", curFilename);
			return false;
		}
		
		if (rbData) delete[] rbData;
		rbChans = nStatus + signalConf.getSavingSelection().getSize();
		rbData = new To[rbMaxSize*rbChans];
		rbSize = rbMaxSize;
		rbWritePos = rbReadPos = 0;
		
		savingEnabled = true;
		return true;
	}
	
	bool start() {
		if (pthread_create(&savingThread, NULL, savingThreadFunction, this)) {
			fprintf(stderr, "Could not spawn GDF saving thread.\n");
			return false;
		}
	
		if (writeHeader()) {
			return false;
		}
		isRunning = true;
		return true;
	}
	
	void stop() {
		if (!isRunning) return;
		isRunning = false;
	}
	
	FtEventList& getEventList() {
		return eventList;
	}
		
	protected:
	
	bool writeHeader() {
		const ChannelSelection& streamSel = signalConf.getStreamingSelection();
		FtBufferRequest req;
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

		req.prepPutHeader(streamSel.getSize(), DATATYPE_FLOAT32, fSample / signalConf.getDownsampling());
		req.prepPutHeaderAddChunk(FT_CHUNK_CHANNEL_NAMES, N, chunk_data);
	
		delete[] chunk_data;
	
		int err = clientrequest(ftSocket, req.out(), resp.in());
		if (err || !resp.checkPut()) {
			fprintf(stderr, "Could not write header to FieldTrip buffer\n.");
			return false;
		}
		return true;
	}
	
	bool handleStreaming() {
		int stride = nStatus + nCont;
		int err;
		const ChannelSelection& streamSel = signalConf.getStreamingSelection();
		int nStream  = streamSel.getSize();
			
		// write events, if any
		if (eventList.count() > 0) {
			err = clientrequest(ftSocket, eventList.asRequest(), resp.in());
			if (err || !resp.checkPut()) {
				fprintf(stderr, "Could not write events to FieldTrip buffer.\n");
				return false;
			}
		}
			
		// write samples, if channels selected
		if (nStream == 0) return true;
		
		int deci        = signalConf.getDownsampling();
		int numThisTime = (nThisBlock - skipSamples + deci - 1)/deci;
	
		Ts *dest = (Ts *) sampleBlock.getMatrix(nStream, numThisTime);
			
		if (lpFilter != 0) {
			for (int j=0;j<nThisBlock;j++) {
				To *src = pBlock + nStatus + j*stride;
				for (int i=0;i<nStream;i++) {
					auxVec[i] = src[streamSel.getIndex(i)];
				}
			
				if (skipSamples == 0) {
					lpFilter->process(dest, auxVec);
					dest += nStream;
				} else {
					lpFilter->process(auxVec);
				}
				if (--skipSamples < 0) skipSamples = deci-1;
			}
		} else {
			for (int j=0;j<nThisBlock;j++) {
				To *src = pBlock + nStatus + j*stride;
				if (skipSamples == 0) {
					for (int i=0;i<nStream;i++) {
						dest[i] = src[streamSel.getIndex(i)];
					}
					dest += nStream;
				}
				if (--skipSamples < 0) skipSamples = deci-1;
			}
		}
		err = clientrequest(ftSocket, sampleBlock.asRequest(), resp.in());
		if (err || !resp.checkPut()) {
			fprintf(stderr, "Could not write samples to FieldTrip buffer.\n");
			return false;
		}
		return true;
	}
	
	bool handleSaving() {
		const ChannelSelection& saveSel = signalConf.getSavingSelection();
		int nSave  = saveSel.getSize();
		int stride = nStatus + nCont;
		
		if (rbWritePos - rbReadPos > rbSize - nThisBlock) {
			fprintf(stderr, "Error: saving thread does not keep up with load\n");
			return false;
		}
		
		int wrappedPos = rbWritePos % rbSize;
		
		for (int j=0;j<nThisBlock;j++) {
			To *dest = rbData + rbChans*wrappedPos;
			To *src  = pBlock  + j*stride;
				
			for (int i=0;i<nStatus;i++) {
				*dest++ = *src++;
			}
			for (int i=0;i<nSave;i++) {
				dest[i] = src[saveSel.getIndex(i)];
			}
			if (++wrappedPos == rbSize) wrappedPos = 0;
		}
		rbWritePos += nThisBlock;
		locPipe.write(sizeof(int64_t), static_cast<void *>(&rbWritePos));
		return true;
	}
	
	void configureSaving() {
		const ChannelSelection& saveSel = signalConf.getSavingSelection();
		
		delete gdfWriter;
		
		int nSave = saveSel.getSize();
		if (nStatus + nSave == 0) return;
		
		gdfWriter = new GDF_Writer(nStatus + nSave, (int) fSample, gdfType);
		for (int i=0;i<nStatus;i++) {
			char lab[16];
			sprintf(lab, "Status_%i", i+1);
			gdfWriter->setLabel(i, lab);
		}
		for (int i=0;i<nSave;i++) {
			gdfWriter->setLabel(nStatus+i, saveSel.getLabel(i));
			// gdfWriter->setPhysicalLimits(1+i, -262144.0, 262143.99987792969);
			// gdfWriter->setPhysDimCode(1+i, GDF_MICRO + GDF_VOLT);
		}
	}
	
	void configureStreaming() {
		delete lpFilter;
		
		if (signalConf.getOrder() > 0) {
			lpFilter = new MultiChannelFilter<Ts,Ts>(signalConf.getStreamingSelection().getSize(), signalConf.getOrder());
			lpFilter->setButterLP(signalConf.getBandwidth() / (0.5*fSample));
		} else {
			lpFilter = 0;
		}
	}
	
	void savingThreadFunc() {
		int64_t fileSize = 256*(1+rbChans);

		while (1) {
			int newSamplesA, newSamplesB;
			To *rbPtr;
			int64_t newSize, newWritePos;
		
			int n = locPipe.read(sizeof(int64_t), &newWritePos);
			if (n!=sizeof(int64_t)) {
				// this should never happen for blocking sockets/pipes
				fprintf(stderr, "Unexpected error in pipe communication\n");
				break;
			}
			if (newWritePos < 0) {
				printf("\nSaving thread received %i - exiting...\n", (int) newWritePos);
				break;
			}
			
			int writePtr = newWritePos % rbSize;
			int readPtr  = rbReadPos % rbSize;
			
			printf("Saving %8lli [%i; %i(\n", newWritePos, readPtr, writePtr);
			
			rbPtr = rbData + readPtr*rbChans;
			if (writePtr > readPtr) {
				newSamplesA = writePtr - readPtr;
				newSamplesB = 0;
			} else {
				newSamplesA = rbSize - readPtr;
				newSamplesB = writePtr;
			}
			
			int64_t addSize = (newSamplesA+newSamplesB) * rbChans * sizeof(To);
		
			newSize = fileSize + addSize;
			if (newSize > maxFileSize) {
				gdfWriter->close();
				fileCounter++;
			
				snprintf(curFilename, sizeof(curFilename), "%s_%i.gdf", baseFilename, fileCounter);
				if (!gdfWriter->createAndWriteHeader(curFilename)) {
					fprintf(stderr, "Error: could not create GDF file %s\n", curFilename);
					break;
				}
				newSize = 256*(1+rbChans) + addSize;
			}
		
			gdfWriter->addSamples(newSamplesA, rbPtr);
			if (newSamplesB > 0) {
				gdfWriter->addSamples(newSamplesB, rbData);
			}
			rbReadPos = newWritePos;
			readPtr   = writePtr;
			fileSize  = newSize;
		}
	}
	
	static void *savingThreadFunction(void *arg) {
		if (arg == 0) return NULL;
		
		OnlineDataManager<To,Ts> *ODM = (OnlineDataManager<To,Ts> *) arg;
		ODM->savingThreadFunc();
		return NULL;
	}
	
	GDF_Writer *gdfWriter;
	MultiChannelFilter<Ts,Ts> *lpFilter;

	int ftType;
	GDF_Type gdfType;
	int nStatus, nCont;
	float fSample;

	int nThisBlock;
	int allocSizeBlock;
	
	To *pBlock;
	Ts *auxVec;
	
	To *rbData;
	int rbSize, rbMaxSize, rbChans;
	int64_t rbReadPos, rbWritePos;
	
	int ftSocket;
	int sampleCounter;
	
	LocalPipe locPipe;
	FtConnection ftConnection;
	FtBufferResponse resp;
	FtEventList eventList;
	FtSampleBlock sampleBlock;
	SignalConfiguration signalConf;
	ft_buffer_server_t *ftServer;
	pthread_t savingThread;
	
	int skipSamples;
	
	int64_t maxFileSize;
	char baseFilename[1024];
	char curFilename[1024];
	int fileCounter;
	
	bool savingEnabled, streamingEnabled;
	bool isValid, isRunning;
};
