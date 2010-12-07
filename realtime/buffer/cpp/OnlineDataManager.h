#include <GDF_BackgroundWriter.h>
#include <SignalConfiguration.h>
#include <stdio.h>
#include <FtBuffer.h>
#include <socketserver.h>
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
		
		curWriter = 0;
		lpFilter = 0;
				
		savingEnabled = false;
	}
	
	virtual ~OnlineDataManager() {
		// TODO: some more of this
	
		// stop buffer server, if spawned
		if (ftServer) ft_stop_buffer_server(ftServer);
		// FtConnection is cleaned up automatically, if needed
		
		// clean up GDF writer, if needed
		if (curWriter) {
			curWriter->stopSync();
			delete curWriter;
		}
		// clean up variables
		delete lpFilter;
		delete[] auxVec;
		delete[] pBlock;
	}
	
	virtual std::string handleStringRequest(const std::string& request) {
		static const std::string ok("OK\n");
		static const std::string unknown("ERROR: UNKNOWN COMMAND\n");
		static const std::string malform("ERROR: MALFORMED COMMAND\n");
		static const std::string stopFirst("ERROR: STOP FIRST\n");
		static const std::string emptyFilename("ERROR: NO FILENAME SET\n");
		static const std::string chanOutOfRange("ERROR: CHANNEL INDEX OUTSIDE HARDWARE LIMITS\n");
		static const std::string cannotSave("ERROR: COULD NOT START SAVING\n");
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
				enableStreaming();
			} else {
				if (savingEnabled) return ok; // silently ignore that it's already running
				if (curFilename.empty()) return emptyFilename;
				if (!enableSaving()) return cannotSave;
			}
			return ok;
		} else if (token2.compare("STOP") == 0) {
			if (!StringServer::getNextToken(request, pos).empty()) return malform;
			if (target == 1) {
				if (!streamingEnabled) return ok; // silently ignore that it's already stopped
				disableStreaming();
			
			} else {
				if (!savingEnabled) return ok; // silently ignore that it's already stopped
				disableSaving();
			}
			return ok;
		} else if (token2.compare("SELECT") == 0) {
			if ((target==1 && streamingEnabled) || (target==2 && savingEnabled)) {
				return stopFirst;
			}
			ChannelSelection cs;
			if (!cs.parseString(request.size() - pos, request.data() + pos)) return malform;
			if (cs.getMaxIndex() >= nCont) return chanOutOfRange;
						
			if (target==1) {
				signalConf.setStreamingSelection(cs);
			} else {
				signalConf.setSavingSelection(cs);
			}
			return ok;
		} else if (target == 2 && token2.compare("FILE") == 0) {
			if (savingEnabled) return stopFirst;
			
			std::string filename = StringServer::getNextToken(request, pos);
			if (!StringServer::getNextToken(request, pos).empty()) return malform;
			
			setFilename(filename);
			return ok;
		} else if (target == 1 && token2.compare("FILTER") == 0) {
			if (streamingEnabled) return stopFirst;
			
			double bandwidth = 0.0;
			int order = 0;
			int factor = 0;
						
			if (convertToDouble(StringServer::getNextToken(request, pos), bandwidth) 
					&& (bandwidth >= 0)
					&& convertToInt(StringServer::getNextToken(request, pos), order)
					&& (order >= 0)
					&& convertToInt(StringServer::getNextToken(request, pos), factor)
					&& (factor >= 1)
					&& StringServer::getNextToken(request, pos).empty()) {
					
				signalConf.setDownsampling(factor);
				if (bandwidth < 0.5*fSample) {
					signalConf.setBandwidth(bandwidth);
					signalConf.setOrder(order);
				} else {
					signalConf.setOrder(0);
				}
				configureStreaming();
				
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
			configureStreaming();
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
		if (streamingEnabled) {
			if (!handleStreaming()) return false;
		}
		if (savingEnabled) {
			return handleSaving();
		}
		return true;
	}
	
	bool enableSaving() {
		if (curFilename.empty()) return false;
		if (curWriter == 0) {
			configureSaving();
			if (curWriter == 0) return false;
			printf("Reconfigured saving\n");
		}
		if (fileCounter == 0) {
			curWriter->start(curFilename.c_str());
		} else {
			char *compName = new char[curFilename.size() + 20];
			sprintf(compName, "%s_S%i", curFilename.c_str(), fileCounter);
			printf("COMPOSED: %s\n", compName);
			curWriter->start(compName);
			delete[] compName;
		}
		fileCounter++;		
		savingEnabled = true;
		return true;
	}
	
	void disableSaving() {
		savingEnabled = false;
		if (!curWriter) return;
		curWriter->stopAsync();
		curWriter = 0;
	}
	
	bool enableStreaming() {
		if (streamingEnabled) return true; // silently ignore      
		if (!writeHeader()) return false;
		streamingEnabled = true;
		return true;
	}
	
	void disableStreaming() {
		// if (!streamingEnabled) return;
		streamingEnabled = false;
	}
	
	FtEventList& getEventList() {
		return eventList;
	}
   
	GDF_Writer *getCurrentGDF() {
		if (curWriter) return curWriter->gdf();
		return 0;
	}
	
	void setFilename(const std::string& filename) {
		curFilename = filename;
		fileCounter = 0;
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
		
		if (curWriter == 0) return false;
		
		if (!curWriter->checkFreeBlock(nThisBlock)) {
			fprintf(stderr, "Error: saving thread does not keep up with load\n");
			return false;
		}
		
		for (int j=0;j<nThisBlock;j++) {
			To *dest = curWriter->getSampleSlot();
			To *src  = pBlock  + j*stride;
				
			for (int i=0;i<nStatus;i++) {
				*dest++ = *src++;
			}
			for (int i=0;i<nSave;i++) {
				dest[i] = src[saveSel.getIndex(i)];
			}
		}
		curWriter->commitBlock();
		return true;
	}
	
	void configureSaving() {
		const ChannelSelection& saveSel = signalConf.getSavingSelection();
		
		int nSave = saveSel.getSize();
		if (nStatus + nSave == 0) return;
		
		curWriter = new GDF_BackgroundWriter<To>(nStatus + nSave, (int) fSample, gdfType);
		for (int i=0;i<nStatus;i++) {
			char lab[16];
			sprintf(lab, "Status_%i", i+1);
			curWriter->gdf().setLabel(i, lab);
		}
		for (int i=0;i<nSave;i++) {
			curWriter->gdf().setLabel(nStatus+i, saveSel.getLabel(i));
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
		
	////////////////////////////////////////////////////////////////
	// Class members
	////////////////////////////////////////////////////////////////
	GDF_BackgroundWriter<To> *curWriter;	/**< currently active GDF writer (in background thread) */
	MultiChannelFilter<Ts,Ts> *lpFilter;	/**< currently active low-pass filter */

	int ftType;			/**< FieldTrip buffer data type */
	GDF_Type gdfType;	/**< GDF data type */
	int nStatus, nCont; /**< Number of status + continuously sampled channels */
	float fSample;		/**< Sampling rate */

	int nThisBlock;		/**< Number of samples in this block (as requested by provideBlock) */
	int allocSizeBlock; /**< Size of buffer that is allocated for providing blocks */
	To *pBlock;			/**< Points to buffer that is allocated for providing blocks */
	Ts *auxVec;			/**< Big enough for keeping one sample of data (streaming format) */
	
	int ftSocket;		/**< The FT buffer socket identifier or 0 for dmarequests, -1 for none */
	int sampleCounter;	/**< Number of samples streamed out since last writeHeader */
	int skipSamples;	/**< Helper variable to keep track of downsampling operation */
	
	FtConnection ftConnection;
	FtBufferResponse resp;
	FtEventList eventList;
	FtSampleBlock sampleBlock;
	ft_buffer_server_t *ftServer;
	
	SignalConfiguration signalConf;
	
	bool savingEnabled, streamingEnabled;
	
	std::string curFilename;
	int fileCounter;
};
