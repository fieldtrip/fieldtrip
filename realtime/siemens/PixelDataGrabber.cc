#include <PixelDataGrabber.h>
#include <stdio.h>
#include <math.h>

#include <message.h>
#include <buffer.h>
#include <FolderWatcher.h>

PixelDataGrabber::PixelDataGrabber() {
	ftbSocket = -1;
	FW = NULL;
	fwEventHandle = INVALID_HANDLE_VALUE;
	
	phaseFOV = readoutFOV = 0.0;
	samplesWritten = 0;
	readResolution = 0;
	numSlices = 0;
	lastAction = Nothing;
	protInfo = NULL;
}

PixelDataGrabber::~PixelDataGrabber() {
	if (fwEventHandle != INVALID_HANDLE_VALUE) CloseHandle(fwEventHandle);
	if (FW) delete FW;
		
	if (ftbSocket > 0) close_connection(ftbSocket);
	sap_destroy(protInfo);
}
	
bool PixelDataGrabber::monitorDirectory(const char *directory) {
	if (directory == 0) return false;

	if (fwEventHandle == INVALID_HANDLE_VALUE) {
		fwEventHandle = CreateEvent(NULL, FALSE, FALSE, NULL);
		if (fwEventHandle == INVALID_HANDLE_VALUE) return false;
	}
	
	if (FW) {
		FW->stopListenForChanges();
		delete FW;
	}
	sourceDir = directory;
	int n = sourceDir.size();
	if (n>0) {
		if (sourceDir[n-1]!='\\' && sourceDir[n-1]!='/') sourceDir += '\\';

		FW = new FolderWatcher(directory);
		FW->startListenForChanges(fwEventHandle);
	}

	return FW->isValid();
}

bool PixelDataGrabber::connectToFieldTrip(const char *hostname, int port) {
	if (ftbSocket > 0) close_connection(ftbSocket);
	ftbSocket = open_connection(hostname, port);
	return (ftbSocket > 0);
}
	
int PixelDataGrabber::run(unsigned int ms) {
	if (ftbSocket <= 0) {
		fprintf(stderr, "No fieldtrip buffer has been connected to.\n");
		return -1;
	}
	if (FW==NULL || !FW->isValid()) {
		fprintf(stderr, "No directory is being monitored.\n");
		return -1;
	}

	DWORD waitResult = WaitForSingleObject(fwEventHandle, ms);

	if (waitResult == WAIT_FAILED) {
		fprintf(stderr, "Error in WaitForSingleObject(...) !\n");
		// TODO: proper error handling
		return -1;
	}
		
	if (waitResult == WAIT_TIMEOUT) return 0;

	if (FW->processChanges() > 0) {
		tryFolderToBuffer();
		FW->startListenForChanges(fwEventHandle);
		return 1;
	}
	return 0;
}

	
void PixelDataGrabber::writeHeader() {
	headerdef_t header_def;

	message_t request;
	messagedef_t request_def;
	message_t *response = NULL;
			
	header_def.nchans = (UINT32_T) numSlices*readResolution*phaseResolution;
	header_def.nsamples = 0;
	header_def.nevents = 0;
	header_def.fsample = 0.5; /* TODO: is this always correct ? */
	header_def.data_type = DATATYPE_UINT16;
	header_def.bufsize = protBuffer.size();
	
	request.def = &request_def;
	request.buf = NULL;
	request_def.version = VERSION;
	request_def.command = PUT_HDR;
	request_def.bufsize = append(&request.buf, 0, &header_def, sizeof(headerdef_t));
	request_def.bufsize = append(&request.buf, request_def.bufsize, protBuffer.data(), protBuffer.size());
	
	// write the request, read the response
	int result = tcprequest(ftbSocket, &request, &response);		
	
	if (result < 0) {
		fprintf(stderr, "Communication error when sending header to fieldtrip buffer\n");
	}
	
	if (!response || !response->def) {
		fprintf(stderr, "PUT_HDR: unknown error in response\n");
	} else {
		if (response->def->command!=PUT_OK) {
			fprintf(stderr, "PUT_HDR: Buffer returned an error (%d)\n", response->def->command);
		}
	}
	
	samplesWritten = 0;
		
	// cleanup_message(response); -- can't even call this in C++ 
	if (response) {
		if (response->def) free(response->def);
		if (response->buf) free(response->buf);
		free(response);
	}
}
	
bool PixelDataGrabber::writePixelData() {
	datadef_t data_def;
	message_t request;
	messagedef_t request_def;
	message_t *response = NULL;
	bool done = false;

	data_def.nchans = (UINT32_T) readResolution*phaseResolution*numSlices;
	data_def.nsamples = 1;
	data_def.data_type = DATATYPE_UINT16; 
	data_def.bufsize = sliceBuffer.size();
	
	// print_datadef(&data_def);
	
	request.def = &request_def;
	request_def.version = VERSION;
	request_def.command = PUT_DAT;
	request.buf = NULL;
	
	request_def.bufsize = sizeof(data_def) + data_def.bufsize;
	
	char *reqbuf = (char *) malloc(request_def.bufsize);
	if (reqbuf == NULL) {
		fprintf(stderr, "Out of memory\n");
		exit(1);
	}
	
	memcpy(reqbuf, &data_def, sizeof(data_def));
	memcpy(reqbuf + sizeof(data_def), sliceBuffer.data(), data_def.bufsize);
	
	request.buf = reqbuf;
	// print_request(request.def);

	/* write the request, read the response */
	DWORD t0 = timeGetTime();
	int result = tcprequest(ftbSocket, &request, &response);
	DWORD t1 = timeGetTime();	
	printf("Time lapsed while writing data: %li ms\n", t1-t0);
	
	if (result < 0) {
		fprintf(stderr, "Communication error when sending pixel data to fieldtrip buffer\n");
	}

	if (!response || !response->def) {
		fprintf(stderr, "PUT_DAT: unknown error in response\n");
	} else {
		if (response->def->command!=PUT_OK) {
			fprintf(stderr, "PUT_DAT: Buffer returned an error (%d)\n", response->def->command);
		} else {
			done = true;
		}
	}
	// cleanup_message(response); -- can't even call this in C++ 
	if (response) {
		if (response->def) free(response->def);
		if (response->buf) free(response->buf);
		free(response);
	}
	free(reqbuf);
	return done;
}	
	
void PixelDataGrabber::writeTimestampEvent(const struct timeval &tv) {
	DWORD t0, t1;
	struct {
		eventdef_t def;
		char buf[40];
	} event;
		
	message_t request;
	messagedef_t request_def;
	message_t *response = NULL;
			
	event.def.type_type = DATATYPE_CHAR;
	event.def.type_numel = 8;
	event.def.value_type = DATATYPE_CHAR;
	event.def.value_numel = 11+1+6;
	event.def.sample = samplesWritten - 1;
	event.def.offset = 0;
	event.def.duration = 0;
	event.def.bufsize = event.def.type_numel + event.def.value_numel;
	sprintf(event.buf, "unixtime%11li.%06li", tv.tv_sec, tv.tv_usec);
			
	request.def = &request_def;
	request_def.version = VERSION;
	request_def.command = PUT_EVT;
	request_def.bufsize = sizeof(eventdef_t) + event.def.bufsize;
	request.buf = &event;

	t0 = timeGetTime();
	int result = tcprequest(ftbSocket, &request, &response);
	t1 = timeGetTime();
	printf("Time lapsed while writing event: %li ms\n", t1-t0);
	
	if (result < 0) {
		fprintf(stderr, "Communication error when sending event to fieldtrip buffer\n");
	}
	
	if (!response || !response->def) {
		fprintf(stderr, "PUT_EVT: unknown error in response\n");
	} else {
		if (response->def->command!=PUT_OK) {
			fprintf(stderr, "PUT_EVT: Buffer returned an error (%d)\n", response->def->command);
		}
	}
	// cleanup_message(response); -- can't even call this in C++ 
	if (response) {
		if (response->def) free(response->def);
		if (response->buf) free(response->buf);
		free(response);
	}
}


void PixelDataGrabber::sendFrameToBuffer(const struct timeval &tv) {
	if (protInfo == NULL) tryReadProtocol(); // look for <WATCH_FOLDER>/mrprot.txt
	
	if (samplesWritten == 0) {
		writeHeader();
	}
	if (writePixelData()) {
		writeTimestampEvent(tv);
		samplesWritten++;
	}
}


void PixelDataGrabber::handleProtocol(const char *info, unsigned int sizeInBytes) {
	const sap_item_t *item;
	sap_item_t *PI = sap_parse(info, sizeInBytes);
	
	if (PI == NULL) return;
	
	/* new info? replace the old one, if present */
	if (protInfo != NULL) {
		sap_destroy(protInfo);
	}

	item = sap_search_deep(PI, "sKSpace.lBaseResolution");
	if (item!=NULL && item->type == SAP_LONG) {
		long res = *((long *) item->value);
		readResolution = (res > 0) ? res : 0;
	}
	
	item = sap_search_deep(PI, "sSliceArray.lSize");
	if (item!=NULL && item->type == SAP_LONG) {
		long slices = *((long *) item->value);
		numSlices = (slices > 0) ? slices : 0;
	}
	
	item = sap_search_deep(PI, "sSliceArray.asSlice[0].dPhaseFOV");
	if (item!=NULL && item->type == SAP_DOUBLE) {
		phaseFOV = *((double *) item->value);
		printf("PhaseFOV: %f\n", phaseFOV);
	}
	
	item = sap_search_deep(PI, "sSliceArray.asSlice[0].dReadoutFOV");
	if (item!=NULL && item->type == SAP_DOUBLE) {
		readoutFOV = *((double *) item->value);
		printf("ReadoutFOV: %f\n", readoutFOV);
	}
	
	if (phaseFOV > 0.0 && readoutFOV > 0.0) {
		phaseResolution = (unsigned int) round(readResolution * phaseFOV / readoutFOV);
	} else {
		phaseResolution = 0;
	}
	printf("Resolution: %i x %i x %i\n", readResolution, phaseResolution, numSlices);
		
	protInfo = PI;
}

bool PixelDataGrabber::tryReadFile(const char *filename, SimpleBuffer &sBuf) {
	HANDLE fHandle;

	fHandle = CreateFile(filename, GENERIC_READ, 0, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
	
	if (fHandle != INVALID_HANDLE_VALUE) {
		DWORD read, size;
			
		size = GetFileSize(fHandle, NULL);  // ignore higher DWORD - our files are not that big
		
		if (!sBuf.resize(size)) {
			fprintf(stderr, "Out of memory!\n");
			CloseHandle(fHandle);
			return false;
		}
			
		ReadFile(fHandle, sBuf.data(), size, &read, NULL);
		CloseHandle(fHandle);
		if (size!=read) {
			fprintf(stderr, "Error reading file contents from disk.\n");
			sBuf.resize(read);
			return false;
		} 
		return true;
	}
	return false;
}


void PixelDataGrabber::tryReadProtocol() {
	std::string defProtFile = sourceDir + "mrprot.txt";
	
	if (!tryReadFile(defProtFile.c_str(), protBuffer)) return;
	
	handleProtocol((char *) protBuffer.data(), protBuffer.size());
}


void PixelDataGrabber::tryFolderToBuffer() {
	const std::vector<std::string>& vfn = FW->getFilenames();
	
	for (unsigned int i=0;i<vfn.size();i++) {
		// both .PixelData and "mrprot.txt" are at least 10 characters long
		if (vfn[i].size() < 10) continue;
		
		fullName = sourceDir + vfn[i];
		
		if (vfn[i].compare(vfn[i].size()-10,10, ".PixelData") == 0) {
			struct timeval tv;
			
			gettimeofday(&tv, NULL);
			if (!tryReadFile(fullName.c_str(), pixBuffer)) continue;
			reshapeToSlices();
			if (sliceBuffer.size()==0) continue;
			sendFrameToBuffer(tv);
		} 
		else if (vfn[i].compare(vfn[i].size()-10,10, "mrprot.txt") == 0) {
			if (!tryReadFile(fullName.c_str(), protBuffer)) return;
	
			handleProtocol((char *) protBuffer.data(), protBuffer.size());
			
			lastAction = ProtocolRead;
		} 
		else {
			// nothing of interest
		}
	}
}

void PixelDataGrabber::reshapeToSlices() {
	unsigned int pixels = pixBuffer.size() >> 1;
	unsigned int root   = (unsigned int) round(sqrt(pixels));
	
	if (readResolution == 0 || phaseResolution == 0 || numSlices == 0) {
		// if for some reason we don't have proper protocol information,
		// we just spit out square images (one slice)
		if (root*root != pixels) {
			// source image not square? don't know what to do
			sliceBuffer.resize(0);
			lastAction = BadPixelData;
		} else	{
			readResolution = phaseResolution = root;
			numSlices = 1;
			if (sliceBuffer.resize(pixels*2)) {
				memcpy(sliceBuffer.data(), pixBuffer.data(), pixels*2);
				lastAction = PixelsTransmitted;
			} else {
				fprintf(stderr, "Out of memory!!!\n");
				sliceBuffer.resize(0);
				lastAction = OutOfMemory;
			}
		}
	} else {
		unsigned int tiles  = pixels / (readResolution * phaseResolution);
		unsigned int mosw = (unsigned int) round(sqrt(tiles));
		
		if ((pixels != readResolution * phaseResolution * tiles) || (mosw*mosw != tiles) || (numSlices > tiles)) {
			// mosaic does not match readResolution, or is too small - do nothing...
			fprintf(stderr, "PixelData (%i) does not match protocol information (%i x %i x %i)\n",pixels,readResolution,phaseResolution,numSlices);
			sliceBuffer.resize(0);
			lastAction = BadPixelData;
		} else if (!sliceBuffer.resize(readResolution*phaseResolution*numSlices*2)) {
			fprintf(stderr, "Out of memory!!!\n");
			sliceBuffer.resize(0);
			lastAction = OutOfMemory;
		} else {
			const uint16_t *src	= (const uint16_t *) pixBuffer.data();
			uint16_t *dest   	= (uint16_t *) sliceBuffer.data();
			
			// row and column in source mosaic
			unsigned int i = 0, j = 0;
			
			for (unsigned int n=0; n<numSlices; n++) {
				uint16_t *dest_n = dest + readResolution*phaseResolution*n;
				const uint16_t *src_n = src + readResolution*(phaseResolution*mosw*i + j);
				
				// copy one line (=readResolution pixels) at a time
				for (unsigned int m=0; m<phaseResolution; m++) {
					memcpy(dest_n + readResolution*m, src_n + mosw*readResolution*m, 2 * readResolution);
				}
				
				// increase source column index
				if (++j == mosw) {
					// and shift to next row eventually
					++i;
					j = 0;
				}
			}
			lastAction = PixelsTransmitted;
		}
	}
}
