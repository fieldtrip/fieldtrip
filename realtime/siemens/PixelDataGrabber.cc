/*
 * Copyright (C) 2010, Stefan Klanke
 * Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
 * Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 */
#include <PixelDataGrabber.h>
#include <stdio.h>
#include <math.h>
#include <FolderWatcher.h>

PixelDataGrabber::PixelDataGrabber() {
	ftbSocket = -1;
	FW = NULL;
	fwEventHandle = INVALID_HANDLE_VALUE;
	
	phaseFOV = readoutFOV = 0.0;
	samplesWritten = 0;
	readResolution = 0;
	numSlices = 0;
	TR = 0;
	lastAction = Nothing;
	protInfo = NULL;
	fwActive = false;
	headerWritten = false;
	verbosity = 1;
}

PixelDataGrabber::~PixelDataGrabber() {
	if (fwEventHandle != INVALID_HANDLE_VALUE) CloseHandle(fwEventHandle);
	if (FW) delete FW;
		
	if (ftbSocket > 0) close_connection(ftbSocket);
	sap_destroy(protInfo);
}
	
bool PixelDataGrabber::monitorDirectory(const char *directory) {
	if (fwEventHandle == INVALID_HANDLE_VALUE) {
		fwEventHandle = CreateEvent(NULL, FALSE, FALSE, NULL);
		if (fwEventHandle == INVALID_HANDLE_VALUE) return false;
	}
	
	if (protInfo) {
		sap_destroy(protInfo);
		protInfo = NULL;
	}
	
	if (FW) {
		FW->stopListenForChanges();
		delete FW;
		FW = NULL;
		fwActive = false;
	}
	if (directory == NULL) return false;

	sourceDir = directory;
	int n = sourceDir.size();
	if (n>0) {
		if (sourceDir[n-1]!='\\' && sourceDir[n-1]!='/') sourceDir += '\\';

		FW = new FolderWatcher(directory);
		fwActive = FW->startListenForChanges(fwEventHandle);
	}

	return fwActive;
}

bool PixelDataGrabber::connectToFieldTrip(const char *hostname, int port) {
	if (ftbSocket > 0) close_connection(ftbSocket);
	if (hostname != NULL) {
		ftbSocket = open_connection(hostname, port);
	} else {
		ftbSocket = -1;
	}
	headerWritten = false;
	return (ftbSocket > 0);
}
	
int PixelDataGrabber::run(unsigned int ms) {
	lastAction = Nothing;
	if (FW==NULL || !FW->isValid()) {
		if (verbosity>3) fprintf(stderr, "No directory is being monitored.\n");
		return -1;
	}

	DWORD waitResult = WaitForSingleObject(fwEventHandle, ms);

	if (waitResult == WAIT_FAILED) {
		if (verbosity>0) fprintf(stderr, "Error in WaitForSingleObject(...) !\n");
		// TODO: proper error handling
		return -1;
	}
		
	if (waitResult == WAIT_TIMEOUT) return 0;

	if (FW->processChanges() > 0) {
		tryFolderToBuffer();
		FW->startListenForChanges(fwEventHandle);
	}
	return 1;
}

	
bool PixelDataGrabber::writeHeader() {
	FtBufferResponse resp;
	UINT32_T nchans = numSlices*readResolution*phaseResolution;
	float    fsamp  = 1.0e6 / (float) TR;
	
	headerWritten = false;
	
	if (ftReq.prepPutHeader(nchans, DATATYPE_INT16, fsamp) &&
		ftReq.prepPutHeaderAddChunk(FT_CHUNK_SIEMENS_AP, protBuffer.size(), protBuffer.data()) &&
		ftReq.prepPutHeaderAddChunk(FT_CHUNK_NIFTI1, sizeof(nifti), &nifti)) {
		
		// write the request, read the response
		int result = tcprequest(ftbSocket, ftReq.out(), resp.in());		
	
		if (result < 0) {
			if (verbosity>0) fprintf(stderr, "Communication error when sending header to fieldtrip buffer\n");
			lastAction = TransmissionError;
			return false;
		}
	
		if (!resp.checkPut()) {
			if (verbosity>0) fprintf(stderr, "PUT_HDR: error from buffer server\n");
			lastAction = TransmissionError;
			return false;
		}
		headerWritten = true;
		samplesWritten = 0;
		return true;
	} else {
		// one of the prepPutHeader* calls failed
		lastAction = OutOfMemory;
		return false;
	}
}
	
bool PixelDataGrabber::writePixelData() {
	FtBufferResponse resp;
	UINT32_T nchans = numSlices*readResolution*phaseResolution;
	
	if (!ftReq.prepPutData(nchans, 1, DATATYPE_INT16, sliceBuffer.data())) {
		if (verbosity > 0) fprintf(stderr, "Out of memory!\n");
		lastAction = OutOfMemory;
		return false;
	}
	
	/* write the request, read the response */
	DWORD t0 = timeGetTime();
	int result = tcprequest(ftbSocket, ftReq.out(), resp.in());
	DWORD t1 = timeGetTime();	
	if (verbosity>2) printf("Time lapsed while writing data: %li ms\n", t1-t0);
	
	if (result < 0) {
		if (verbosity>0) fprintf(stderr, "Communication error when sending pixel data to fieldtrip buffer\n");
		lastAction = TransmissionError;
		return false;
	}
	if (!resp.checkPut()) {
		if (verbosity>0) fprintf(stderr, "PUT_DAT: error from buffer server\n");
		lastAction = TransmissionError;
		return false;
	}
	samplesWritten++;
	lastAction = PixelsTransmitted;
	return true;
}	
	
bool PixelDataGrabber::writeTimestampEvent(const struct timeval &tv) {
	FtBufferResponse resp;
	char ts[20];
	DWORD t0, t1;
	
	sprintf(ts, "%11li.%06li", tv.tv_sec, tv.tv_usec);
	
	if (!ftReq.prepPutEvent(samplesWritten-1, 0, 0, "unixtime", ts)) {
		if (verbosity>0) fprintf(stderr, "Out of memory!\n");
		lastAction = OutOfMemory;
		return false;
	}

	t0 = timeGetTime();
	int result = tcprequest(ftbSocket, ftReq.out(), resp.in());
	t1 = timeGetTime();
	if (verbosity>2) printf("Time lapsed while writing event: %li ms\n", t1-t0);
	
	if (result < 0) {
		if (verbosity>0) fprintf(stderr, "Communication error when sending pixel data to fieldtrip buffer\n");
		lastAction = TransmissionError;
		return false;
	}
	if (!resp.checkPut()) {
		if (verbosity>0) fprintf(stderr, "PUT_EVT: error from buffer server\n");
		lastAction = TransmissionError;
		return false;
	}
	return true;
}

bool PixelDataGrabber::sendFrameToBuffer(const struct timeval &tv) {
	if (!headerWritten) {
		if (!writeHeader()) return false;
	}
	if (!writePixelData()) return false;
	if (!writeTimestampEvent(tv)) return false;
	return true;
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
		if (verbosity>2) printf("PhaseFOV: %f\n", phaseFOV);
	}
	
	item = sap_search_deep(PI, "sSliceArray.asSlice[0].dReadoutFOV");
	if (item!=NULL && item->type == SAP_DOUBLE) {
		readoutFOV = *((double *) item->value);
		if (verbosity>2) printf("ReadoutFOV: %f\n", readoutFOV);
	}
	
	if (phaseFOV > 0.0 && readoutFOV > 0.0) {
		phaseResolution = (unsigned int) round(readResolution * phaseFOV / readoutFOV);
	} else {
		phaseResolution = 0;
	}
	if (verbosity>2) printf("Resolution: %i x %i x %i\n", readResolution, phaseResolution, numSlices);
	
	item = sap_search_deep(PI, "alTR");
	if (item!=NULL && item->type == SAP_LONG) {
		TR = ((long *) item->value)[0];
		if (verbosity>2) printf("TR: %i microsec.\n", TR);
	} else {
		TR = 2000000;
	}
	
	memset(&nifti, 0, sizeof(nifti));
	nifti.sizeof_hdr = 348;
	nifti.dim[0] = 3;
	nifti.dim[1] = (short) readResolution;
	nifti.dim[2] = (short) phaseResolution;
	nifti.dim[3] = (short) numSlices;
	nifti.datatype = DT_INT16;
	nifti.bitpix = 16;
	nifti.slice_start = 0;
	nifti.scl_slope = 1.0f;
	nifti.slice_end = (short) numSlices-1;
	nifti.xyzt_units = NIFTI_UNITS_MM | NIFTI_UNITS_SEC;
	nifti.slice_duration = (TR/numSlices)*1e-6f;	// in seconds
	
	nifti.pixdim[0] = 1.0f;		// TODO: this is qfac - check if this needs to be -1 
	nifti.pixdim[1] = (float) (readoutFOV / readResolution);
	nifti.pixdim[2] = (float) (phaseFOV / phaseResolution);
	nifti.pixdim[4] = TR*1e-6f;	// repetition time in seconds 
	
	item = sap_search_deep(PI, "sSliceArray.asSlice[0].dThickness");
	if (item!=NULL && item->type == SAP_DOUBLE) {
		nifti.pixdim[3] = (float) *((double *) item->value);
	}

	long ucMode = 1;
	item = sap_search_deep(PI, "sSliceArray.ucMode");
	if (item!=NULL && item->type == SAP_LONG) {
		ucMode = *((long *) item->value);
	}
	switch(ucMode) {
		case 1:
			nifti.slice_code = NIFTI_SLICE_SEQ_INC;
			break;
		case 2:
			nifti.slice_code = NIFTI_SLICE_SEQ_DEC;
			break;
		case 4:
			// in case of an odd number of slices, start at zero (ALT), 
			// otherwise slice[1] (=second!) has been acquired first (ALT2)
			nifti.slice_code = (numSlices & 1) ? NIFTI_SLICE_ALT_INC : NIFTI_SLICE_ALT_INC2;
			break;
	}
	
	float px = 0.0;
	float py = 0.0;
	float pz = 0.0;
	item = sap_search_deep(PI, "sSliceArray.asSlice[0].sPosition.dSag");
	if (item!=NULL && item->type == SAP_DOUBLE) {
		px = (float) *((double *) item->value);
	}	
	item = sap_search_deep(PI, "sSliceArray.asSlice[0].sPosition.dCor");
	if (item!=NULL && item->type == SAP_DOUBLE) {
		py = (float) *((double *) item->value);
	}	
	item = sap_search_deep(PI, "sSliceArray.asSlice[0].sPosition.dTra");
	if (item!=NULL && item->type == SAP_DOUBLE) {
		pz = (float) *((double *) item->value);
	}	
	// TODO: squeeze out correct rotation from SAP info
	nifti.qform_code = NIFTI_XFORM_SCANNER_ANAT;
	nifti.qoffset_x = px;
	nifti.qoffset_y = py;
	nifti.qoffset_z = pz;
	
	nifti.sform_code = NIFTI_XFORM_SCANNER_ANAT;
	nifti.srow_x[0] = 1.0;
	nifti.srow_x[3] = px;
	nifti.srow_y[1] = 1.0;
	nifti.srow_y[3] = py;
	nifti.srow_z[2] = 1.0;
	nifti.srow_z[3] = pz;
	
	nifti.magic[0] = 'n';
	nifti.magic[1] = 'i';
	nifti.magic[2] = '1';
		
	protInfo = PI;
}

bool PixelDataGrabber::tryReadFile(const char *filename, SimpleStorage &sBuf) {
	HANDLE fHandle;

	fHandle = CreateFile(filename, GENERIC_READ, 0, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
	
	if (fHandle != INVALID_HANDLE_VALUE) {
		DWORD read, size;
			
		size = GetFileSize(fHandle, NULL);  // ignore higher DWORD - our files are not that big
		
		if (!sBuf.resize(size)) {
			if (verbosity>0) fprintf(stderr, "Out of memory in tryReadFile !!!\n");
			CloseHandle(fHandle);
			return false;
		}
			
		ReadFile(fHandle, sBuf.data(), size, &read, NULL);
		CloseHandle(fHandle);
		if (size==read) return true;

		if (verbosity>0) fprintf(stderr, "Error reading file contents from disk.\n");
		sBuf.resize(read);
	}
	return false;
}


bool PixelDataGrabber::tryReadProtocol() {
	std::string defProtFile = sourceDir + "mrprot.txt";
	
	if (!tryReadFile(defProtFile.c_str(), protBuffer)) return false;
	
	handleProtocol((char *) protBuffer.data(), protBuffer.size());
	return true;
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
			
			if (protInfo == NULL) tryReadProtocol();
			reshapeToSlices();
			if (sliceBuffer.size()==0) continue;
			if (ftbSocket == -1) continue;
			sendFrameToBuffer(tv);
		} 
		else if (vfn[i].compare(vfn[i].size()-10,10, "mrprot.txt") == 0) {
			if (!tryReadFile(fullName.c_str(), protBuffer)) return;
	
			handleProtocol((char *) protBuffer.data(), protBuffer.size());
			if (ftbSocket != -1) writeHeader();
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
			if (sliceBuffer.resize(pixels*sizeof(INT16_T))) {
				memcpy(sliceBuffer.data(), pixBuffer.data(), pixels*sizeof(INT16_T));
			} else {
				if (verbosity>0) fprintf(stderr, "Out of memory in reshapeToSlices !!!\n");
				sliceBuffer.resize(0);
				lastAction = OutOfMemory;
			}
		}
	} else {
		unsigned int tiles  = pixels / (readResolution * phaseResolution);
		unsigned int mosw = (unsigned int) round(sqrt(tiles));
		
		if ((pixels != readResolution * phaseResolution * tiles) || (mosw*mosw != tiles) || (numSlices > tiles)) {
			// mosaic does not match readResolution, or is too small - do nothing...
			if (verbosity>0) fprintf(stderr, "PixelData (%i) does not match protocol information (%i x %i x %i)\n",pixels,readResolution,phaseResolution,numSlices);
			sliceBuffer.resize(0);
			lastAction = BadPixelData;
		} else if (!sliceBuffer.resize(readResolution*phaseResolution*numSlices*sizeof(INT16_T))) {
			if (verbosity>0) fprintf(stderr, "Out of memory in reshapeToSlices !!!\n");
			sliceBuffer.resize(0);
			lastAction = OutOfMemory;
		} else {
			const int16_t *src	= (const int16_t *) pixBuffer.data();
			int16_t *dest   	= (int16_t *) sliceBuffer.data();
			
			// row and column in source mosaic
			unsigned int i = 0, j = 0;
			
			for (unsigned int n=0; n<numSlices; n++) {
				int16_t *dest_n = dest + readResolution*phaseResolution*n;
				const int16_t *src_n = src + readResolution*(phaseResolution*mosw*i + j);
				
				// copy one line (=readResolution pixels) at a time
				for (unsigned int m=0; m<phaseResolution; m++) {
					memcpy(dest_n + readResolution*m, src_n + mosw*readResolution*m, sizeof(INT16_T) * readResolution);
				}
				
				// increase source column index
				if (++j == mosw) {
					// and shift to next row eventually
					++i;
					j = 0;
				}
			}
		}
	}
}
