/*
 * Copyright (C) 2010, Stefan Klanke
 * Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
 * Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 */

#include <vector>
#include <string>
#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <stdlib.h>

#ifdef WIN32
#include <windows.h>
#else
#include <unistd.h>
#include <sys/time.h>
#endif

#include <nifti1.h>
#include <FtBuffer.h>

nifti_1_header commonHeader;
std::vector<std::string> niiFiles;
void *dataBuffer;
unsigned int dataSize;
int nChans;
int ftSocket;
FtBufferRequest ftReq;
int interval = 2000; // in milliseconds

bool isRelative(const char *fn) {
	while (*fn != 0) {
		if (*fn=='/' || *fn=='\\') return false;
		fn++;
	}
	return true;
}

bool areEqual(const nifti_1_header &A, const nifti_1_header &B) {
	if (A.datatype != B.datatype) return false;
	if (A.bitpix != B.bitpix) return false;
	for (int i=0;i<8;i++) {
		if (A.dim[i] != B.dim[i]) return false;
		if (A.pixdim[i] != B.pixdim[i]) return false;
	}
	// everything else, we don't care
	return true;
}

int readAndCheckScannerFile(char *filename) {
	nifti_1_header thisHeader;
	std::string path;
	FILE *f;
	int len = 0;
	char *last;
	
	last = strrchr(filename,'\\');
	if (last != NULL) {
		len = last-filename;
	}
	last = strrchr(filename,'/');
	if (last != NULL) {
		int lf = last-filename;
		if (lf>len) len = lf;
	}
	if (len>0) {
		path.assign(filename, len+1);
	}
	
	f = fopen(filename, "r");
	if (f==NULL) return -1;
	
	while (!feof(f)) {
		char line[512];
		
		fgets(line, 512, f);
		len = strlen(line);
		while (len>0 && isspace(line[len-1])) len--;
		
		if (len==0) continue;
		
		std::string fullname = path;
		fullname.append(line,len);
		
		FILE *nf = fopen(fullname.c_str(), "rb");
		if (nf==NULL) {
			fprintf(stderr,"Cannot open %s\n", fullname.c_str());
			continue;
		}
		
		int read = fread(&thisHeader, 1, sizeof(thisHeader), nf);
		if (read == sizeof(thisHeader)) {
			if (niiFiles.size() == 0 && strcmp(thisHeader.magic, "n+1")==0) {
				commonHeader = thisHeader;
				niiFiles.push_back(fullname);
			} else {
				if (!areEqual(thisHeader, commonHeader)) {
					fprintf(stderr, "Header of %s does not equal header of %s\n",fullname.c_str(), niiFiles[0].c_str());
				} else {
					niiFiles.push_back(fullname);
				}
			}
		} else {
			fprintf(stderr,"Could not read NIFTI-1 header from %s\n", fullname.c_str());
		}
		fclose(nf);
	}
	
	fclose(f);
	return niiFiles.size();
}

bool grabScan(int i) {
	FILE *nf = fopen(niiFiles[i].c_str(),"rb");
	if (nf == NULL) return false;
	
	fseek(nf, 352, SEEK_SET);
	if (fread(dataBuffer, 1, dataSize, nf) != dataSize) {
		printf("Warning: Could not read all voxels from %s\n", niiFiles[i].c_str());
	}
	
	fclose(nf);
	return true;
}

bool writeHeader() {
	FtBufferResponse resp;
	
	if (ftReq.prepPutHeader(nChans, DATATYPE_INT16, 1000.0 / (float) interval) && 
		ftReq.prepPutHeaderAddChunk(FT_CHUNK_NIFTI1, sizeof(commonHeader), &commonHeader)) {

		tcprequest(ftSocket, ftReq.out(), resp.in());
	} else {
		fprintf(stderr, "Out of memory!\n");
		exit(1);
	}
	return resp.checkPut();
}

bool writeScan(int i, const struct timeval &tv) {
	FtBufferResponse resp;
	//char ts[20];
	
	if (!ftReq.prepPutData(nChans, 1, DATATYPE_INT16, dataBuffer)) {
		fprintf(stderr, "Out of memory!\n");
		exit(1);
	}
	
	tcprequest(ftSocket, ftReq.out(), resp.in());
	
	if (!resp.checkPut()) return false;
	
	//sprintf(ts, "%11li.%06li", tv.tv_sec, tv.tv_usec);
	
	if (!ftReq.prepPutEvent(i, 0, 0, "scan","ready")) {
		fprintf(stderr, "Out of memory!\n");
		exit(1);
	}
	tcprequest(ftSocket, ftReq.out(), resp.in());
	return resp.checkPut();
}

int main(int argc, char *argv[]) {
	char hostname[256];
	int port;
	struct timeval tvBefore, tvAfter;
	
	if (sizeof(nifti_1_header) != 348) {
		fprintf(stderr, "Struct nifti_1_header has bad size (%i bytes, should be 348)\n", sizeof(nifti_1_header));
	}
	
	#ifdef WIN32
	timeBeginPeriod(1);
	#endif
	
	if (argc<2) {
		printf("Usage:\n  nii_to_buffer  path_to/niscanner.txt  [deltaMilliSec=2000 [hostname=localhost [port=1972]]]");
		return 1;
	}
	
	printf("Reading %s and checking listed files...\n", argv[1]);
	int numFiles = readAndCheckScannerFile(argv[1]);
	
	if (numFiles == -1) {
		fprintf(stderr, "Could not open scanner description file %s\n",argv[1]);
		return 1;
	}
	if (numFiles == 0) {
		fprintf(stderr, "Could not read a single NII file :-(\n");
	}
	printf("Will transmit %i scans\n", numFiles);
	
	if (argc>=3) {
		interval = atoi(argv[2]);
	}
	
	if (argc>=4) {
		strncpy(hostname, argv[3], 256);
	} else {
		strncpy(hostname, "localhost", 256);
	}
	
	if (argc>=5) {
		port = atoi(argv[4]);
	} else {
		port = 1972;
	}
	
	if (commonHeader.datatype != DT_INT16) { 
		fprintf(stderr, "Sorry, datatype != DT_INT16, don't know what to do.\n");
		return 1;
	}
	nChans = commonHeader.dim[1]*commonHeader.dim[2]*commonHeader.dim[3];
	if (nChans == 0) {
		fprintf(stderr, "Sorry, one of the dimensions is 0 - don't know what to do.\n");
		return 1;
	}
	dataSize = sizeof(UINT16_T) * nChans;
	dataBuffer = calloc(dataSize,1);
	if (dataBuffer == NULL) {
		fprintf(stderr, "Sorry, could not allocate memory for caching voxel data.\n");
		return 1;
	}
	
	printf("Trying to connect to fieldtrip buffer on %s:%d\n", hostname, port);
	ftSocket = open_connection(hostname, port);
	if (ftSocket == -1) {
		fprintf(stderr, "Sorry, connection to FieldTrip failed.\n");
		free(dataBuffer);
		return 1;
	}
	printf("Ok!\n\n");
				
	gettimeofday(&tvBefore, NULL);
	if (writeHeader()) {
		printf("Wrote header at unixtime = %li.%06li\n", tvBefore.tv_sec, tvBefore.tv_usec);
		gettimeofday(&tvAfter, NULL);
		
		for (unsigned int i=0;i<niiFiles.size();i++) {
			#ifdef WIN32
			int milliSec = interval - 1000*(tvAfter.tv_sec - tvBefore.tv_sec) - (tvAfter.tv_usec - tvBefore.tv_usec)/1000;
			if (milliSec>0) Sleep(milliSec);
			#else
			int usec = 1000*(interval - 1000*(tvAfter.tv_sec - tvBefore.tv_sec)) - (tvAfter.tv_usec - tvBefore.tv_usec);
			if (usec>0) usleep(usec);
			#endif
		
			gettimeofday(&tvBefore, NULL);
			grabScan(i);
			if (!writeScan(i,tvBefore)) continue;
			printf("Wrote scan %i at unixtime = %li.%06li\n", i, tvBefore.tv_sec, tvBefore.tv_usec);
			
			gettimeofday(&tvAfter, NULL);
		}
	}
	
	close_connection(ftSocket);
	free(dataBuffer);
	return 0;
}
