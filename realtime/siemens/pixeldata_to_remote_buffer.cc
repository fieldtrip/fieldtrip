#include <stdio.h>
#include <FolderWatcher.h>
#include <math.h>

#include <message.h>
#include <buffer.h>


class PixelDataReader {
	protected:
		
	HANDLE fwEventHandle;
	int ftbSocket;
	int thisNumPixels;
	int lastWrittenSample;
	std::string sourceDir;
	FolderWatcher *FW;
	bool isValid;
	
	public: 
	
	PixelDataReader() {
		ftbSocket = -1;
		thisNumPixels = 0;
		FW = NULL;
		isValid = false;
		
		fwEventHandle = CreateEvent(NULL, FALSE, FALSE, NULL);
		// if (fwEventHandle == INVALID_HANDLE_VALUE) return;
	}
	
	~PixelDataReader() {
		printf("Cleaning up...\n");
		if (fwEventHandle != INVALID_HANDLE_VALUE) CloseHandle(fwEventHandle);
		if (FW) delete FW;
		
		if (ftbSocket > 0) close_connection(ftbSocket);
	}
	
	int monitor(const char *directory) {
		if (directory == 0) return 0;
		
		//mutex.lock();
		
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
		// mutex.unlock();
		return FW->isValid();
	}
	
	int connect(const char *hostname, int port) {
		if (ftbSocket > 0) close_connection(ftbSocket);
		ftbSocket = open_connection(hostname, port);
		return (ftbSocket > 0);
	}
	
	void run() {
		if (ftbSocket <= 0) {
			fprintf(stderr, "No fieldtrip buffer has been connected to.\n");
			return;
		}
		if (FW==NULL || !FW->isValid()) {
			fprintf(stderr, "No directory is being monitored.\n");
			return;
		}

		DWORD waitResult = WaitForSingleObject(fwEventHandle, 200);
			
		if (waitResult == WAIT_TIMEOUT) return;
			
		if (waitResult == WAIT_FAILED) {
			fprintf(stderr, "Error in WaitForSingleObject(...) !\n");
			// TODO: proper error handling
			return;
		}

		if (FW->processChanges() > 0) {
			tryFolderToBuffer();
			FW->startListenForChanges(fwEventHandle);
		}
	}

	protected:	
	
	void writeHeader(int numPixels) {
		headerdef_t header_def;

		message_t request;
		messagedef_t request_def;
		message_t *response = NULL;
				
		header_def.nchans = (UINT32_T) numPixels;
		header_def.nsamples = 1;
		header_def.nevents = 1;
		header_def.fsample = 2.0; /* TODO: is this always correct ? */
		header_def.data_type = DATATYPE_INT16;
		header_def.bufsize = 0;
		
		request.def = &request_def;
		request.buf = &header_def;
		request_def.version = VERSION;
		request_def.command = PUT_HDR;
		request_def.bufsize = sizeof(headerdef_t);
		
		// print_headerdef(&header_def);
		
		// print_request(request.def);

		/* write the request, read the response */
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
		
		thisNumPixels = numPixels;
		lastWrittenSample = -1;
		
		// cleanup_message(response); -- can't even call this in C++ 
		if (response) {
			if (response->def) free(response->def);
			if (response->buf) free(response->buf);
			free(response);
		}
	}
	
	void writePixelData(void *pixelData) {
		datadef_t data_def;
		message_t request;
		messagedef_t request_def;
		message_t *response = NULL;

		data_def.nchans = (UINT32_T) thisNumPixels;
		data_def.nsamples = 1;
		data_def.data_type = DATATYPE_INT16;  // it's UINT16 actually, but the buffer can't do that yet
		data_def.bufsize = sizeof(UINT16_T)*thisNumPixels;
		
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
		memcpy(reqbuf + sizeof(data_def), pixelData, data_def.bufsize);
		
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
				lastWrittenSample++;
			}
		}
		// cleanup_message(response); -- can't even call this in C++ 
		if (response) {
			if (response->def) free(response->def);
			if (response->buf) free(response->buf);
			free(response);
		}
		free(reqbuf);
	}	
	
	void writeTimestampEvent(const struct timeval &tv) {
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
		event.def.sample = lastWrittenSample;
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
	
	void sendToBuffer(void *pixelData, int sizeInBytes, const struct timeval &tv) {
		int numPixels = sizeInBytes / 2;
		
		if (thisNumPixels != numPixels) {
			writeHeader(numPixels);
		}
		writePixelData(pixelData);
		writeTimestampEvent(tv);
	}
	
	void tryFolderToBuffer() {
		const std::vector<std::string>& vfn = FW->getFilenames();
		
		for (unsigned int i=0;i<vfn.size();i++) {
			if ((vfn[i].size() > 9) && 	(vfn[i].compare(vfn[i].size()-9,9, "PixelData") == 0)) {
				std::string fullName = sourceDir + vfn[i];
			
				HANDLE fHandle = CreateFile(fullName.c_str(), GENERIC_READ, 0, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
				if (fHandle != INVALID_HANDLE_VALUE) {
					DWORD read, size;
					void *pixelData;
					struct timeval tv;
					
					gettimeofday(&tv, NULL);
					
					size = GetFileSize(fHandle, NULL);  // ignore higher DWORD - images are not that big
					
					pixelData = malloc(size);
					if (pixelData == NULL) {
						fprintf(stderr, "Out of memory!\n");
						exit(1);
					}
				
					ReadFile(fHandle, pixelData, size, &read, NULL);
					if (size == read) {
						sendToBuffer(pixelData, size, tv);
					} else {
						fprintf(stderr, "Error reading pixel data from disk.\n");
					}
					// printf("\nR %16.4f   %s  %i / %i Byte\n", timestamp, fullName.c_str(), read, size);
					CloseHandle(fHandle);
					free(pixelData);
				
				} else {
					// Can't open without sharing: this means the scanner is still writing
					// Do nothing until the next notification
				}
			}
		}
	}
};



int main(int argc, char *argv[]) {
	PixelDataReader pdr;
	char hostname[256];
	int port;
	char directory[256];
	
	timeBeginPeriod(1);
	
	if (argc>=2) {
		strncpy(hostname, argv[1], 256);
	} else {
		strncpy(hostname, "localhost", 256);
	}
	
	if (argc>=3) {
		port = atoi(argv[2]);
	} else {
		port = 1972;
	}
	
	if (argc>=4) {
		strncpy(directory, argv[3], 256);
	} else {
		strncpy(directory, "E:\\IMAGE", 256);
	}
	
	printf("Trying to connect to fieldtrip buffer on %s:%d\n", hostname, port);
	
	if (pdr.connect(hostname, port)) {
		printf("OK\n");
	} else {
		printf("Failure!\n");
		exit(1);
	}
	
	printf("Trying to open directory %s for monitoring...\n", directory);
	if (pdr.monitor(directory)) {
		printf("OK\n");
	} else {
		printf("Failure!\n");
		exit(1);
	}
	
	printf("Starting to monitor... - press CTRL-C to cancel\n");
	
	while (1) {
		pdr.run();
	}
	
	return 0;
}
