/** Simple C++ class for managing ASCII request from a TCP port.
	(C) 2010 S. Klanke
*/
#ifndef __ControlServer_h
#define __ControlServer_h

#include <stdio.h>
#include <stdint.h>
#include <string.h>

#ifdef WIN32
#include <windows.h>
#else
#define SOCKET          int
#define INVALID_SOCKET  -1
#endif

class ControlServer {
	public:
	
	ControlServer(int bufSize = 2048) {
		bufRead  = (char *) malloc(bufSize);
		sizeRead = (bufRead == NULL) ? 0 : bufSize;

		bufWrite = (char *) malloc(bufSize);
		sizeWrite= (bufWrite == NULL) ? 0 : bufSize;

		numRead = numWrite = state = 0;
	}
	
	~ControlServer() {
		close();
	}
	
	boolean listen(int port) {
		return true;
	}
	
	void close() {
	}
	
	const char *checkRequest() {int milliSeconds = 0) {
		struct timeval tv = {0, 1000*milliSeconds};
		return NULL;
	}
	
	boolean respond(const char *msg) {
		int n = strlen(msg) + 1;
		
		if (state != 2) return false;
		
		if (sizeWrite < n) {
			if (bufWrite != NULL) free(bufWrite);
			bufWrite = (char *) malloc(n);
			if (bufWrite == NULL) {
				fprintf(stderr, "ControlServer.writeResponse: out of memory!\n");
				bufSize = 0;
				return false;
			}
			bufSize = n;
		}
		memcpy(bufWrite, msg, n);
		numWrite = 0;
		state = 3;
	}
	
	protected:
	
	void checkWrite() {
		struct timeval tv = {0, 0};
		FD_ZERO(&fdSett);
		if (state < 2) {
			FD_SET(sock, &readSet);
		} else {
			FD_SET(sock, &writeSet);
		}
		sel = select((int) sock+1, &readSet, &writeSet, NULL, &tv);
	}
		
	SOCKET servSocket, clientSocket;
	
	char *bufRead, *bufWrite;   // buffer to store input and response
	int sizeRead, sizeWrite;    // size of allocated buffers
	int numRead, numWrite;		// number of bytes processed so far
	int state;        			// 0 = nothing, 1 = reading, 2 = read ok, 3 = writing
	FD_SET readSet, writeSet;
};

#endif
