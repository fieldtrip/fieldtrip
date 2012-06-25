#ifndef __LocalPipe_h
#define __LocalPipe_h

#ifdef WIN32
#include <windows.h>
#else
#include <unistd.h>
#include <sys/socket.h>
#endif


class LocalPipe {
	public:
	
	LocalPipe() {
		#ifdef WIN32
		if (!CreatePipe(&pipeRead, &pipeWrite, NULL, 0)) {
			pipeRead = pipeWrite = NULL;
			fprintf(stderr, "Warning: could not create pipes.\n");
		}
		#else
		if (socketpair(AF_UNIX, SOCK_STREAM, 0, sockets)) {
			sockets[0] = sockets[1] = -1;
			fprintf(stderr, "Warning: could not create socket pair.\n");
		}
		#endif
	}
	
	~LocalPipe() {
		#ifdef WIN32
		CloseHandle(pipeRead);
		CloseHandle(pipeWrite);
		#else
		close(sockets[0]);
		close(sockets[1]);
		#endif
	}
	
	int read(unsigned int size, void *data) {
		#ifdef WIN32
			DWORD numRead;
			if (!ReadFile(pipeRead, data, size, &numRead, NULL)) return -1;
			return numRead;
		#else
			return recv(sockets[0], data, size,0);
		#endif
	}

	int write(unsigned int size, void *data) {
		#ifdef WIN32
			DWORD numWritten;
			if (!WriteFile(pipeWrite, data, size, &numWritten, NULL)) return -1;
			return numWritten;
		#else
			return send(sockets[1], data, size,0);
		#endif
	}
	
	protected:
	#ifdef WIN32
		HANDLE pipeRead, pipeWrite;
	#else
		int sockets[2];
	#endif
};

#endif
