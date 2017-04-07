/*
 * Copyright (C) 2010, Stefan Klanke
 * Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
 * Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 */

#include <FtBuffer.h>

#ifdef WIN32
WSADATA FtConnection::wsa = {0,0};
#endif
int FtConnection::numConnections = 0;


FtConnection::~FtConnection() {
	disconnect();
	if (--numConnections == 0) {
		#ifdef WIN32
		if (wsa.wVersion > 0) {
			WSACleanup();
			wsa.wVersion = 0;
		}
		#endif
	}
}

bool FtConnection::connectTcp(const char *hostname, int port) {
	struct sockaddr_in sa;
	struct hostent *host;

	disconnect();

#ifdef WIN32
	if (wsa.wVersion == 0) {
		// We only need to do this once
		if(WSAStartup(MAKEWORD(1, 1), &wsa)) {
			fprintf(stderr, "open_connection: cannot start sockets\n");
			return false;
		}
	}
#endif

	if ((host = gethostbyname(hostname)) == NULL) {
		fprintf(stderr, "connectTcp: nslookup1 failed on '%s'\n", hostname);
		return false;
	}

	if (host->h_length == 0) {
		fprintf(stderr, "connectTcp: nslookup2 failed on '%s'\n", hostname);
		return false;
	}

	bzero(&sa, sizeof(sa));
	sa.sin_family = AF_INET;
	sa.sin_port = htons(port);
	memcpy(&(sa.sin_addr.s_addr), host->h_addr_list[0], sizeof(sa.sin_addr.s_addr));

	sock = socket(PF_INET, SOCK_STREAM, 0);
	if (sock < 0) {
		perror("connectTcp: socket");
		return false;
	}

	for (int r = 0;r<=retry;r++) {
		if (::connect(sock, (struct sockaddr *)&sa, sizeof(sa))==0) {
			type = 1;
			return true;
		}
		if (r==retry) break;
		/* wait 5 miliseconds and try again */
		perror("connectTcp:connect");
		usleep(5000);
	}
	closesocket(sock);
	sock = -1;
	return false;
}

bool FtConnection::connectUnix(const char *pathname) {
	disconnect();
#ifndef WIN32
	struct sockaddr_un sa;

	bzero(&sa, sizeof(sa));
	sa.sun_family = AF_UNIX;
	strncpy(sa.sun_path, pathname, sizeof(sa.sun_path));

	sock = socket(AF_UNIX, SOCK_STREAM, 0);
	if (sock < 0) {
		perror("connectUnix: socket");
		return false;
	}

	for (int r = 0;r<=retry;r++) {
		if (::connect(sock, (struct sockaddr *)&sa, sizeof(sa))==0) {
			type = 2;
			return true;
		}
		if (r==retry) break;
		/* wait 5 miliseconds and try again */
		perror("connectUnix:connect");
		usleep(5000);
	}
	closesocket(sock);
	sock = -1;
#endif
	return false;
}

bool FtConnection::connect(const char *address) {
	const char *colPos = strchr(address, ':');
	if (colPos != NULL) {
		int len = colPos - address;
		char *hostname = new char[len+1];
		memcpy(hostname, address, len);
		hostname[len] = 0;

		int port = atoi(colPos+1);
		if (port == 0) return false;
		bool result = connectTcp(hostname, port);
		delete[] hostname;
		return result;
	} else {
		return connectUnix(address);
	}
}
