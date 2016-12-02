/** Simple C++ class for managing ASCII requests from a TCP port.

	(C) 2010 S. Klanke
*/

#include <stdio.h>
#include <stdint.h>
#include <string.h>

#include <StringServer.h>

#ifdef WIN32

#ifndef SD_BOTH
#define SD_BOTH   0x02
#endif

#include <windows.h>
#define socklen_t       int
#define shutdown_rw(s)  shutdown(s, SD_BOTH)

#else // Linux and OS X

#include <arpa/inet.h>
#include <sys/socket.h>
#include <sys/select.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>
#include <fcntl.h>

#define SOCKET          int
#define INVALID_SOCKET  -1
#define closesocket     close
#define shutdown_rw(s)  shutdown(s, SHUT_RDWR)

#endif

#define MAXCLIENTS 32

static int numStringServers = 0;

struct StringServerClientCtrl {
	SOCKET sock;
	std::string req;
	std::string resp;
};

struct StringServerCtrl {
	SOCKET sock;
};

StringServer::StringServer(int bufSize) {
	defaultBufSize = bufSize;
	server = new StringServerCtrl;
	client = new StringServerClientCtrl*[MAXCLIENTS];
	listening = false;
	numClients = 0;
	maxFD = 0;

	if (numStringServers++ == 0) {
		#ifdef WIN32
		WSADATA wsa = {0,0};
		if(WSAStartup(MAKEWORD(1, 1), &wsa)) {
			fprintf(stderr, "StringServer: cannot start WIN32 sockets.\n");
		}
		#endif
	}
}

StringServer::~StringServer() {
	stopListening();
	delete[] client;
	delete server;

	if (--numStringServers == 0) {
		#ifdef WIN32
		WSACleanup();
		#endif
	}
}

void StringServer::stopListening() {
	if (!listening) return;

	for (int i=0;i<numClients;i++) {
		shutdown_rw(client[i]->sock);
		closesocket(client[i]->sock);
		delete client[i];
		client[i] = NULL;
	}
	numClients = 0;
	closesocket(server->sock);
	listening = false;
}


bool StringServer::startListening(int port) {
	if (listening) return false;

	// setup TCP socket
	struct sockaddr_in sa;
	server->sock = socket(PF_INET, SOCK_STREAM, 0);
	if (server->sock == INVALID_SOCKET) {
		perror("StringServer::listen -> socket");
		return false;
	}
	// prevend "bind: address already in use"
	int optval = 1;
	if (setsockopt(server->sock, SOL_SOCKET, SO_REUSEADDR, (const char*)&optval, sizeof(optval)) < 0) {
		perror("StringServer::listen -> setsockopt");
		closesocket(server->sock);
		return false;
	}

	memset(&sa, 0, sizeof(sa));
	sa.sin_family = AF_INET;
	sa.sin_port   = htons(port);
	sa.sin_addr.s_addr = htonl(INADDR_ANY);
	if (bind(server->sock, (struct sockaddr *) &sa, sizeof(sa)) < 0) {
		perror("StringServer::listen -> bind");
		closesocket(server->sock);
		return false;
	}

	/* place the socket in non-blocking mode, required to do thread cancelation */
#ifdef WIN32
	{
		unsigned long enable = 0;
		ioctlsocket(server->sock, FIONBIO, &enable);
	}
#else
	optval = fcntl(server->sock, F_GETFL, NULL);
	optval = optval | O_NONBLOCK;
	if (fcntl(server->sock, F_SETFL, optval)<0) {
		perror("StringServer::listen -> fcntl");
		closesocket(server->sock);
		return false;
	}
	maxFD = server->sock;
#endif

	if (listen(server->sock, 5)<0) {
		perror("StringServer::listen -> listen");
		closesocket(server->sock);
		return false;
	}

	listening = true;

	return true;
}


int StringServer::checkRequests(StringRequestHandler& handler, int milliSeconds) {
	fd_set readSet;
	fd_set writeSet;
	struct timeval tv;
	int requests = 0;

	if (!listening) return -1;

	tv.tv_sec = milliSeconds / 1000;
	tv.tv_usec = (milliSeconds % 1000) * 1000;

	FD_ZERO(&readSet);
	FD_ZERO(&writeSet);

	if (numClients < MAXCLIENTS) {
		FD_SET(server->sock, &readSet);
	}
	for (int i=0;i<numClients;i++) {
		FD_SET(client[i]->sock, &readSet);
		if (client[i]->resp.size() > 0) FD_SET(client[i]->sock, &writeSet);
	}

	int n = select(maxFD+1, &readSet, &writeSet, NULL, &tv);

	if (n==0) return 0;

	// first iterate through client connections
	int idx = 0;
	while (idx<numClients) {
		StringServerClientCtrl *cli = client[idx];
		bool startResponse = cli->resp.empty();

		if (FD_ISSET(cli->sock, &readSet)) {
			int r = recv(cli->sock, rcvBuf, sizeof(rcvBuf), 0);
			if (r<=0) {
				// printf("Closing client %i\n", idx);
				closeClient(idx);
				continue;
			}

			cli->req.append(rcvBuf, r);
			// startResponse = "empty before" and "non-empty now"
			startResponse &= checkCompletion(handler, cli);
		}

		if (FD_ISSET(cli->sock, &writeSet) || startResponse) {
			int r = send(cli->sock, cli->resp.data(), cli->resp.size(), 0);
			if (r>0) {
				cli->resp.erase(0,r);
			}
		}

		++idx;
	}

	// accept new connection, if any
	if (FD_ISSET(server->sock, &readSet)) {
		struct sockaddr_in sa;
		socklen_t size_sa = sizeof(sa);

		SOCKET c = accept(server->sock, (struct sockaddr *)&sa, &size_sa);

		if (c == INVALID_SOCKET) {
			perror("StringServer::checkRequests -> accept");
		} else {
			client[numClients] = new StringServerClientCtrl;
			client[numClients]->req.reserve(defaultBufSize);
			client[numClients]->sock = c;
			numClients++;

			#ifndef WIN32
			if (c > maxFD) maxFD = c;
			#endif
		}
	}

	return requests;
}

bool StringServer::checkCompletion(StringRequestHandler& handler, StringServerClientCtrl *cli) {
	bool newResp = false;
	while (1) {
		size_t p0 = cli->req.find('\n');
		if (p0 == cli->req.npos) return newResp;

		std::string out(cli->req, 0, p0);
		std::string in = handler.handleStringRequest(out);
		cli->req.erase(0, p0+1);
		if (!in.empty()) {
			newResp = true;
			cli->resp.append(in);
		}
	}
}

void StringServer::closeClient(int idx) {
	shutdown_rw(client[idx]->sock);
	closesocket(client[idx]->sock);
	delete client[idx];

	if (idx < numClients - 1) {
		client[idx] = client[numClients-1];
	} else {
		client[idx] = NULL;
	}
	numClients--;
}
