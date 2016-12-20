/** Simple C++ class for managing ASCII requests from a TCP port.

	(C) 2010 S. Klanke
*/

#ifndef __StringServer_h
#define __StringServer_h

#include <string>
#include <ctype.h>

class StringRequestHandler {
	public:
	virtual std::string handleStringRequest(const std::string& request)=0;
};

struct StringServerClientCtrl;
struct StringServerCtrl;

class StringServer {
	public:

	StringServer(int bufSize = 2048);
	~StringServer();

	bool startListening(int port);
	void stopListening();

	int checkRequests(StringRequestHandler& handler, int milliSeconds = 0);

	static std::string getNextToken(const std::string& in, unsigned int& pos) {
		while (isspace(in[pos])) {
			if (++pos == in.size()) {
				return std::string();
			}
		}
		unsigned int p0 = pos;
		while (!isspace(in[pos])) {
			if (++pos == in.size()) break;
		}
		return in.substr(p0, pos-p0);
	}

	protected:

	void closeClient(int idx);
	bool checkCompletion(StringRequestHandler& handler, StringServerClientCtrl *cli);

	int defaultBufSize;
	bool listening;
	int numClients;
	int maxFD;
	StringServerClientCtrl **client;
	StringServerCtrl *server;
	char rcvBuf[16*1024];
};

#endif
