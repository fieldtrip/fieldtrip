/** Simple C++ class for managing ASCII requests from a TCP port.
	(C) 2010 S. Klanke
*/
#ifndef __StringServer_h
#define __StringServer_h

#include <string>

class StringRequestHandler {
	public:
	virtual void handleStringRequest(const std::string& request, std::string& response)=0;
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
