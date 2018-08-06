/*
 * Copyright (C) 2010, Stefan Klanke
 * 	Modified by Tim van Mourik 2015
 * Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
 * Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 */

#include <stdio.h>
#include <math.h>

#include <PixelDataGrabber.h>

#include <Fl/Fl.H>
#include <Fl/Fl_Window.H>
#include <Fl/Fl_Input.H>
#include <Fl/Fl_Int_Input.H>
#include <Fl/Fl_Browser.H>
#include <Fl/Fl_Button.H>
#include <Fl/Fl_Box.H>

Fl_Window *window;
Fl_Input *inpHostname, *inpDirectory, *inpUdpTargetHost;
Fl_Int_Input *inpPort, *inpUdpTargetPort;
Fl_Button *butConnect, *butListen, *butUdpEnable;
Fl_Browser *msgBrowser;
Fl_Box *piBox;
PixelDataGrabber pdg;
char hostname[256];
int port;
char directory[256];

SOCKET udpSocket;
char udpTargetHost[256];
int udpTargetPort;


SOCKET createSocketUDP(const char *address, int port) {
	struct hostent *host;
	struct sockaddr_in addr;
	unsigned int one = 1;
	SOCKET sock;
	
	host = gethostbyname(address);
	
	if (host == NULL || host->h_length == 0) {
		fprintf(stderr, "DNS look up failed on '%s'\n", address);
		return INVALID_SOCKET;
	}


	memset(&addr, 0, sizeof(struct sockaddr_in));
	addr.sin_family = AF_INET;
	addr.sin_port = htons(port);
	memcpy(&(addr.sin_addr.s_addr), host->h_addr_list[0], sizeof(addr.sin_addr.s_addr));

	sock = socket(PF_INET, SOCK_DGRAM, IPPROTO_UDP);
	if (sock == INVALID_SOCKET) {
		perror("socket(...): ");
		return INVALID_SOCKET;
	}
	setsockopt(sock, SOL_SOCKET, SO_BROADCAST, (const char *) &one, sizeof(one));
	if (connect(sock, (struct sockaddr *) &addr, sizeof(struct sockaddr_in))) {
		perror("connect(...): ");
		closesocket(sock);
		return INVALID_SOCKET;
	}

	return sock;
}


void updatePiBox() {
	char msg[200];
	int ne = pdg.getNumEchos();
	int nt = pdg.getNumWritten();
	
	snprintf(msg, 200, "Format: %i x %i x %i  -  Echos: %i  -  Scans: %i", 
					pdg.getReadoutResolution(), pdg.getPhaseResolution(), 
					pdg.getNumSlices(), ne, nt);
	piBox->copy_label(msg);
	piBox->redraw();
}

void addTimedMsg(const char *s) {
	char msg[300];
	SYSTEMTIME sysT;
	
	GetLocalTime(&sysT);
	snprintf(msg, 300, "%02u:%02u:%02u %s", sysT.wHour, sysT.wMinute, sysT.wSecond, s);
	msgBrowser->add(msg);
	msgBrowser->bottomline(msgBrowser->size());
}


void changeHostPortDisplay(bool connected) {
	Fl_Color col = connected ? FL_DARK_GREEN : FL_RED;
	const char *lab = connected ? "Disconnect" : "Connect";
	int ro = connected ? 1 : 0;
	
	inpHostname->readonly(ro);
	inpPort->readonly(ro);
	butConnect->label(lab);
	inpHostname->textcolor(col);
	inpHostname->redraw();
	inpPort->textcolor(col);
	inpPort->redraw();
}


void changeUdpDisplay(bool connected) {
	Fl_Color col = connected ? FL_DARK_GREEN : FL_RED;
	const char *lab = connected ? "Disable" : "Enable";
	int ro = connected ? 1 : 0;
	
	inpUdpTargetHost->readonly(ro);
	inpUdpTargetPort->readonly(ro);
	butUdpEnable->label(lab);
	inpUdpTargetHost->textcolor(col);
	inpUdpTargetHost->redraw();
	inpUdpTargetPort->textcolor(col);
	inpUdpTargetPort->redraw();
}

void disconnect() {
	pdg.connectToFieldTrip(NULL,0);
	addTimedMsg("Disconnected from FieldTrip");
	changeHostPortDisplay(false);
}	


bool connect() {
	strncpy(hostname, inpHostname->value(), 256);
	port = atoi(inpPort->value());
	if (pdg.connectToFieldTrip(hostname, port)) {
		addTimedMsg("Connected to FieldTrip");
		changeHostPortDisplay(true);
		updatePiBox();
		return true;
	} else {
		addTimedMsg("Could not connect");
		return false;
	}
}

void connectCallback(Fl_Widget *widget) {
	if (pdg.isConnected()) {
		disconnect();
	} else {
		connect();
	}
}


void udpCallback(Fl_Widget *widget) {
	if (udpSocket != INVALID_SOCKET) {
		closesocket(udpSocket);
		udpSocket = INVALID_SOCKET;
		addTimedMsg("Disabled UDP messages");
		changeUdpDisplay(false);
	} else {
		strncpy(udpTargetHost, inpUdpTargetHost->value(), 256);
		udpTargetPort = atoi(inpUdpTargetPort->value());
		udpSocket = createSocketUDP(udpTargetHost, udpTargetPort);
		if (udpSocket != INVALID_SOCKET) {
			addTimedMsg("Set UDP message target");
			changeUdpDisplay(true);
		} else {
			addTimedMsg("Could not set UDP target");
			changeUdpDisplay(false);
		}
	} 
}

bool startMonitor() {
    strncpy(directory, inpDirectory->value(), 256);
    if (pdg.monitorDirectory(directory))
    {
		addTimedMsg("Starting to monitor");
		butListen->label("Stop");
		inpDirectory->textcolor(FL_DARK_GREEN);
		inpDirectory->readonly(1);
		inpDirectory->redraw();
		return true;
	}
    else
    {
        addTimedMsg("Could not start monitoring");
        return false;
    }
}

void stopMonitor() {
	pdg.monitorDirectory(NULL);
	addTimedMsg("Stopped monitoring");
	butListen->label("Start");
	inpDirectory->readonly(0);
	inpDirectory->textcolor(FL_RED);	
	inpDirectory->redraw();
}

void listenCallback(Fl_Widget *widget) {
	if (pdg.isListening()) {
		stopMonitor();
	} else {
		startMonitor();
	}
}

int main(int argc, char *argv[]) {
	bool autoLogin;
	char portAsChar[8];
	WSADATA wsa;

	if (WSAStartup(MAKEWORD(1, 1), &wsa)) {
		fprintf(stderr, "WSAStartup failed!\n");
		exit(1);
	}	
	
	udpSocket = INVALID_SOCKET;
	
	timeBeginPeriod(1);
	
	pdg.setVerbosity(3);
	
	if (argc>=2) {
		strncpy(hostname, argv[1], 256);
		autoLogin = true;
	} else {
		strncpy(hostname, "localhost", 256);
		autoLogin = false;
	}
	
	if (argc>=3) {
		port = atoi(argv[2]);
	} else {
		port = 1972;
	}
	snprintf(portAsChar,7, "%i",port);
	
	if (argc>=4) {
		strncpy(directory, argv[3], 256);
	} else {
		strncpy(directory, "E:\\IMAGE", 256);
	}
	
	Fl::visual(FL_RGB);
	window = new Fl_Window(100,100,400,350,"Realtime fMRI streamer, (C) Stefan Klanke & Tim van Mourik");
	inpHostname = new Fl_Input(20,30,200,25,"Hostname");
	inpPort = new Fl_Int_Input(230,30,60,25,"Port");
	butConnect = new Fl_Button(300,30,80,25,"Connect");
	inpDirectory = new Fl_Input(20,80,270,25,"Directory");
	butListen = new Fl_Button(300,80,80,25,"Start");
	piBox = new Fl_Box(20, 120, 360, 25);
	msgBrowser = new Fl_Browser(20, 160, 360, 120);
	inpUdpTargetHost = new Fl_Input(20,300,200,25,"UDP target hostname");
	inpUdpTargetPort = new Fl_Int_Input(230,300,60,25,"UDP port");
	butUdpEnable = new Fl_Button(300,300,80,25,"Enable");
	
	inpHostname->align(FL_ALIGN_TOP);
	inpPort->align(FL_ALIGN_TOP);
	inpHostname->value(hostname);
	inpPort->value(portAsChar);
	inpDirectory->align(FL_ALIGN_TOP);
	inpDirectory->value(directory);
	butConnect->callback(connectCallback);
	butListen->callback(listenCallback);
	inpUdpTargetHost->value("lab-mri001");
	inpUdpTargetHost->align(FL_ALIGN_TOP);
	inpUdpTargetPort->value("1990");
	inpUdpTargetPort->align(FL_ALIGN_TOP);
	butUdpEnable->callback(udpCallback);
	
	piBox->box(FL_DOWN_BOX);
	
	window->end();
	window->show();
	
	if (autoLogin) {
		connect();
	} else {
		changeHostPortDisplay(false);
	}
	startMonitor();
	
	while (Fl::check()) {
		Fl::wait(0);
		if (pdg.run(20) < 0) {
			Fl::wait();
		} else {
			switch(pdg.getLastAction()) {
				case PixelDataGrabber::Nothing:
					break;
				case PixelDataGrabber::OutOfMemory:
					addTimedMsg("Out of memory!");
					break;
				case PixelDataGrabber::ProtocolRead:
					addTimedMsg("Read protocol");
					if (udpSocket != INVALID_SOCKET) {
						int n = send(udpSocket, "RESET", 5, 0);
						printf("Send out %i bytes over UDP\n", n);
					}
					updatePiBox();
					break;
				case PixelDataGrabber::PixelsTransmitted:			
					if (pdg.getNumWritten()==1) addTimedMsg("Transmitted first sample");
					updatePiBox();
					break;
				case PixelDataGrabber::BadPixelData:
					addTimedMsg("Pixel data does not match protocol");
					break;
				case PixelDataGrabber::TransmissionError:
					addTimedMsg("Communication error - disconnecting");
					disconnect();
					break;
			}
		}
	}
	
	delete window;
	
	WSACleanup();
	
	return 0;
}
