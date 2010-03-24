/*
 * Copyright (C) 2010, Stefan Klanke
 * Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
 * Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 */
#include <Fl/Fl.h>
#include <Fl/Fl_Window.h>
#include <Fl/Fl_Input.h>
#include <Fl/Fl_Int_Input.h>
#include <Fl/Fl_Browser.h>
#include <Fl/Fl_Button.h>
#include <Fl/Fl_Box.h>
#include <stdio.h>
#include <PixelDataGrabber.h>
#include <math.h>

Fl_Window *window;
Fl_Input *inpHostname, *inpDirectory;
Fl_Int_Input *inpPort;
Fl_Button *butConnect, *butListen;
Fl_Browser *msgBrowser;
Fl_Box *piBox;
PixelDataGrabber pdg;
char hostname[256];
int port;
char directory[256];


void updatePiBox() {
	char msg[200];
	snprintf(msg, 200, "Format: %i x %i x %i  -  Scans: %4i", 
				pdg.getReadoutResolution(), pdg.getPhaseResolution(), 
				pdg.getNumSlices(), pdg.getNumScansWritten());
	piBox->copy_label(msg);
	piBox->redraw();
}

void addTimedMsg(const char *s) {
	char msg[300];
	SYSTEMTIME sysT;
	
	GetSystemTime(&sysT);
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

bool startMonitor() {
	strncpy(directory, inpDirectory->value(), 256);
	if (pdg.monitorDirectory(directory)) {
		addTimedMsg("Starting to monitor");
		butListen->label("Stop");
		inpDirectory->textcolor(FL_DARK_GREEN);
		inpDirectory->readonly(1);
		inpDirectory->redraw();
		return true;
	}
	addTimedMsg("Could not start monitoring");
	return false;
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
	
	
	timeBeginPeriod(1);
	
	pdg.setVerbosity(0);
	
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
	window = new Fl_Window(100,100,400,300,"Realtime fMRI streamer");
	inpHostname = new Fl_Input(20,30,200,25,"Hostname");
	inpPort = new Fl_Int_Input(230,30,60,25,"Port");
	butConnect = new Fl_Button(300,30,80,25,"Connect");
	inpDirectory = new Fl_Input(20,80,270,25,"Directory");
	butListen = new Fl_Button(300,80,80,25,"Start");
	piBox = new Fl_Box(20, 120, 360, 25);
	msgBrowser = new Fl_Browser(20, 160, 360, 120);
	
	inpHostname->align(FL_ALIGN_TOP);
	inpPort->align(FL_ALIGN_TOP);
	inpHostname->value(hostname);
	inpPort->value(portAsChar);
	inpDirectory->align(FL_ALIGN_TOP);
	inpDirectory->value(directory);
	butConnect->callback(connectCallback);
	butListen->callback(listenCallback);
	
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
		if (pdg.run(200) < 0) {
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
					updatePiBox();
					break;
				case PixelDataGrabber::PixelsTransmitted:			
					if (pdg.getNumScansWritten()==1) addTimedMsg("Transmitted first scan");
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
	
	return 0;
}
