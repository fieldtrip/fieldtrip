/*
 * Graphical signal generator that streams data to online buffer.
 *
 * Copyright (C) 2010, Stefan Klanke
 * Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
 * Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 */

#include <Fl/Fl.H>
#include <Fl/Fl_Window.H>
#include <Fl/Fl_Input.H>
#include <Fl/Fl_Int_Input.H>
#include <Fl/Fl_Float_Input.H>
#include <Fl/Fl_Button.H>
#include <Fl/Fl_Box.H>
#include <Fl/Fl_Value_Slider.H>

#include <stdio.h>
#include <math.h>
#include <buffer.h>
#include <FtBuffer.h>
#include <Clock.h>

#define MAX_CHAN  300
#define MAX_BLOCK 10000

Clock T;
Fl_Window *window;
Fl_Input *inpHostname;
Fl_Int_Input *inpPort, *inpNumChannels, *inpBlockSize;
Fl_Float_Input *inpFSample;
Fl_Button *butStartStop;
Fl_Box *infoBox;
Fl_Value_Slider *sliFreq, *sliAmp;
char hostname[256];
int port, numChannels, blockSize;
float fSample;
int blocksWritten;
int ftSocket = -1;
int curSample;
bool active = false;
float phase = 0.0;
float amplitude, frequency;
float samples[MAX_CHAN*MAX_BLOCK];
double tBlock;

void updateInfoBox() {
	char msg[200];
	snprintf(msg, 200, "Blocks / samples written: %i / %i", blocksWritten, blocksWritten*blockSize);
	infoBox->copy_label(msg);
	infoBox->redraw();
}

void changeWidgets(bool connected) {
	const char *lab = connected ? "Stop" : "Start";
	int ro = connected ? 1 : 0;
	Fl_Color col = connected ? FL_GREEN : FL_BLACK;

	inpHostname->readonly(ro);
	inpPort->readonly(ro);
	inpNumChannels->readonly(ro);
	inpBlockSize->readonly(ro);
	inpFSample->readonly(ro);

	inpHostname->textcolor(col);
	inpPort->textcolor(col);
	inpNumChannels->textcolor(col);
	inpBlockSize->textcolor(col);
	inpFSample->textcolor(col);

	inpHostname->redraw();
	inpPort->redraw();
	inpNumChannels->redraw();
	inpBlockSize->redraw();
	inpFSample->redraw();

	butStartStop->label(lab);
	butStartStop->redraw();
}

void startStopCallback(Fl_Widget *widget) {
	char lab[20];

	if (active) {
		close_connection(ftSocket);
		infoBox->label("Not connected\n");
		infoBox->redraw();

		ftSocket = -1;
		active = false;
	} else {
		port = atoi(inpPort->value());
		strncpy(hostname, inpHostname->value(), sizeof(hostname));
		hostname[255] = 0;

		blockSize = atoi(inpBlockSize->value());
		if (blockSize < 1) blockSize = 1;
		if (blockSize > MAX_BLOCK) blockSize = MAX_BLOCK;

		snprintf(lab, 20, "%i", blockSize);
		inpBlockSize->value(lab);
		inpBlockSize->redraw();

		numChannels = atoi(inpNumChannels->value());
		if (numChannels < 1) numChannels = 1;
		if (numChannels > MAX_CHAN) numChannels = MAX_CHAN;

		snprintf(lab, 20, "%i", numChannels);
		inpNumChannels->value(lab);
		inpNumChannels->redraw();

		fSample = atof(inpFSample->value());
		if (fSample < 1) fSample = 1;
		if (fSample > 100000) fSample = 100000;
		// sliFreq->range(0.01, 0.5*fSample);

		snprintf(lab, 20, "%.2f", fSample);
		inpFSample->value(lab);
		inpFSample->redraw();


		ftSocket = open_connection(hostname, port);
		if (ftSocket == -1) {
			infoBox->label("Could not connect\n");
			infoBox->redraw();
		} else {
			FtBufferRequest req;
			FtBufferResponse resp;
			int r;

			req.prepPutHeader(numChannels, DATATYPE_FLOAT32, fSample);

			r = tcprequest(ftSocket, req.out(), resp.in());
			if (r<0 || !resp.checkPut()) {
				infoBox->label("Could not write header\n");
				infoBox->redraw();
			} else {
				infoBox->label("Connected + wrote header\n");
				infoBox->redraw();
				T.reset();
				blocksWritten = 0;
				phase = 0.0;
				curSample = 0;
				tBlock = blockSize / fSample;
				active = true;
			}
		}
	}

	changeWidgets(active);
}

void addSample() {
	if (curSample >= blockSize) return;

	float v = sinf(phase)*amplitude;
	phase += (frequency/fSample)*2.0*M_PI;
	if (phase > 2*M_PI) phase-=2*M_PI;

	float *s = samples + curSample * numChannels;

	for (int i=0;i<numChannels;i++) {
		s[i] = v;
	}

	++curSample;
}


int main(int argc, char *argv[]) {
	double t, tb;
	FtBufferRequest req;
	FtBufferResponse resp;


	Fl::visual(FL_RGB);
	window = new Fl_Window(100,100,300,300,"Sinewave to FieldTrip buffer");
	inpHostname = new Fl_Input(20,25,190,25,"Hostname");
	inpHostname->align(FL_ALIGN_TOP);
	inpHostname->value("localhost");
	inpPort = new Fl_Int_Input(220,25,60,25,"Port");
	inpPort->align(FL_ALIGN_TOP);
	inpPort->value("1972");

	inpNumChannels = new Fl_Int_Input(20,75,80,25,"#Channels");
	inpNumChannels->align(FL_ALIGN_TOP);
	inpNumChannels->value("16");

	inpBlockSize = new Fl_Int_Input(110,75,80,25,"Block size");
	inpBlockSize->align(FL_ALIGN_TOP);
	inpBlockSize->value("32");

	inpFSample = new Fl_Float_Input(200,75,80,25,"Sampl.freq");
	inpFSample->align(FL_ALIGN_TOP);
	inpFSample->value("256");

	butStartStop = new Fl_Button(100,110,100,25,"Start");
	butStartStop->callback(startStopCallback);


	sliFreq = new Fl_Value_Slider(20, 170, 260, 25, "Signal frequency");
	sliFreq->type(FL_HOR_NICE_SLIDER);
	sliFreq->align(FL_ALIGN_TOP);
	sliFreq->range(0.01, 50.0);
	sliFreq->value(1.0);

	sliAmp  = new Fl_Value_Slider(20, 220, 260, 25, "Signal amplitude");
	sliAmp->type(FL_HOR_NICE_SLIDER);
	sliAmp->align(FL_ALIGN_TOP);
	sliAmp->value(1.0);

	infoBox = new Fl_Box(5, 270, 290, 25);
	infoBox->label("Not connected");
	infoBox->box(FL_THIN_DOWN_BOX);


	window->end();
	window->show();

	tb = 0.0;

	while (Fl::check()) {
		Fl::wait(0.01);
		if (active) {
			amplitude = sliAmp->value();
			frequency = sliFreq->value();

			// add samples roughly corresponding to current time
			t = T.getRel();
			int n = (int) ((t - tb)/tBlock);

			for (int i=0;i<n;i++) addSample();

			// enough time passed for one block ? write it
			if (t >= tBlock*(1+blocksWritten)) {
				tb -= tBlock;

				// add remaing samples
				n = blockSize - curSample;
				for (int i=0;i<n;i++) addSample();

				bool ok = req.prepPutData(numChannels, blockSize, DATATYPE_FLOAT32, samples);

				if (!ok) {
					infoBox->label("Out of memory\n");
					infoBox->redraw();
				} else {
					int r = tcprequest(ftSocket, req.out(), resp.in());
					if (r < 0 || !resp.checkPut()) {
						startStopCallback(NULL); // this will disconnect
					} else {
						++blocksWritten;
						curSample = 0;
						updateInfoBox();
					}
				}
			}
		}
	}

	if (ftSocket != -1) close_connection(ftSocket);

	delete window;

	return 0;
}
