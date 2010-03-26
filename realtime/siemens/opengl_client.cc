/*
 * Copyright (C) 2010, Stefan Klanke
 * Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
 * Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 */
#include <Fl/Fl.h>

#include <Fl/Fl_Gl_Window.h>

#include <GL/gl.h>
#include <GL/glu.h>

#include <Brain3dWindow.h>

#include <Fl/Fl_Value_Slider.h>

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <buffer.h>
#include <siemensap.h>
#include <SimpleStorage.h>
#include <FtBuffer.h>


class PixelData2Texture {
	public:
	
	PixelData2Texture() {
		numAlloc = 0;
		NS = 0;
		NT = 0;
		image = 0;
	}
	
	~PixelData2Texture() {
		if (image) free(image);
		if (NT>0) {
			glDeleteTextures(NT, texture);
		}
	}
		
	void getFromBuffer(void *pixelData, int w, int h, int ns) {
		uint16_t *pixels;
		int numPixels = w*h*ns;
		
		if (NT < ns) {
			for (int i=NT;i<ns;i++) {
				glGenTextures(1, &texture[i]);
				int err = glGetError();
				
				if (err!=GL_NO_ERROR) printf("Error in glGenTextures %i -> %i!\n",i,err);
			}
			NT = ns;
		}
		
		printf("To tonvert: %i\n", numPixels);
		if (numPixels < 1) return;
		
		pixels = (uint16_t *) pixelData;
				
		if (numPixels > numAlloc) {
			if (image != 0) free(image);
			image = (unsigned char *) malloc(2*numPixels);
			if (image == 0) {
				numAlloc = 0;
				W = H = NS = 0;
				return;
			}
			numAlloc = numPixels;
			W = w; 
			H = h;
			NS = ns;
		}
		
		for (int i=0;i<numPixels;i++) {
			int v32 = pixels[i]*255;
			image[2*i] = v32 / 1500; /* use proper scaling sometime */
			image[2*i+1] = v32 / 1500;
		}
		
		for (int i=0;i<ns; i++) {
			// Generate The Texture
			glBindTexture(GL_TEXTURE_2D, texture[i]);
			if (glGetError()!=GL_NO_ERROR) printf("Error in bind %i!\n",i);
			
			glPixelTransferf(GL_RED_BIAS, 0.0);
			glPixelTransferf(GL_GREEN_BIAS, 0.0);
			glPixelTransferf(GL_BLUE_BIAS, 0.0);
			
			glTexImage2D(GL_TEXTURE_2D, 0, 4, w, h, 0, GL_LUMINANCE_ALPHA, GL_UNSIGNED_BYTE, image + i*w*h*2);
			if (glGetError()!=GL_NO_ERROR) printf("Error in texImage %i!\n",i);
			glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);	// Linear Filtering
			glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);	// Linear Filtering
		}
		
		printf("Converted: %i\n", numPixels);
	}
			
	const unsigned char *getImage() const {	return image; }
	int getW() const { return W; }
	int getH() const {	return H; }
	int getNumSlices() const {	return NS;	}	
	const GLuint *getTextures() const {	return texture;	}
	
	protected:
	
	GLuint texture[64];	
	
	unsigned char *image;
	int numAlloc;
	int W, H, NS, NT;
};



Brain3dWindow *BW;
Fl_Value_Slider *sliFactor, *sliOffset;
PixelData2Texture px2tex;
unsigned int prevSamples = 0;
int ftbSocket = -1;
int readoutResolution = 0, phaseResolution = 0, numSlices = 0;
double phaseFOV = 0.0, readoutFOV = 0.0;


void getProtInfo(const char *ap, unsigned int length) {
	const sap_item_t *item;
	sap_item_t *PI = sap_parse(ap, length);
	
	item = sap_search_deep(PI, "sKSpace.lBaseResolution");
	if (item!=NULL && item->type == SAP_LONG) {
		long res = *((long *) item->value);
		readoutResolution = (res > 0) ? res : 0;
	}
	
	item = sap_search_deep(PI, "sSliceArray.lSize");
	if (item!=NULL && item->type == SAP_LONG) {
		long slices = *((long *) item->value);
		numSlices = (slices > 0) ? slices : 0;
	}
	
	item = sap_search_deep(PI, "sSliceArray.asSlice[0].dPhaseFOV");
	if (item!=NULL && item->type == SAP_DOUBLE) {
		phaseFOV = *((double *) item->value);
		printf("PhaseFOV: %f\n", phaseFOV);
	}
	
	item = sap_search_deep(PI, "sSliceArray.asSlice[0].dReadoutFOV");
	if (item!=NULL && item->type == SAP_DOUBLE) {
		readoutFOV = *((double *) item->value);
		printf("ReadoutFOV: %f\n", readoutFOV);
	}
	
	if (phaseFOV > 0.0 && readoutFOV > 0.0) {
		phaseResolution = (unsigned int) round(readoutResolution * phaseFOV / readoutFOV);
	} else {
		phaseResolution = 0;
	}
	printf("Resolution: %i x %i x %i\n", readoutResolution, phaseResolution, numSlices);
	sap_destroy(PI);
}

bool readHeader() {
	SimpleStorage protBuffer;
	headerdef_t header_def;
	FtBufferRequest request;
	FtBufferResponse response;
	
	request.prepGetHeader();
	
	if (tcprequest(ftbSocket, request.out(), response.in()) < 0) {
		fprintf(stderr, "Error in communication. Buffer server aborted??\n");
		return false;
	}
	
	if (!response.checkGetHeader(header_def, &protBuffer)) {
		fprintf(stderr, "Error in received packet.\n");
		return false;
	}
	
	getProtInfo((char *) protBuffer.data(), protBuffer.size());
		
	if (header_def.data_type != DATATYPE_UINT16) {
		fprintf(stderr, "Data type != uint16\n");
		return false;
	}
	
	printf("GET_HDR: samples / channels: %i / %i\n", header_def.nsamples, header_def.nchans);
	return true;
}


void idleCall(void *dummy) {
	SimpleStorage pixBuffer;
	datadef_t data_def;
	FtBufferRequest request;
	FtBufferResponse response;
	unsigned newSamples;
	
	if (numSlices == 0) {
		if (!readHeader()) return;
	}
	
	request.prepWaitData(prevSamples, 50);
	
	if (tcprequest(ftbSocket, request.out(), response.in()) < 0) {
		fprintf(stderr, "Error in communication. Buffer server aborted??\n");
		return;
	}
	if (!response.checkWait(newSamples)) {
		fprintf(stderr, "Error in received packet.\n");
		return;
	}
	
	if (newSamples == prevSamples) return; // nothing new
	if (newSamples < prevSamples) {
		// oops ? do we have a new header?
		if (!readHeader()) return;
	}
	prevSamples = newSamples;
	
	if (prevSamples == 0) return;
	
	request.prepGetData(prevSamples-1, prevSamples-1);
	
	if (tcprequest(ftbSocket, request.out(), response.in()) < 0) {
		fprintf(stderr, "Error in communication. Buffer server aborted??\n");
		return;
	}
	if (!response.checkGetData(data_def, &pixBuffer)) {
		fprintf(stderr, "Error in received packet.\n");
		return;
	}
	
	px2tex.getFromBuffer(pixBuffer.data(), readoutResolution, phaseResolution, numSlices);
	BW->setSliceTextures(px2tex.getNumSlices(), px2tex.getTextures());
	BW->redraw();
}


void sliderCallback(Fl_Widget *widget) {
	BW->setWarp(sliFactor->value(), sliOffset->value());
	BW->redraw();
}

int main(int argc, char *argv[]) {
	const char *hostname = "localhost";
	int port = 1972;

	if (argc>1) hostname = argv[1];
	if (argc>2) port = atoi(argv[2]);
	
	printf("Trying to connect to fieldtrip buffer at %s:%d\n",hostname,port);
	
	ftbSocket = open_connection(hostname, port);
	if (ftbSocket<=0) {
		fprintf(stderr,"Failed\n");
		exit(1);
	}
   
   printf("Ok.\n");
		
	Fl::visual(FL_RGB);
	Fl_Window *window = new Fl_Window(100,100,600,700,"fMRI 3D client");

	BW = new Brain3dWindow(0,0,600,600);
	
	sliFactor = new Fl_Value_Slider(20,615,280,25);
	sliOffset = new Fl_Value_Slider(20,660,280,25);
	sliFactor->type(FL_HORIZONTAL);
	sliOffset->type(FL_HORIZONTAL);
	
	sliFactor->label("Z-axis distortion");
	sliOffset->label("Z-focus at");
	
	sliFactor->callback(sliderCallback);
	sliOffset->callback(sliderCallback);
	
	sliFactor->step(0.001);
	sliOffset->step(0.001);
	
	window->resizable(BW);
	window->end();
	
	window->show();
	
	BW->show();
	BW->setFanSize(0);
	
	Fl::wait(0);
	
	Fl::add_idle(idleCall);
	
	BW->redraw();
	
	Fl::run();
	//delete window;
	delete BW;
	
	printf("Closing connection...\n");
	close_connection(ftbSocket);
}
