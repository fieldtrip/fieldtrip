/*
 * Copyright (C) 2010, Stefan Klanke
 * Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
 * Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <platform.h>
#if defined (PLATFORM_OSX)
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include <FL/Fl.H>
#include <FL/Fl_Gl_Window.H>
#include <FL/Fl_Value_Slider.H>

#include <Brain3dWindow.h>

#include <buffer.h>
#include <siemensap.h>
#include <SimpleStorage.h>
#include <FtBuffer.h>
#include <nifti1.h>

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
		int16_t *pixels;
		int numPixels = w*h*ns;
		
		if (numPixels < 1) return;
		
		if (NT < ns) {
			for (int i=NT;i<ns;i++) {
				glGenTextures(1, &texture[i]);
				int err = glGetError();
				
				if (err!=GL_NO_ERROR) printf("Error in glGenTextures %i -> %i!\n",i,err);
			}
			NT = ns;
		}
			
		pixels = (int16_t *) pixelData;
				
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

sap_essentials_t essProtInfo;


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
		
	if (header_def.data_type != DATATYPE_INT16) {
		fprintf(stderr, "Data type != int16\n");
		return false;
	}
	
	printf("\nHeader information: %i samples / %i channels\n\n", header_def.nsamples, header_def.nchans);
	
	
	const ft_chunk_t *chunk = find_chunk(protBuffer.data(), 0, protBuffer.size(), FT_CHUNK_NIFTI1);
	if (chunk != NULL && chunk->def.size == sizeof(nifti_1_header)) {
		nifti_1_header *NH = (nifti_1_header *) chunk->data;
		if (!strcmp(NH->magic, "ni1") || !strcmp(NH->magic,"n+1")) {
			printf("Got NIFTI-1 header!\n");
			essProtInfo.readoutPixels  = NH->dim[1];
			essProtInfo.phasePixels    = NH->dim[2];
			essProtInfo.numberOfSlices = NH->dim[3];
			printf("Resolution (px)...: %i x %i x %i\n", essProtInfo.readoutPixels, essProtInfo.phasePixels, essProtInfo.numberOfSlices);
			printf("Voxel size (mm) ..: %f x %f x %f\n", NH->pixdim[1], NH->pixdim[2], NH->pixdim[2]);
			return true;
		}
	}
	chunk = find_chunk(protBuffer.data(), 0, protBuffer.size(), FT_CHUNK_SIEMENS_AP);
	if (chunk != NULL) {
		sap_item_t *PI = sap_parse(chunk->data, chunk->def.size);
		printf("Got Siemens ASCII protocol information!\n");
		if (sap_get_essentials(PI, &essProtInfo) != SAP_NUM_ESSENTIALS) {
			printf("Not all information could be parsed :-(\n");
		}
		printf("Resolution (px)...: %i x %i x %i\n", essProtInfo.readoutPixels, essProtInfo.phasePixels, essProtInfo.numberOfSlices);
		printf("FOV (mm)..........: %f x %f\n", essProtInfo.readoutFOV, essProtInfo.phaseFOV);
		printf("Slice thickness...: %f\n", essProtInfo.sliceThickness);
		printf("TR (microsec.)....: %li\n", essProtInfo.TR);
		printf("#Contrasts........: %i\n", essProtInfo.numberOfContrasts);
		sap_destroy(PI);
		return true;
	}
	printf("No meta information (e.g. resolution) found.\n");
	return false;
}


// this will get called repeatedly from the GUI loop, use this to poll for new data
void idleCall(void *dummy) {
	SimpleStorage pixBuffer;
	datadef_t data_def;
	FtBufferRequest request;
	FtBufferResponse response;
	unsigned int newSamples, newEvents;
	
	if (essProtInfo.numberOfSlices == 0) {
		if (!readHeader()) return;
	}
	
	request.prepWaitData(prevSamples, 0xFFFFFFFF, 50);
	
	if (tcprequest(ftbSocket, request.out(), response.in()) < 0) {
		fprintf(stderr, "Error in communication. Buffer server aborted??\n");
		return;
	}
	if (!response.checkWait(newSamples, newEvents)) {
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
	
	px2tex.getFromBuffer(pixBuffer.data(), essProtInfo.readoutPixels, essProtInfo.phasePixels, essProtInfo.numberOfSlices);
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
	
	// this will trigger reading the header during the idleCall
	essProtInfo.numberOfSlices = 0;
	
	printf("Trying to connect to fieldtrip buffer at %s:%d\n",hostname,port);
	
	ftbSocket = open_connection(hostname, port);
	if (ftbSocket<=0) {
		fprintf(stderr,"Failed\n");
		exit(1);
	}
		
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
	
	// "Run" the FLTK message loop once: This is needed to create 
	// an OpenGL context before trying to generate textures etc.
	Fl::wait(0);		
	
	Fl::add_idle(idleCall);
	
	BW->redraw();
	
	Fl::run();
	//delete window;
	delete BW;
	
	printf("Closing connection...\n");
	close_connection(ftbSocket);
}
