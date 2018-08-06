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
#include <FtBuffer.h>
#include <siemensap.h>
#include <nifti1.h>

#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <FL/fl_draw.H>

class PixelData2Image {
	public:
	
	PixelData2Image() {
		numAlloc = 0;
		image = 0;
	}
	
	~PixelData2Image() {
		if (image) free(image);
	}
	
	void getFromBuffer(void *pixelData, int w, int h, int ns) {
		int16_t *pixels;
		int numPixels = w*h*ns;
		
		if (numPixels < 1) return;
		
		pixels = (int16_t *) pixelData;
				
		if (numPixels > numAlloc) {
			if (image != 0) free(image);
			image = (unsigned char *) malloc(numPixels);
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
			image[i] = v32 / 1500; /* use proper scaling sometime */
		}
	}
			
	const unsigned char *getImage() const {
		return image;
	}
		
	int getW() const {
		return W;
	}
	
	int getH() const {
		return H;
	}
	
	int getNS() const {
		return NS;
	}
	
	protected:
	
	unsigned char *image;
	int numAlloc;
	int W, H, NS;
};

class ImageWidget : public Fl_Widget {
	public:
	
	ImageWidget(int x, int y, int w, int h, const PixelData2Image *px, const char *label = NULL) : Fl_Widget(x,y,w,h,label) {
		this->px = px;
	}
	
	virtual void draw() {
		if (px) {
			const unsigned char *image = px->getImage();
			if (image) {
				int i = 0, j = 0;
				int w = px->getW();
				int h = px->getH();
				int ns = px->getNS();
								
				for (int n=0;n<ns;n++) {
					fl_draw_image_mono(image + n*w*h, 10+j*(w+10), 10+i*(h+10), w, h);
					if (++j==6) {
						i++;
						j=0;
					}
				}
			}
		}
	}
	
	protected:
	
	const PixelData2Image *px;
};


Fl_Window *window;
ImageWidget *IW;
PixelData2Image px2i;
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
	
	px2i.getFromBuffer(pixBuffer.data(), essProtInfo.readoutPixels, essProtInfo.phasePixels, essProtInfo.numberOfSlices);
	IW->redraw();
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
	
	Fl::visual(FL_RGB);
	window = new Fl_Window(200,100,454,454,"FMRI buffer client");
	IW = new ImageWidget(0,0,454,454, &px2i);
	window->resizable(IW);
	window->end();
	
	window->show();
	
	Fl::add_idle(idleCall);
	
	Fl::run();
	delete window;
	
	printf("Closing connection...\n");
	close_connection(ftbSocket);
}
