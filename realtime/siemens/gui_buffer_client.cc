#include <Fl/Fl.h>
#include <Fl/Fl_Window.h>
#include <Fl/fl_draw.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <buffer.h>
#include <siemensap.h>

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
		uint16_t *pixels;
		int numPixels = w*h*ns;
		
		printf("To tonvert: %i\n", numPixels);
		if (numPixels < 1) return;
		
		pixels = (uint16_t *) pixelData;
				
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
		printf("Converted: %i\n", numPixels);
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



void idleCall(void *dummy) {
	const char *protocol;
	message_t request;
	messagedef_t request_def;
	message_t *response = NULL;
	headerdef_t header_def;
	datasel_t datasel;
		
	request.def = &request_def;
	request.buf = NULL;
	request_def.version = VERSION;
	request_def.command = GET_HDR;
	request_def.bufsize = 0;
	
	if (tcprequest(ftbSocket, &request, &response) < 0) {
		fprintf(stderr, "Error in communication. Buffer server aborted??\n");
		return;
	}
	
	if (!response) {
		fprintf(stderr, "GET_HDR: unknown error in response\n");
		return;
	}
	if (!response->def) {
		fprintf(stderr, "GET_HDR: unknown error in response\n");
		goto cleanup;
	}
	if (response->def->command!=GET_OK) {
		fprintf(stderr, "GET_HDR: Buffer returned an error (%d)\n", response->def->command);
		Sleep(50);
		goto cleanup;
	}
	
	memcpy(&header_def, response->buf, sizeof(header_def));
	
	protocol = (const char *)response->buf + sizeof(header_def);
	getProtInfo(protocol, header_def.bufsize);
	
	free(response->buf);
	free(response->def);
	free(response);
		
	if (header_def.data_type != DATATYPE_UINT16) {
		fprintf(stderr, "Data type != uint16\n");
		return;
	}
	
	printf("GET_HDR: samples / channels: %i / %i\n", header_def.nsamples, header_def.nchans);
	if (header_def.nsamples == prevSamples) {
		Sleep(50);
		return;
	} 
	prevSamples = header_def.nsamples;

    datasel.begsample = prevSamples-1;
    datasel.endsample = prevSamples-1;
	
	request_def.version = VERSION;
	request_def.command = GET_DAT;
	request_def.bufsize = sizeof(datasel_t);
	request.buf = &datasel;
	
	if (tcprequest(ftbSocket, &request, &response) < 0) {
		fprintf(stderr, "Error in communication. Buffer server aborted??\n");
		return;
	}
	
	if (!response) {
		fprintf(stderr, "GET_DAT: unknown error in response\n");
		return;
	}
	if (!response->def) {
		fprintf(stderr, "GET_DAT: unknown error in response\n");
		goto cleanup;
	}
	if (response->def->command!=GET_OK) {
		fprintf(stderr, "GET_DAT: Buffer returned an error (%d)\n", response->def->command);
		goto cleanup;
	}
	
	void *data_buf = (void *)((char *)response->buf + sizeof(datadef_t));
	
	px2i.getFromBuffer(data_buf, readoutResolution, phaseResolution, numSlices);
	IW->redraw();
	
cleanup:
	if (response) {
		if (response->buf) free(response->buf);
		if (response->def) free(response->def);
		free(response);
	}
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
	window = new Fl_Window(800,100,454,454,"FMRI buffer client");
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
