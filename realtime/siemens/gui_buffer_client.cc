#include <Fl/Fl.h>
#include <Fl/Fl_Window.h>
#include <Fl/fl_draw.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unixtime.h>
#include <buffer.h>

class PixelData2Image {
	public:
	
	PixelData2Image() {
		numEl = numAlloc = 0;
		image = 0;
	}
	
	~PixelData2Image() {
		if (image) free(image);
	}
	
	/* return value used to be the timestamp */
	double getFromBuffer(void *pixelData, int numPixels) {
		uint16_t *pixels;
		
		if (numPixels < 1) return -1.0;  
		
		pixels = (uint16_t *) pixelData;
		
		W = (int) sqrt(numPixels);
		H = W;
				
		numEl = W*H;
		
		if (numEl > numAlloc) {
			if (image != 0) free(image);
			image = (unsigned char *) malloc(numEl);
			if (image == 0) {
				numAlloc = numEl = 0;
				return -1.0;
			}
			numAlloc = numEl;
		}
		
		for (int i=0;i<numEl;i++) {
			int v32 = pixels[i]*255;
			image[i] = v32 / 1500; /* use proper scaling sometime */
		}
		
		return 0.0;
	}
			
	const unsigned char *getImage() const {
		return image;
	}
	
	int getLength() const {
		return numEl;
	}
	
	int getW() const {
		return W;
	}
	
	int getH() const {
		return H;
	}
	
	protected:
	
	unsigned char *image;
	int numEl, numAlloc;
	int W, H;
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
				fl_draw_image_mono(image, 0, 0, px->getW(), px->getH());
			}
		}
	}
	
	protected:
	
	const PixelData2Image *px;
};


Fl_Window *window;
ImageWidget *IW;
PixelData2Image px2i;
int prevSamples = -1;
int ftbSocket = -1;


void idleCall(void *dummy) {
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
	free(response->buf);
	free(response->def);
	free(response);
		
	if (header_def.data_type != DATATYPE_INT16) {
		fprintf(stderr, "Data type != int16\n");
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
	
	datadef_t *datadef = (datadef_t *) response->buf;
	void *data_buf = (void *)((char *)response->buf + sizeof(datadef_t));
	
	px2i.getFromBuffer(data_buf, datadef->nchans);
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
	window = new Fl_Window(800,100,400,400,"FMRI buffer client");
	IW = new ImageWidget(0,0,400,400, &px2i);
	window->resizable(IW);
	window->end();
	
	window->show();
	
	Fl::add_idle(idleCall);
	
	Fl::run();
	delete window;
	
	printf("Closing connection...\n");
	close_connection(ftbSocket);
}
