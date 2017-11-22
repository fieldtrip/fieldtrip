/*
 * Graphical viewer for displaying online data that is streamed to the buffer.
 *
 * Copyright (C) 2010, Stefan Klanke
 * Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
 * Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 */

#include <FL/Fl.H>
#include <FL/Fl_Double_Window.H>
//#include <FL/Fl_Menu.H>
#include <FL/Fl_Multi_Browser.H>
#include <FL/Fl_Slider.H>
#include <FL/Fl_Scrollbar.H>
#include <FL/Fl_Check_Button.H>
#include <FL/Fl_Input.H>
#include <FL/Fl_Box.H>
#include <FL/fl_draw.H>
#include <FL/fl_ask.H>

#include <stdlib.h>
#include <math.h>
//#include <sys/time.h>
#include <stdio.h>
#include <string.h>
#include <FtBuffer.h>
#include <MultiChannelFilter.h>

#define HPFILTORD 2
#define HPFREQ  4.0

#define LPFILTORD 7
#define LPFREQ 70.0

#define HIDDEN  16  // color code for hiding particular channels

class OnlineDataDisplay;
class ColorBrowser;

char **labels;
int *colorTable;

int numChannels = 0;
unsigned int numSamples = 0;
OnlineDataDisplay *oddis = 0;
Fl_Slider *slideSpace, *slideScale;
Fl_Scrollbar *scrollbar;
ColorBrowser *browser;
SimpleStorage rawStore,floatStore;
Fl_Check_Button *hpButton, *lpButton;
Fl_Button *conButton;
Fl_Input *addrField;
Fl_Box *scaleBox;
bool useHighpass;
bool useLowpass;
MultiChannelFilter<float,float> *hpFilter = NULL;
MultiChannelFilter<float,float> *lpFilter = NULL;
FtConnection ftCon;

bool readHeader();

/** This is the widget on which the channels are drawn.
  It maintains its own little ringbuffer and includes a
  scrollbar widget.
 */
class OnlineDataDisplay : public Fl_Widget {
  public:

    OnlineDataDisplay(int x, int y, int w, int h, const char *label = NULL);
    virtual ~OnlineDataDisplay();

    void draw();
    // int handle(int event);
    void setSize(int nChans, int nSamples, int skip);
    void addSamples(int N, const float *data);

    // catch window resize events to recalibrate the scroll bar etc.
    void resize(int x, int y, int w, int h) {
      Fl_Widget::resize(x,y,w,h);
      if (height == h) return;
      height = h;
      do_callback();
    }

    int getHeight() const {
      return height;
    }

    // set scrollbar parameter
    void setScrollSize(int sPos, float scale, int space) {
      yPos = sPos;
      yScale = scale;
      ySpace = space;
      redraw();
    }

  protected:
    float *data;
    int nChans;
    int nSamp;
    int numTotal, pos;
    int yPos, ySpace;
    float yScale;
    int height;
    int skipped;
};

OnlineDataDisplay::OnlineDataDisplay(int x,int y,int w,int h,const char *label) : Fl_Widget(x,y,w,h,label) {
  data = NULL;
  pos = 0;
  yPos = 0;
  yScale = 1;
  ySpace = 10;
  height = h-4;
  nChans = 0;
}

void OnlineDataDisplay::setSize(int nChans, int nSamples, int skip) {
  delete[] data;

  this->nSamp = nSamples;
  this->nChans = nChans;

  data = new float[nChans * nSamp];
  numTotal = 0;
  pos = 0;
  skipped = skip;
}

OnlineDataDisplay::~OnlineDataDisplay() {
  delete[] data;
}

void OnlineDataDisplay::draw() {
  fl_draw_box(FL_DOWN_BOX, x(), y(), w(), h(), FL_WHITE);

  fl_push_clip(x()+2,y()+2,w()-4,h()-4);

  int ix = 0;
  int numDraw = (nSamp < numTotal) ? nSamp : numTotal;

  float xScale = (float) (w()-44)/nSamp;

  int oldFont = fl_font();
  int oldSize = fl_size();
  fl_font(FL_HELVETICA, 10);

  int ypos = y() + 2 - yPos;

  fl_color(FL_GRAY);
  for (int yp = ypos; yp < y()+h(); yp+=ySpace) {
    if (yp<0) continue;
    fl_line(x()+2, yp, x()+w()-4, yp);
  }

  for (int i=0;i<nChans;i++) {
    if (colorTable[i]==HIDDEN) continue;

    ix++;
    ypos += ySpace;

    if (ypos < -50) continue;
    if (ypos > height+50) continue;

    fl_color(colorTable[i]);
    fl_draw(labels[i], x()+2, ypos);

    fl_push_matrix();
    fl_translate(x()+42, ypos);
    fl_scale(xScale, -yScale*ySpace);

    fl_begin_line();
    for (int j=0;j<numDraw;j++) {
      fl_vertex(j,data[i+j*nChans]);
    }
    fl_end_line();
    fl_pop_matrix();
  }

  int xpos = x()+42 + (int) (xScale*pos);
  fl_color(FL_GRAY);
  fl_rectf(xpos,y()+2,4,h()-4);

  fl_pop_clip();
  fl_font(oldFont, oldSize);
}

/*
   int OnlineDataDisplay::handle(int event) {

   if (event == FL_PUSH) {
   int mouse_x = Fl::event_x();
   int mouse_y = Fl::event_y();

   if (Fl::event_state() & FL_BUTTON3) {
   Fl_Menu_Item m[] = {{"Aaa"},{"Bbb"},0};

   m[0].popup(mouse_x, mouse_y, 0, m, 0);

   redraw();
   }

   return 1;
   }
   return 0;
   }
 */

void OnlineDataDisplay::addSamples(int N, const float *sdata) {
  /*
     if (N>nSamp) {
     data += nChans * (N-nSamp);
     N = nSamp;
     }
   */
  for (int j=0;j<N;j++) {
    float *d = data + ((pos + j) % nSamp)*nChans;
    const float *s = sdata + j*nChans;
    for (int i=0;i<nChans;i++) d[i] = s[i];
  }

  pos=(pos + N) % nSamp;
  numTotal+=N;
  redraw();
}

class ColorBrowser : public Fl_Multi_Browser {
  public:
    ColorBrowser(int x,int y,int w, int h, const char *label=NULL) : Fl_Multi_Browser(x,y,w,h,label) {
      sizeTable = 0;
      colorTable = NULL;
      labels = NULL;
      numActive= 0;
      //color2(fl_rgb_color(200,200,255));
    }

    virtual ~ColorBrowser() {}

    void setupTable(int N, const char * const *labels, int *colorTable) {
      char tmp[80];
      clear();

      this->colorTable = colorTable;
      this->labels = labels;
      sizeTable = N;
      for (int n=0;n<N;n++) {
	setColor(n,0, tmp);
	add(tmp);
      }
      numActive = N;
    }

    int handle(int event) {
      if (event != FL_KEYBOARD) return Fl_Multi_Browser::handle(event);

      int code = -1;
      switch(Fl::event_key()) {
	case ' ': deselect(); return 1;
	case 'r': code = 1; break; // red
	case 'g': code = 2; break; // green
	case 'b': code = 4; break; // blue
	case 'k': code = 0; break; // black
	case 'n': code = HIDDEN; break; // hide=none
	case 'h': code = HIDDEN; break; // hide=none
	case 'y': code = 94; break; // yellow=3
	case 'l': code = 220; break; // lime
	case 'p': code = 89; break; // purple
	case 'a': code = 8; break; // gray
	case 'o': code = 92; break; // orange
	case 'c': code = 6; break; // cyan
	case 'm': code = 5; break; // magenta
	case 't': code = 181; break; // turquoise
	default: return 0;
      }
      int count = 0;
      for (int n=0;n<size();n++) {
	if (selected(1+n)) setColor(n, code);
	if (colorTable[n] != HIDDEN) count++;
      }
      numActive = count;
      redraw();
      do_callback();
      return 1;
    }

    int getNumActive() {
      return numActive;
    }

  protected:

    void setColor(int n, int code) {
      char tmp[80];
      setColor(n,code,tmp);
    }

    void setColor(int n, int code, char *tmp) {
      colorTable[n] = code;
      snprintf(tmp, 80, "@C%i@.%s", code, labels[n]);
      text(1+n, tmp);
    }

    int sizeTable, numActive;
    int *colorTable;
    const char * const *labels;
};

void hpCallback(Fl_Widget*, void*) {
  useHighpass = hpButton->value();
  if (useHighpass && hpFilter) hpFilter->clear();
}

void lpCallback(Fl_Widget*, void*) {
  useLowpass = lpButton->value();
  if (useLowpass && lpFilter) lpFilter->clear();
}

void generalCallback(Fl_Widget*, void*)  {
  int   ySpace = (int) slideSpace->value();
  float sv = slideScale->value();
  float yScale = powf(10.0, sv);

  int   sPos = scrollbar->value();
  int   oddisHeight = oddis->getHeight();
  int   numSelected = browser->getNumActive();

  scrollbar->value(sPos, oddisHeight, 0, (numSelected+2)*ySpace);

  char scLab[12];


  int ex = (int) floorf(sv);
  float ma = powf(10.0f, sv-ex);

  if (ex > 0) {
    snprintf(scLab, 11, "%.2ge+%i", ma, ex);
  } else if (ex < 0) {
    snprintf(scLab, 11, "%.2ge%i", ma, ex);
  } else {
    snprintf(scLab, 11, "%.3g", ma);
  }
  scaleBox->copy_label(scLab);

  oddis->setScrollSize(sPos, yScale, ySpace);
}

void connectCallback(Fl_Widget*, void*) {
  if (ftCon.isOpen()) {
    ftCon.disconnect();
    addrField->textcolor(FL_BLACK);
    addrField->redraw();
    conButton->label("Connect");
  } else {
    if (ftCon.connect(addrField->value())) {
      addrField->textcolor(FL_DARK_GREEN);
      addrField->redraw();
      conButton->label("Disconnect");
      readHeader();
    }
  }
}

void freeTables() {
  if (colorTable != NULL) free(colorTable);
  if (labels != NULL) {
    for (int i=0;i<numChannels;i++) {
      if (labels[i] != NULL) free(labels[i]);
    }
    free(labels);
  }
  numChannels = 0;
  labels = NULL;
}



void errDisconnect(const char *msg) {
  ftCon.disconnect();
  numChannels = 0; // reset to enforce reading header next time
  addrField->textcolor(FL_BLACK);
  addrField->redraw();
  conButton->label("Connect");
  fl_alert(msg);
}

bool readHeader() {
  SimpleStorage chunkBuffer;
  headerdef_t header_def;
  FtBufferRequest request;
  FtBufferResponse response;

  request.prepGetHeader();

  if (tcprequest(ftCon.getSocket(), request.out(), response.in()) < 0) {
    errDisconnect("Error in communication - check buffer server");
    return false;
  }

  if (!response.checkGetHeader(header_def, &chunkBuffer)) {
    fprintf(stderr, "Could not read header.\n");
    return false;
  }

  freeTables();
  numChannels = header_def.nchans;
  numSamples = header_def.nsamples;

  labels = (char **) calloc(numChannels, sizeof(char *));
  colorTable = (int *) calloc(numChannels, sizeof(int));

  const ft_chunk_t *cnc = find_chunk(chunkBuffer.data(), 0, chunkBuffer.size(), FT_CHUNK_CHANNEL_NAMES);
  if (cnc == NULL) {
    printf("No channel names found\n");
    for (int n=0;n<numChannels;n++) {
      labels[n] = (char *) malloc(8);
      snprintf(labels[n],7,"#%i",n+1);
    }
  } else {
    const char *s = (const char *) cnc->data;
    for (int n=0;n<numChannels;n++) {
      int ln = strlen(s);
      if (ln==0) {
	labels[n] = (char *) malloc(8);
	snprintf(labels[n],7,"#%i",n+1);
      } else {
	labels[n] = strdup(s);
      }
      s+=ln+1;
    }
  }
  browser->setupTable(numChannels, labels, colorTable);
  oddis->setSize(numChannels, (int) (4.0*header_def.fsample), numSamples);

  if (hpFilter != NULL) {
    delete hpFilter;
    hpFilter = NULL;
  }
  if (lpFilter != NULL) {
    delete lpFilter;
    lpFilter = NULL;
  }
  hpFilter = new MultiChannelFilter<float,float>(numChannels, HPFILTORD);
  hpFilter->setButterHP(HPFREQ/header_def.fsample);
  lpFilter = new MultiChannelFilter<float,float>(numChannels, LPFILTORD);
  lpFilter->setButterLP(LPFREQ/header_def.fsample);
  generalCallback(NULL, 0);

  return true;
}



template<typename T>
void convertToFloat(float *dest, const void *src, unsigned int nsamp, unsigned int nchans) {
  const T *srcT = static_cast<const T *>(src);
  for (unsigned int j=0;j<nsamp;j++) {
    for (unsigned int i=0;i<nchans;i++) {
      dest[i] = (float) srcT[i];
    }
    dest += nchans;
    srcT += nchans;
  }
}

// this will get called repeatedly from the GUI loop, use this to poll for new data
void idleCall(void *dummy) {
  datadef_t ddef;
  FtBufferRequest request;
  FtBufferResponse response;
  unsigned int newSamples, newEvents;

  if (!ftCon.isOpen()) return;

  if (numChannels == 0) {
    if (!readHeader()) {
#ifdef WIN32
      Sleep(50);
#else
      usleep(50000);
#endif
      return;
    }
  }

  // wait for new samples, don't care about events, up to 40ms
  request.prepWaitData(numSamples, 0xFFFFFFFF, 40);

  if (tcprequest(ftCon.getSocket(), request.out(), response.in()) < 0) {
    errDisconnect("Error in communication. Buffer server aborted??");
    return;
  }
  if (!response.checkWait(newSamples, newEvents)) {
    errDisconnect("Error in received packet - disconnecting...");
    return;
  }

  if (newSamples == numSamples) return; // nothing new
  if (newSamples < numSamples) {
    // oops ? do we have a new header?
    if (!readHeader()) return;
    if (numSamples == 0) return; // no data yet
    if (numSamples > 1024 || numChannels > 512) {
      // "lots" of data already in the buffer
      // -> don't do anything with that data
      //    continue next idleCall
      return;
    }
    // read data from the start of the buffer up to newSamples right away
    newSamples = numSamples;
    numSamples = 0;
  }

  request.prepGetData(numSamples, newSamples-1);

  if (tcprequest(ftCon.getSocket(), request.out(), response.in()) < 0) {
    errDisconnect("Error in communication. Buffer server aborted??");
    return;
  }
  if (!response.checkGetData(ddef, &rawStore)) {
    errDisconnect("Error in received packet - disconnecting...");
    return;
  }

  floatStore.resize(sizeof(float) * ddef.nsamples * ddef.nchans);

  float *fdata = (float *) floatStore.data();
  switch(ddef.data_type) {
    case DATATYPE_UINT8:
      convertToFloat<uint8_t>(fdata, rawStore.data(), ddef.nsamples, ddef.nchans);
      break;
    case DATATYPE_INT8:
      convertToFloat<int8_t>(fdata, rawStore.data(), ddef.nsamples, ddef.nchans);
      break;
    case DATATYPE_UINT16:
      convertToFloat<uint16_t>(fdata, rawStore.data(), ddef.nsamples, ddef.nchans);
      break;
    case DATATYPE_INT16:
      convertToFloat<int16_t>(fdata, rawStore.data(), ddef.nsamples, ddef.nchans);
      break;
    case DATATYPE_UINT32:
      convertToFloat<uint32_t>(fdata, rawStore.data(), ddef.nsamples, ddef.nchans);
      break;
    case DATATYPE_INT32:
      convertToFloat<int32_t>(fdata, rawStore.data(), ddef.nsamples, ddef.nchans);
      break;
    case DATATYPE_UINT64:
      convertToFloat<uint64_t>(fdata, rawStore.data(), ddef.nsamples, ddef.nchans);
      break;
    case DATATYPE_INT64:
      convertToFloat<int64_t>(fdata, rawStore.data(), ddef.nsamples, ddef.nchans);
      break;
  	case DATATYPE_FLOAT32:
	  convertToFloat<float>(fdata, rawStore.data(), ddef.nsamples, ddef.nchans);
      break;
    case DATATYPE_FLOAT64:
      convertToFloat<double>(fdata, rawStore.data(), ddef.nsamples, ddef.nchans);
      break;
  }
  if (useHighpass) {
    hpFilter->process(ddef.nsamples, fdata, fdata); // in place
  }
  if (useLowpass) {
    lpFilter->process(ddef.nsamples, fdata, fdata); // in place
  }
  oddis->addSamples(ddef.nsamples, fdata);

  numSamples = newSamples;
}


int main(int argc, char **argv) {
  Fl_Double_Window *window = new Fl_Double_Window(10,10,800,430,"FieldTrip Realtime Buffer Viewer");

  addrField = new Fl_Input(100,10,210,20,"Buffer address");
  addrField->textsize(12);
  addrField->textfont(FL_COURIER);
  addrField->labelsize(12);

  conButton = new Fl_Button(320,10,100,20,"Connect");
  conButton->labelsize(12);
  conButton->callback(connectCallback);

  hpButton = new Fl_Check_Button(440,10,100,20,"Highpass filter");
  hpButton->labelsize(12);
  hpButton->callback(hpCallback);

  lpButton = new Fl_Check_Button(540,10,100,20,"Lowpass filter");
  lpButton->labelsize(12);
  lpButton->callback(lpCallback);

  oddis = new OnlineDataDisplay(10,40,620,380);
  oddis->callback(generalCallback);

  browser = new ColorBrowser(700,40,90,380);
  browser->textsize(12);

  scrollbar = new Fl_Scrollbar(630,40,16,380);
  scrollbar->value(0,380,0,16*40);
  scrollbar->callback(generalCallback);

  slideScale = new Fl_Slider(660,40,28,150);
  slideScale->labelsize(12);
  slideScale->label("scale");
  slideScale->range(6, -6); // log scale
  slideScale->step(0.04);
  slideScale->value(0);
  slideScale->type(FL_VERT_NICE_SLIDER);
  slideScale->callback(generalCallback);

  scaleBox = new Fl_Box(650,220,46,20,"1.0");
  scaleBox->box(FL_THIN_DOWN_BOX);
  scaleBox->labelsize(10);

  slideSpace = new Fl_Slider(660,250,28,150);
  slideSpace->labelsize(12);
  slideSpace->label("space");
  slideSpace->range(100, 10);
  slideSpace->value(40);
  slideSpace->type(FL_VERT_NICE_SLIDER);
  slideSpace->callback(generalCallback);

  window->resizable(oddis);
  window->end();
  window->show();

  generalCallback(NULL, 0);
  Fl::add_idle(idleCall);

  if (argc > 1) {
    addrField->value(argv[1]);
    conButton->do_callback();
  } else {
    addrField->value("localhost:1972");
  }

  Fl::run();
  delete window;
}
