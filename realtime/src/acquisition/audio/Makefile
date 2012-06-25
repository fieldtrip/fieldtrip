# Copyright (C) 2010 Bart Niehuis, Stefan Klanke
# Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
# Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
# Kapittelweg 29, 6525 EN Nijmegen, The Netherlands

# Use MinGW for compilation on Windows!
# On other platforms, make sure you have PortAudio installed.

PLATFORM = $(shell gcc -dumpmachine)
CPPDIR  = ../../buffer/cpp

ifeq "$(PLATFORM)" "mingw32"
   LDFLAGS = pthreadGC2.dll portaudio.dll ../../buffer/src/libbuffer.a -lws2_32
   INCPATH = -I. -I../../buffer/src -I$(CPPDIR) -I../../buffer/pthreads-win32/include -Ipa_include
   SUFFIX = .exe
   RM = del
else 
   LDFLAGS = -lpthread -lportaudio ../../buffer/src/libbuffer.a -ldl 
   INCPATH =  -I.  -I../../buffer/src -I$(CPPDIR)
   SUFFIX = 
   RM = rm -f
endif

CFLAGS = $(INCPATH) -Wunused -Wall -pedantic -O3 -g
CPPOBJS = $(CPPDIR)/GdfWriter.o $(CPPDIR)/SignalConfiguration.o $(CPPDIR)/FtConnection.o $(CPPDIR)/StringServer.o

all: audio2ft$(SUFFIX)

audio2ft.o: audio2ft.cc 
	g++ -c $^ $(CFLAGS) 

audio2ft$(SUFFIX): audio2ft.o 
	g++ -o $@ $^ $(CFLAGS) $(CPPOBJS) $(LDFLAGS)

clean:
	$(RM) *.o