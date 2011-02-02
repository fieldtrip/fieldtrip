# Copyright (C) 2010 Stefan Klanke
# Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
# Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
# Kapittelweg 29, 6525 EN Nijmegen, The Netherlands

# Use MinGW for compilation on Windows!

PLATFORM = $(shell gcc -dumpmachine)

#defaults, might be over-written further down
CPPDIR  = ../../buffer/cpp
LDFLAGS = -lpthread ../../buffer/src/libbuffer.a -ldl 
INCPATH = -I../../buffer/src -I. -I$(CPPDIR)
SUFFIX = 

ifeq "$(PLATFORM)" "mingw32"
   LDFLAGS = pthreadGC2.dll ../../buffer/src/libbuffer.a -LC:/msys/1.0/local/lib -lws2_32 -lwinmm 
   INCPATH = -I../../buffer/src -I. -I$(CPPDIR) -I../../buffer/pthreads-win32/include
   SUFFIX = .exe
endif

CFLAGS = $(INCPATH) -Wunused -Wall -pedantic -O3 -g $(FLAGS)
LDFLAGS += $(FLAGS)

CPPOBJS = $(CPPDIR)/SignalConfiguration.o $(CPPDIR)/GdfWriter.o $(CPPDIR)/FtConnection.o $(CPPDIR)/StringServer.o

all: biosemi2ft$(SUFFIX) fakebiosemi2ft$(SUFFIX)
	
biosemi2ft.o: biosemi2ft.cc
	g++ -c biosemi2ft.cc $(CFLAGS)

BioSemiClient.o: BioSemiClient.cc BioSemiClient.h
	g++ -c BioSemiClient.cc $(CFLAGS)
	
FakeBioSemiClient.o: FakeBioSemiClient.cc BioSemiClient.h
	g++ -c FakeBioSemiClient.cc $(CFLAGS)

#TODO: For now we compile and pull in the C++ modules during linking
#Think about having a separate library like libbuffer
biosemi2ft$(SUFFIX): biosemi2ft.o BioSemiClient.o 
	g++ -o $@ $^ $(CPPOBJS) $(LDFLAGS)

fakebiosemi2ft$(SUFFIX): biosemi2ft.o FakeBioSemiClient.o 
	g++ -o $@ $^ $(CPPOBJS) $(LDFLAGS)
