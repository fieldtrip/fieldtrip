# Copyright (C) 2010 Stefan Klanke
# Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
# Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
# Kapittelweg 29, 6525 EN Nijmegen, The Netherlands

# Use MinGW for compilation on Windows!

PLATFORM = $(shell gcc -dumpmachine)

ifeq "$(PLATFORM)" "mingw32"
   LDFLAGS = pthreadGC2.dll ../../buffer/src/libbuffer.a -LC:/msys/1.0/local/lib -lws2_32 -lwinmm 
   INCPATH = -I../../buffer/src -I../../buffer/cpp -I. -I../../buffer/pthreads-win32/include
   SUFFIX = .exe
else
   LDFLAGS = -lpthread ../../buffer/src/libbuffer.a
   INCPATH = -I../../buffer/src -I../../buffer/cpp -I. 
   SUFFIX = 
endif

CFLAGS = $(INCPATH) -Wunused -Wall -pedantic -O3 -g

all: biosemi2ft$(SUFFIX)
	
biosemi2ft.o: biosemi2ft.cc
	g++ -c biosemi2ft.cc $(CFLAGS)

BioSemiClient.o: BioSemiClient.cc BioSemiClient.h
	g++ -c BioSemiClient.cc $(CFLAGS)
	
biosemi2ft$(SUFFIX): biosemi2ft.o BioSemiClient.o
	g++ -o $@ $^ $(LDFLAGS) 
