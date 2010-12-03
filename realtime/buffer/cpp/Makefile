# Copyright (C) 2010 Stefan Klanke
# Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
# Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
# Kapittelweg 29, 6525 EN Nijmegen, The Netherlands

# Use MinGW for compilation on Windows!

PLATFORM = $(shell gcc -dumpmachine)

#defaults, might be over-written further down
LDFLAGS = -lpthread ../src/socketserver.o ../src/libbuffer.a -ldl 
INCPATH = -I../src -I. 
SUFFIX = 

ifeq "$(PLATFORM)" "mingw32"
   LDFLAGS = pthreadGC2.dll ../src/socketserver.o ../src/libbuffer.a -LC:/msys/1.0/local/lib -lws2_32 -lwinmm 
   INCPATH = -I../src -I. -I../pthreads-win32/include
   SUFFIX = .exe
endif

CFLAGS = $(INCPATH) -Wunused -Wall -pedantic -O3 -g $(FLAGS)
LDFLAGS += $(FLAGS)

all: odmTest$(SUFFIX)

FtConnection.o: FtConnection.cc FtBuffer.h
	g++ -c FtConnection.cc $(CFLAGS)	
	
SignalConfiguration.o: SignalConfiguration.cc SignalConfiguration.h
	g++ -c SignalConfiguration.cc $(CFLAGS)	
	
StringServer.o: StringServer.cc StringServer.h
	g++ -c StringServer.cc $(CFLAGS)	
	
GdfWriter.o: GdfWriter.cc GdfWriter.h 
	g++ -c  GdfWriter.cc $(CFLAGS)	
	
odmTest.o: odmTest.cc OnlineDataManager.h StringServer.o
	g++ -c odmTest.cc $(CFLAGS)

odmTest$(SUFFIX): odmTest.o SignalConfiguration.o GdfWriter.o FtConnection.o StringServer.o
	g++ -o $@ $^ $(LDFLAGS) 