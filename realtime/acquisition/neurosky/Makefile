# Copyright (C) 2010 Stefan Klanke
# Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
# Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
# Kapittelweg 29, 6525 EN Nijmegen, The Netherlands

# Use MinGW for compilation on Windows!

PLATFORM = $(shell gcc -dumpmachine)
CPPDIR   = ../../buffer/cpp

ifeq "$(PLATFORM)" "mingw32"
   LDFLAGS = ../../buffer/pthreads-win32/lib/pthreadGC2.dll ../../buffer/src/libbuffer.a -LC:/msys/1.0/local/lib -lws2_32 -lwinmm 
   INCPATH = -I../../buffer/src -I. -I$(CPPDIR) -I../../buffer/pthreads-win32/include
   SUFFIX = .exe
else 
   LDFLAGS = -lpthread ../../buffer/src/libbuffer.a -ldl 
   INCPATH = -I../../buffer/src -I. -I$(CPPDIR)
   SUFFIX = 
endif

CFLAGS = $(INCPATH) -Wunused -Wall -pedantic -O3 -g

all: thinkgear2ft$(SUFFIX) 

serial.o: serial.c
	gcc -c serial.c -I.
	
thinkgear2ft.o: thinkgear2ft.cc
	g++ -c thinkgear2ft.cc $(CFLAGS)
	
thinkgear2ft$(SUFFIX): thinkgear2ft.o serial.o
	g++ -o $@ $^ $(CFLAGS) $(CPPDIR)/SignalConfiguration.cc $(CPPDIR)/GdfWriter.o $(CPPDIR)/FtConnection.o $(CPPDIR)/StringServer.o $(LDFLAGS)
