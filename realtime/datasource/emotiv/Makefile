# Copyright (C) 2010 Stefan Klanke
# Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
# Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
# Kapittelweg 29, 6525 EN Nijmegen, The Netherlands

# Use MinGW for compilation on Windows!

# Put the EMOTIV SDK into subdirectories include, or modify the FLAGS below

PLATFORM = $(shell gcc -dumpmachine)

ifeq "$(PLATFORM)" "mingw32"
   LDFLAGS = pthreadGC2.dll ../../buffer/src/libbuffer.a -lws2_32 -lwinmm edk.dll
   INCPATH = -I../../buffer/src -I../../buffer/cpp -I. -I../../buffer/pthreads-win32/include -Iinclude
   SUFFIX = .exe
else
   LDFLAGS = -lpthread ../../buffer/src/libbuffer.a
   INCPATH = -I../buffer/src -I../../buffer/cpp -I. 
   SUFFIX = 
endif

CFLAGS = $(INCPATH) -Wunused -Wall -pedantic -O3 -g

all: emotiv2ft$(SUFFIX)
	
emotiv2ft.o: emotiv2ft.cc
	g++ -c emotiv2ft.cc $(CFLAGS)
	
emotiv2ft$(SUFFIX): emotiv2ft.o ../../buffer/src/socketserver.o
	g++ -o $@ $^ $(LDFLAGS) 
