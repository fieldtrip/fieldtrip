# Copyright (C) 2010 Stefan Klanke
# Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
# Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
# Kapittelweg 29, 6525 EN Nijmegen, The Netherlands

# Use MinGW for compilation on Windows!

# Put the EMOTIV SDK header files into the subdirectory "include", or modify the FLAGS below
# The MinGW compiler will directly link against edk.dll, so you'll need to copy that file,
# as well as edk_utils.dll to the same directory as emotiv2ft.[cc/exe]

PLATFORM = $(shell gcc -dumpmachine)

ifeq "$(PLATFORM)" "mingw32"
CPPDIR  = ../../buffer/cpp
LDFLAGS = pthreadGC2.dll ../../buffer/src/libbuffer.a -lws2_32 -lwinmm edk.dll
INCPATH = -I../../buffer/src -I$(CPPDIR) -I../../buffer/pthreads-win32/include -I. -Iinclude
SUFFIX = .exe
endif

CFLAGS = $(INCPATH) -Wunused -Wall -pedantic -O3 -g

CPPOBJS = $(CPPDIR)/SignalConfiguration.o $(CPPDIR)/GdfWriter.o $(CPPDIR)/FtConnection.o $(CPPDIR)/StringServer.o

all: emotiv2ft$(SUFFIX)
	
emotiv2ft.o: emotiv2ft.cc
	g++ -c emotiv2ft.cc $(CFLAGS)
	
emotiv2ft$(SUFFIX): emotiv2ft.o ../../buffer/src/socketserver.o $(CPPOBJS)
	g++ -o $@ $^ $(LDFLAGS)

