# Copyright (C) 2010 Stefan Klanke
# Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
# Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
# Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
#
# Use MinGW for compilation on Windows!
# This Makefile assumes that FLTK has been installed to a place where the compiler can find it.
# On Windows, this location is C:/msys/1.0/local -- change this if need be

PLATFORM = $(shell gcc -dumpmachine)
UNAME   = $(shell uname)
MACHINE = $(shell uname -m)

FIELDTRIP = ../../../..
FTBUFFER = $(FIELDTRIP)/realtime/src/buffer

# defaults, might be overwritten further down
CC = gcc
CXX = g++
CXXFLAGS = $(FLFLAGS) -I../../buffer/src -I../../buffer/cpp -I. -Wunused -Wall -pedantic
LDFLAGS = -lpthread ../../buffer/src/libbuffer.a
FLFLAGS = $(shell fltk-config --cxxflags)
LDFL    = $(shell fltk-config --ldflags)
SUFFIX  = 

# Windows doesn't know "uname", so we need to check for the GCC platform first 
ifeq "$(PLATFORM)" "mingw32"
	FLPREFIX = C:/msys/1.0/local
	LDFLAGS  = pthreadGC2.dll ../../buffer/src/libbuffer.a	-mwindows
	FLFLAGS  = -I$(FLPREFIX)/include -mwindows
	LDFL     = -mno-cygwin -lfltk -lole32 -luuid -lcomctl32 -lgdi32 -lwsock32 -L$(FLPREFIX)/lib
	SUFFIX   = .exe
endif

ifeq "$(UNAME)" "Linux"
    ifeq "$(MACHINE)" "i686"
        BINDIR = $(FIELDTRIP)/realtime/bin/glnxa64
    endif
    ifeq "$(MACHINE)" "x86_64"
        BINDIR = $(FIELDTRIP)/realtime/bin/glnx86
    endif
    ifeq "$(MACHINE)" "armv6l"
        BINDIR = $(FIELDTRIP)/realtime/bin/raspberrypi
    endif
endif

ifeq "$(UNAME)" "Darwin"
    ifeq "$(MACHINE)" "i386"
        BINDIR = $(FIELDTRIP)/realtime/bin/maci
    endif
    ifeq "$(MACHINE)" "x86_64"
        BINDIR = $(FIELDTRIP)/realtime/bin/maci64
    endif
endif

ifndef BINDIR
$(warning Unsupported platform $(PLATFORM) :/.)
endif

TARGETS = $(BINDIR)/bufferViewer$(SUFFIX)

##############################################################################
all: $(TARGETS)

bufferViewer.o: bufferViewer.cc
	g++ -c bufferViewer.cc $(CXXFLAGS)

$(BINDIR)/bufferViewer$(SUFFIX): bufferViewer.o ../../buffer/cpp/FtConnection.o
	g++ -o $@ $^ $(LDFLAGS) $(LDFL)

clean:
	$(RM) *.o $(call fixpath, $(TARGETS))
