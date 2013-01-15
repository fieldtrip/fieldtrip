# Copyright (C) 2010, Stefan Klanke
# Copyright (C) 2012, Robert Oostenveld
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
FLFLAGS = `fltk-config --cxxflags`
LDFL = `fltk-config --ldflags`
SUFFIX = 

ifeq "$(PLATFORM)" "mingw32"
   LDFLAGS = pthreadGC2.dll ../../buffer/src/libbuffer.a -LC:/msys/1.0/local/lib -lws2_32 -lwinmm  -mwindows
   FLFLAGS = -IC:/msys/1.0/local/include -mwindows
   LDFL = -mno-cygwin -lfltk -lole32 -luuid -lcomctl32 -lgdi32 -lwsock32
   SUFFIX = .exe
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

TARGETS = $(BINDIR)/event2ft$(SUFFIX)

##############################################################################
all: $(TARGETS)

event2ft.o: event2ft.cc
	g++ -c event2ft.cc $(CXXFLAGS)

$(BINDIR)/event2ft$(SUFFIX): event2ft.o
	g++ -o $@ event2ft.o $(LDFLAGS) $(LDFL)

clean:
	$(RM) *.o $(call fixpath, $(TARGETS))
