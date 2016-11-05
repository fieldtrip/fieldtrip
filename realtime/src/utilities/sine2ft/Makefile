# Copyright (C) 2010 Stefan Klanke
# Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
# Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
# Kapittelweg 29, 6525 EN Nijmegen, The Netherlands

# Use MinGW for compilation on Windows!
# This Makefile assumes that FLTK has been installed to a place where the compiler can find it.
# On Windows, this location is C:/msys/1.0/local -- change this if need be

PLATFORM = $(shell gcc -dumpmachine)
UNAME = $(shell uname)
MACHINE = $(shell uname -m)

FIELDTRIP = ../../../..
FTBUFFER = $(FIELDTRIP)/realtime/src/buffer

# defaults, might be overwritten further down
LDFLAGS = `fltk-config --ldflags` -lpthread ../../buffer/src/libbuffer.a
CXXFLAGS = `fltk-config --cxxflags` -I../../buffer/src -I../../buffer/cpp -I. -Wunused -Wall -pedantic
SUFFIX = 

ifeq "$(PLATFORM)" "mingw32"
	LDFLAGS = pthreadGC2.dll ../../buffer/src/libbuffer.a -LC:/msys/1.0/local/lib -lws2_32 -lwinmm -mwindows -mno-cygwin -lfltk -lole32 -luuid -lcomctl32 -lgdi32 -lwsock32
	CXXFLAGS = -IC:/msys/1.0/local/include -mwindows
	SUFFIX = .exe
endif

ifeq "$(UNAME)" "Linux"
	LDLIBS += -lrt
	ifeq "$(MACHINE)" "i686"
		BINDIR = $(FIELDTRIP)/realtime/bin/glnx86
	endif
	ifeq "$(MACHINE)" "x86_64"
		BINDIR = $(FIELDTRIP)/realtime/bin/glnxa64
	endif
	ifeq "$(MACHINE)" "armv6l"
		BINDIR = $(FIELDTRIP)/realtime/bin/raspberrypi
	endif
	ifeq "$(MACHINE)" "armv7l"
		BINDIR = $(FIELDTRIP)/realtime/bin/raspberrypi
	endif
endif

ifeq "$(UNAME)" "Darwin"
	ifeq "$(MACHINE)" "i386"
		BINDIR = $(FIELDTRIP)/realtime/bin/maci
		CFLAGS += -m32
		CXXFLAGS += -m32
		LDFLAGS += -m32
	endif
	ifeq "$(MACHINE)" "x86_64"
		BINDIR = $(FIELDTRIP)/realtime/bin/maci64
		CFLAGS += -m64
		CXXFLAGS += -m64
		LDFLAGS += -m64
	endif
	ifeq "$(MACHINE)" "Power Macintosh"
		BINDIR = $(FIELDTRIP)/realtime/bin/mac
	endif
endif

ifndef BINDIR
$(error Unsupported platform: $(PLATFORM) :/.)
endif

TARGETS = $(BINDIR)/sine2ft$(SUFFIX)

###############################################################################
all: $(TARGETS)

sine2ft.o: sine2ft.cc
	$(CXX) -c sine2ft.cc $(CXXFLAGS)

$(BINDIR)/sine2ft$(SUFFIX): sine2ft.o
	$(CXX) -o $@ sine2ft.o $(LDFLAGS) $(LDLIBS)

clean:
	$(RM) *.o $(call fixpath, $(TARGETS))
