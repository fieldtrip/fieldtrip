# Copyright (C) 2010 Bart Niehuis, Stefan Klanke
# Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
# Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
# Kapittelweg 29, 6525 EN Nijmegen, The Netherlands

# Use MinGW for compilation on Windows!

PLATFORM = $(shell gcc -dumpmachine)
UNAME   = $(shell uname)
MACHINE = $(shell uname -m)

FIELDTRIP = ../../../..
FTBUFFER = $(FIELDTRIP)/realtime/src/buffer
CPPDIR  = $(FTBUFFER)/cpp

# better would be to make a proper library of the C++ framework
CPPOBJ = $(patsubst %,$(CPPDIR)/%,\
	SignalConfiguration.o GdfWriter.o FtConnection.o StringServer.o)

# defaults, might be overwritten further down
LDFLAGS = -L$(FTBUFFER)/src
LDLIBS = -lbuffer $(CPPOBJ)
CFLAGS = -I$(FTBUFFER)/src -I. -I$(CPPDIR) -Wunused -Wall -pedantic -O3 -g 
CXXFLAGS = $(CFLAGS)

ifeq "$(PLATFORM)" "mingw32"
	BINDIR = $(FIELDTRIP)/realtime/bin/win32
	SUFFIX = .exe
	CXXFLAGS += -I$(FIELDTRIP)/realtime/src/external/pthreads-win32/include
	# LDFLAGS += -L$(BINDIR) does not seem to work for DLLs. Hence:
	LDLIBS += -lws2_32 -lwinmm $(BINDIR)/pthreadGC2.dll
	# make clean target work on windows:
	fixpath = $(subst /,\,$1)
	RM = del
endif

ifndef BINDIR
$(error Unsupported platform: $(PLATFORM) :/.)
endif

TARGETS = $(patsubst %,$(BINDIR)/%$(SUFFIX), tmsidriver tmsi2ft)

##############################################################################
all: $(TARGETS)

$(BINDIR)/%$(SUFFIX): %.o
	$(CXX) -o $@ $^ $(LDFLAGS) $(LDLIBS)

$(BINDIR)/tmsidriver$(SUFFIX): Feature.o RTDevice.o

$(BINDIR)/tmsi2ft$(SUFFIX): RTDevice.o

RTDevice.o: RTDevice.cpp RTDevice.h Sadio.h

clean:
	$(RM) *.o $(call fixpath, $(TARGETS))
