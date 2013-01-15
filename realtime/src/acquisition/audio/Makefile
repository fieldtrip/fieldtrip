# Copyright (C) 2010 Bart Niehuis, Stefan Klanke
# Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
# Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
# Kapittelweg 29, 6525 EN Nijmegen, The Netherlands

# Use MinGW for compilation on Windows!
# On other platforms, make sure you have PortAudio installed.

PLATFORM = $(shell gcc -dumpmachine)
FIELDTRIP = ../../../..
FTBUFFER = $(FIELDTRIP)/realtime/src/buffer
CPPDIR  = $(FTBUFFER)/cpp
PORTAUDIO = ../../external/portaudio

# defaults, might be overwritten further down
INCPATH = -I$(PORTAUDIO) -I$(FTBUFFER)/src -I$(CPPDIR) -I. 
LDFLAGS = -L$(FTBUFFER)/src
CFLAGS = $(INCPATH) -Wunused -Wall -pedantic -O3
CXXFLAGS = $(CFLAGS)

# better would be to make a proper library of the C++ framework
CPPOBJ = $(patsubst %,$(CPPDIR)/%,\
	SignalConfiguration.o GdfWriter.o FtConnection.o StringServer.o)

ifeq "$(PLATFORM)" "mingw32"
	BINDIR = $(FIELDTRIP)/realtime/bin/win32
	EXTENSION = .exe
	fixpath = $(subst /,\,$1)
	# LDFLAGS += -L$(BINDIR) does not seem to work for DLLs. Hence:
	LDLIBS += $(CPPOBJ) -lbuffer -lws2_32 $(BINDIR)/pthreadGC2.dll $(BINDIR)/PortAudio.dll 
	RM = del
endif

ifndef BINDIR
$(error Unsupported platform: $(PLATFORM) :/.)
endif

TARGETS = $(BINDIR)/audio2ft$(EXTENSION)

##############################################################################
all: $(TARGETS)

$(BINDIR)/%$(EXTENSION): %.cc
	$(CXX) -o $@ $(CFLAGS) $(LDFLAGS) $^ $(LDLIBS)

clean:
	$(RM) *.o $(call fixpath, $(TARGETS))
