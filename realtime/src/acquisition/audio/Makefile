# Copyright (C) 2010 Bart Niehuis, Stefan Klanke
# Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
# Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
# Kapittelweg 29, 6525 EN Nijmegen, The Netherlands

# Use MinGW for compilation on Windows!
# On other platforms, make sure you have PortAudio installed.

PLATFORM = $(shell gcc -dumpmachine)
UNAME = $(shell uname)
MACHINE = $(shell uname -m)

FIELDTRIP = ../../../..
FTBUFFER = $(FIELDTRIP)/realtime/src/buffer
CPPDIR	= $(FTBUFFER)/cpp

# defaults, might be overwritten further down
INCPATH = -I$(FTBUFFER)/src -I$(CPPDIR) -I. 
LDFLAGS = -L$(FTBUFFER)/src
CFLAGS = $(INCPATH) -Wunused -Wall -pedantic -O3
CXXFLAGS = $(CFLAGS)

# better would be to make a proper library of the C++ framework
CPPOBJS = $(patsubst %,$(CPPDIR)/%,\
	SignalConfiguration.o GdfWriter.o FtConnection.o StringServer.o)

ifeq "$(PLATFORM)" "mingw32"
	BINDIR = $(FIELDTRIP)/realtime/bin/win32
	SUFFIX = .exe
	fixpath = $(subst /,\,$1)
	# LDFLAGS += -L$(BINDIR) does not seem to work for DLLs. Hence:
	LDLIBS += $(CPPOBJS) -lbuffer -lws2_32 $(BINDIR)/pthreadGC2.dll $(BINDIR)/PortAudio.dll 
	INCPATH += ../../external/portaudio
	RM = del
endif

ifeq "$(UNAME)" "Linux"
	fixpath = $1
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
	LDFLAGS += -L$(FTBUFFER)/src
 	LDLIBS += $(CPPOBJS) -lbuffer -lpthread -ldl -lportaudio
endif

ifeq "$(UNAME)" "Darwin"
	fixpath = $1
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
	LDLIBS += $(CPPOBJS) -lpthread -lportaudio
endif

ifndef BINDIR
$(error Unsupported platform: $(PLATFORM) :/.)
endif

TARGETS = $(BINDIR)/audio2ft$(SUFFIX)

##############################################################################
all: $(TARGETS)

$(BINDIR)/%$(SUFFIX): %.cc
	$(CXX) -o $@ $(CFLAGS) $(LDFLAGS) $	 $(LDLIBS)

clean:
	$(RM) *.o $(call fixpath, $(TARGETS))

