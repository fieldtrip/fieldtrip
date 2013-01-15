# Copyright (C) 2010 Stefan Klanke
# Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
# Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
# Kapittelweg 29, 6525 EN Nijmegen, The Netherlands

# Use MinGW for compilation on Windows!

PLATFORM = $(shell gcc -dumpmachine)
FIELDTRIP = ../../../..
FTBUFFER = $(FIELDTRIP)/realtime/src/buffer
CPPDIR = $(FTBUFFER)/cpp

# better would be to make a proper library of the C++ framework
CPPOBJ = $(patsubst %,$(CPPDIR)/%,\
	SignalConfiguration.o GdfWriter.o FtConnection.o StringServer.o)

# defaults, might be overwritten further down
CFLAGS = -Wunused -Wall -pedantic -O3 -g -I.
CXXFLAGS = $(CFLAGS) -I$(FTBUFFER)/src -I$(CPPDIR)
LDFLAGS = -L$(FTBUFFER)/src
LDLIBS = -lbuffer $(CPPOBJ)

ifeq "$(PLATFORM)" "mingw32"
	BINDIR = $(FIELDTRIP)/realtime/bin/win32
	EXTENSION = .exe
	
	CXXFLAGS += -I$(FIELDTRIP)/realtime/src/external/pthreads-win32/include
	# LDFLAGS += -L$(BINDIR) does not seem to work for DLLs. Hence:
	LDLIBS += -lws2_32 -lwinmm $(BINDIR)/pthreadGC2.dll
	
	# Make clean target work on windows:
	fixpath = $(subst /,\,$1)
	RM = del
endif

ifeq "$(PLATFORM)" "i686-apple-darwin10"
	BINDIR = $(FIELDTRIP)/realtime/bin/maci64
	CFLAGS += -m64
	CXXFLAGS += -m64
	LDFLAGS += -m64
	fixpath = $1
endif

ifeq "$(PLATFORM)" "x86_64-linux-gnu"
	LDLIBS += -lpthread -ldl 
	BINDIR = $(FIELDTRIP)/realtime/bin/glnxa64
endif

ifeq "$(PLATFORM)" "x86_64-redhat-linux"
	LDLIBS += -lpthread -ldl 
	BINDIR = $(FIELDTRIP)/realtime/bin/glnxa64
endif

ifeq "$(PLATFORM)" "arm-linux-gnueabi"
	LDLIBS += -lpthread -ldl 
	BINDIR = $(FIELDTRIP)/realtime/bin/raspberrypi
endif

ifndef BINDIR
$(error Unsupported platform: $(PLATFORM) :/.)
endif

TARGETS = $(BINDIR)/thinkgear2ft$(EXTENSION)

all: $(TARGETS)

thinkgear2ft$(EXTENSION): thinkgear2ft.cc serial.o

$(BINDIR)/%$(EXTENSION): %.o
	$(CXX) -o $@ $^ $(LDFLAGS) $(LDLIBS)

$(TARGETS): serial.o

clean:
	$(RM) *.o $(call fixpath, $(TARGETS))
