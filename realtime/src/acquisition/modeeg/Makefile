# Copyright (C) 2010 Stefan Klanke
# Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
# Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
# Kapittelweg 29, 6525 EN Nijmegen, The Netherlands

# Use MinGW for compilation on Windows!

PLATFORM = $(shell gcc -dumpmachine)
FIELDTRIP = ../../../..
FTBUFFER = $(FIELDTRIP)/realtime/src/buffer
CPPDIR  = $(FTBUFFER)/cpp

# TODO: make proper library of C++ framework.
CPPOBJ = $(patsubst %,$(CPPDIR)/%,\
	SignalConfiguration.o GdfWriter.o FtConnection.o StringServer.o)

#defaults, might be over-written further down
CFLAGS = -I$(FTBUFFER)/src -I. -I$(CPPDIR)
CFLAGS += -Wunused -Wall -pedantic -O3 -g 
CXXFLAGS = $(CFLAGS)
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
endif

ifeq "$(PLATFORM)" "x86_64-linux-gnu"
	LDLIBS += -lpthread -ldl -lstdc++ 
	BINDIR = $(FIELDTRIP)/realtime/bin/glnxa64
endif

ifeq "$(PLATFORM)" "x86_64-redhat-linux"
	LDLIBS += -lpthread -ldl -lstdc++
	BINDIR = $(FIELDTRIP)/realtime/bin/glnxa64
endif

ifndef BINDIR
$(error Unsupported platform: $(PLATFORM) :/.)
endif

TARGETS = $(patsubst %,$(BINDIR)/%$(EXTENSION), \
	modeeg2ft modeeg2ft_2chn modeeg_setupbluesmirf)

all: $(TARGETS)
	
$(BINDIR)/%$(EXTENSION): %.o
	$(CXX) -o $@ $^ $(LDFLAGS) $(LDLIBS)

$(TARGETS): serial.o
$(BINDIR)/modeeg_setupbluesmirf$(EXTENSION): LDLIBS=

clean:
	$(RM) *.o $(call fixpath, $(TARGETS))
