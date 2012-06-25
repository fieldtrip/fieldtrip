# Copyright (C) 2010 Bart Niehuis, Stefan Klanke
# Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
# Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
# Kapittelweg 29, 6525 EN Nijmegen, The Netherlands

# Use MinGW for compilation on Windows!
PLATFORM = $(shell gcc -dumpmachine)
FIELDTRIP = ../../../..
FTBUFFER = $(FIELDTRIP)/realtime/src/buffer
CPPDIR  = $(FTBUFFER)/cpp

#defaults, might be over-written further down
LDFLAGS = -L$(FTBUFFER)/src
CFLAGS = -I$(FTBUFFER)/src -I. -I$(CPPDIR) -Wunused -Wall -pedantic -O3 -g 
CXXFLAGS = $(CFLAGS)

# TODO: make proper library of C++ framework.
CPPOBJ = $(patsubst %,$(CPPDIR)/%,\
	SignalConfiguration.o GdfWriter.o FtConnection.o StringServer.o)

ifeq "$(PLATFORM)" "mingw32"
	BINDIR = $(FIELDTRIP)/realtime/bin/win32
	EXTENSION = .exe
	CXXFLAGS += -I$(FTBUFFER)/pthreads-win32/include
	LDFLAGS += -L$(FTBUFFER)/pthreads-win32
	# LDFLAGS += -L$(BINDIR) does not seem to work for DLLs. Hence:
	LDLIBS += $(CPPOBJ) -lbuffer -lws2_32 -lwinmm $(BINDIR)/pthreadGC2.dll

	# Make clean target work on windows:
	fixpath = $(subst /,\,$1)
	RM = del
endif
ifndef BINDIR
$(warning Unsupported platform: $(PLATFORM) :/.)
endif


TARGETS = $(patsubst %,$(BINDIR)/%$(EXTENSION), tmsidriver tmsi2ft)

all: $(TARGETS)

$(BINDIR)/%$(EXTENSION): %.o
	$(CXX) -o $@ $^ $(LDFLAGS) $(LDLIBS)

$(BINDIR)/tmsidriver$(EXTENSION): Feature.o RTDevice.o
$(BINDIR)/tmsi2ft$(EXTENSION): RTDevice.o
RTDevice.o: RTDevice.cpp RTDevice.h Sadio.h

clean:
	$(RM) *.o $(call fixpath, $(TARGETS))
