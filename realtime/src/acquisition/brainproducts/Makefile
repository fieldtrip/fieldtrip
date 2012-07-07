# Copyright (C) 2010 Stefan Klanke
# Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
# Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
# Kapittelweg 29, 6525 EN Nijmegen, The Netherlands

# Use MinGW for compilation on Windows!

PLATFORM = $(shell gcc -dumpmachine)
FIELDTRIP = ../../../..
FTBUFFER = $(FIELDTRIP)/realtime/src/buffer

#defaults, might be over-written further down
CFLAGS = -Wunused -Wall -pedantic
CFLAGS += -I$(FTBUFFER)/src -I. 
LDFLAGS = -L$(FTBUFFER)/src
LDLIBS = -lbuffer -lpthread

ifeq "$(PLATFORM)" "mingw32"
	fixpath = $(subst /,\,$1)
	BINDIR = $(FIELDTRIP)/realtime/bin/win32
	EXTENSION = .exe

	CFLAGS += -I$(FIELDTRIP)/realtime/src/external/pthreads-win32/include
	# LDFLAGS += -L$(BINDIR) does not seem to work for DLLs. Hence:
	LDLIBS += -lws2_32 -lwinmm $(BINDIR)/pthreadGC2.dll
	RM = del
endif

ifeq "$(PLATFORM)" "i686-apple-darwin10"
	BINDIR = $(FIELDTRIP)/realtime/bin/maci64
	CFLAGS += -m64
	CXXFLAGS += -m64
	LDFLAGS += -m64
	fixpath = $1
endif

ifeq "$(PLATFORM)" "x86_64-redhat-linux"
	BINDIR = $(FIELDTRIP)/realtime/bin/glnxa64
	fixpath = $1
endif

ifeq "$(PLATFORM)" "arm-linux-gnueabi"
	BINDIR = $(FIELDTRIP)/realtime/bin/raspberrypi
	fixpath = $1
endif

ifndef BINDIR
$(error Unsupported plaform: $(PLATFORM) :/.)
endif

TARGETS = $(BINDIR)/rda2ft$(EXTENSION)

all: $(TARGETS)
	
$(BINDIR)/rda2ft$(EXTENSION): rda2ft.c
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(LDLIBS)

clean:
	$(RM) $(call fixpath, *.o $(TARGETS))
