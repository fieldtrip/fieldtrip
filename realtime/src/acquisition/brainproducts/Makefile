# Copyright (C) 2010 Stefan Klanke
# Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
# Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
# Kapittelweg 29, 6525 EN Nijmegen, The Netherlands

# Use MinGW for compilation on Windows!

PLATFORM = $(shell gcc -dumpmachine)
UNAME = $(shell uname)
MACHINE = $(shell uname -m)

FIELDTRIP = ../../../..
FTBUFFER = $(FIELDTRIP)/realtime/src/buffer

# defaults, might be overwritten further down
CFLAGS = -Wunused -Wall -pedantic
CFLAGS += -I$(FTBUFFER)/src -I. 
LDFLAGS = -L$(FTBUFFER)/src
LDLIBS = -lbuffer -lpthread

ifeq "$(PLATFORM)" "mingw32"
	BINDIR = $(FIELDTRIP)/realtime/bin/win32
	SUFFIX = .exe
	CFLAGS += -I$(FIELDTRIP)/realtime/src/external/pthreads-win32/include
	# LDFLAGS += -L$(BINDIR) does not seem to work for DLLs. Hence:
	LDLIBS += -lws2_32 -lwinmm $(BINDIR)/pthreadGC2.dll
	# make clean target work on windows:
	fixpath = $(subst /,\,$1)
	RM = del
endif

ifeq "$(UNAME)" "Linux"
	LDLIBS += -lrt
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
	ifeq "$(MACHINE)" "Power Macintosh"
		BINDIR = $(FIELDTRIP)/realtime/bin/mac
	endif
endif

ifndef BINDIR
$(error Unsupported plaform: $(PLATFORM) :/.)
endif

TARGETS = $(BINDIR)/rda2ft$(SUFFIX)

all: $(TARGETS)
	
$(BINDIR)/rda2ft$(SUFFIX): rda2ft.c
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(LDLIBS)

clean:
	$(RM) $(call fixpath, *.o $(TARGETS))
