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
CFLAGS = -I../../buffer/src -I. $(ADDINC) -Wunused -Wall -pedantic
LDFLAGS =
LDLIBS = -lpthread ../../buffer/src/libbuffer.a 
ADDINC = 
SUFFIX = 

ifeq "$(PLATFORM)" "mingw32"
	BINDIR = $(FIELDTRIP)/realtime/bin/win32
	LDFLAGS = -LC:/msys/1.0/local/lib
	LDLIBS = pthreadGC2.dll ../buffer/src/libbuffer.a -lws2_32 -lwinmm 
	ADDINC = -I../../buffer/pthreads-win32/include
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

TARGETS = $(BINDIR)/serial_event$(SUFFIX)

###############################################################################
all: $(TARGETS)

serial.o: serial.c serial.h
	$(CC) $(CFLAGS) -c serial.c -I.
		
serial_event.o: serial_event.c serial.h
	$(CC) $(CFLAGS) -c serial_event.c $(CFLAGS)
	
$(BINDIR)/serial_event$(SUFFIX): serial_event.o serial.o
	$(CC) $(CFLAGS) -o $@ serial_event.o serial.o $(LDFLAGS) $(LDLIBS)

clean:
	$(RM) *.o $(call fixpath, $(TARGETS))
