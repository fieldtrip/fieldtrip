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
CFLAGS = $(INCPATH) -Wunused -Wall -pedantic -g
LDFLAGS = -lpthread ../../buffer/src/libbuffer.a
LDLIBS =
INCPATH = -I../../buffer/src -I. 
SUFFIX = 

ifeq "$(PLATFORM)" "mingw32"
	LDFLAGS = pthreadGC2.dll ../../buffer/src/libbuffer.a -LC:/msys/1.0/local/lib -lws2_32 -lwinmm 
	LDLIBS =
	INCPATH = -I../../buffer/src -I. -I../../buffer/pthreads-win32/include
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

TARGETS = $(BINDIR)/playback$(SUFFIX) $(BINDIR)/recording$(SUFFIX)

###############################################################################
all: $(TARGETS)

playback.o: playback.c
	$(CC) -c playback.c $(CFLAGS)
	
recording.o: recording.c
	$(CC) -c $^ $(CFLAGS)
	
ft_storage.o: ft_storage.c
	$(CC) -c $^ $(CFLAGS)
			
$(BINDIR)/playback$(SUFFIX): playback.o
	$(CC) -o $@ $^ $(LDFLAGS) $(LDLIBS)
		
$(BINDIR)/recording$(SUFFIX): recording.o ft_storage.o ../../buffer/src/socketserver.o
	$(CC) -o $@ $^ $(LDFLAGS) $(LDLIBS)

clean:
	$(RM) *.o $(call fixpath, $(TARGETS))
