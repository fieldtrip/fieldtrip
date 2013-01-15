# Copyright (C) 2010 Stefan Klanke
# Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
# Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
# Kapittelweg 29, 6525 EN Nijmegen, The Netherlands

# Use MinGW for compilation on Windows!

PLATFORM = $(shell gcc -dumpmachine)
UNAME   = $(shell uname)
MACHINE = $(shell uname -m)

FIELDTRIP = ../../../..
FTBUFFER = $(FIELDTRIP)/realtime/src/buffer

# defaults, might be overwritten further down
CFLAGS = $(INCPATH) -Wunused -Wall -pedantic -g
LDFLAGS = -lpthread ../../buffer/src/libbuffer.a
INCPATH = -I../../buffer/src -I. 
SUFFIX = 

ifeq "$(PLATFORM)" "mingw32"
   LDFLAGS = pthreadGC2.dll ../../buffer/src/libbuffer.a -LC:/msys/1.0/local/lib -lws2_32 -lwinmm 
   INCPATH = -I../../buffer/src -I. -I../../buffer/pthreads-win32/include
   SUFFIX = .exe
endif

ifeq "$(UNAME)" "Linux"
    ifeq "$(MACHINE)" "i686"
        BINDIR = $(FIELDTRIP)/realtime/bin/glnxa64
    endif
    ifeq "$(MACHINE)" "x86_64"
        BINDIR = $(FIELDTRIP)/realtime/bin/glnx86
    endif
    ifeq "$(MACHINE)" "armv6l"
        BINDIR = $(FIELDTRIP)/realtime/bin/raspberrypi
    endif
endif

ifeq "$(UNAME)" "Darwin"
    ifeq "$(MACHINE)" "i386"
        BINDIR = $(FIELDTRIP)/realtime/bin/maci
    endif
    ifeq "$(MACHINE)" "x86_64"
        BINDIR = $(FIELDTRIP)/realtime/bin/maci64
    endif
endif

ifndef BINDIR
$(warning Unsupported platform $(PLATFORM) :/.)
endif

TARGETS = $(BINDIR)/playback$(SUFFIX) $(BINDIR)/recording$(SUFFIX)

###############################################################################
all: $(TARGETS)

playback.o: playback.c
	gcc -c playback.c $(CFLAGS)
	
recording.o: recording.c
	gcc -c $^ $(CFLAGS)
	
ft_offline.o: ft_offline.c
	gcc -c $^ $(CFLAGS)
			
$(BINDIR)/playback$(SUFFIX): playback.o
	gcc -o $@ $^ $(LDFLAGS)
		
$(BINDIR)/recording$(SUFFIX): recording.o ft_offline.o ../../buffer/src/socketserver.o
	gcc -o $@ $^ $(LDFLAGS)	

clean:
	$(RM) *.o $(call fixpath, $(TARGETS))
