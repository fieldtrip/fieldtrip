# Copyright (C) 2010 Stefan Klanke
# Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
# Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
# Kapittelweg 29, 6525 EN Nijmegen, The Netherlands

# Use MinGW for compilation on Windows

ifneq "$(OS)" "Windows_NT"
	OS = $(shell uname -s)
	MACHINE = $(shell uname -m)
endif

PLATFORM = $(shell gcc -dumpmachine)

FIELDTRIP = ../../../..
FTBUFFER  = $(FIELDTRIP)/realtime/src/buffer

# defaults, might be overwritten further down
CC       = gcc
CXX      = g++
CFLAGS   = -O2 -w # -Wunused -Wall -pedantic
INCPATH  = -I$(FTBUFFER)/src -I.
LIBPATH  = -L$(FTBUFFER)/src
LDLIBS   = -lpthread -lbuffer

ifeq "$(PLATFORM)" "mingw32"
	BINDIR    = $(FIELDTRIP)/realtime/bin/win32
	CFLAGS   += -I$(FIELDTRIP)/realtime/src/external/pthreads-win32/include
	# LDFLAGS += -L$(BINDIR) does not seem to work for DLLs. Hence:
	LDLIBS   += -lws2_32 -lwinmm $(BINDIR)/pthreadGC2.dll
	SUFFIX    = .exe
	# make clean target work on windows:
	fixpath   = $(subst /,\,$1)
	RM        = del
endif

ifeq "$(PLATFORM)" "i686-w64-mingw32"
	# prevent dependency on libwinpthread-1.dll by linking statically
	# see http://stackoverflow.com/questions/13768515/how-to-do-static-linking-of-libwinpthread-1-dll-in-mingw
	LDLIBS   = -lbuffer -lm -lpthread -lws2_32 -static
	BINDIR   = $(FIELDTRIP)/realtime/bin/win32
	SUFFIX   = .exe
	# make clean target work on windows
	fixpath  = $(subst /,\,$1)
	RM       = del
endif

ifeq "$(PLATFORM)" "x86_64-w64-mingw32"
	# prevent dependency on libwinpthread-1.dll by linking statically
	# see http://stackoverflow.com/questions/13768515/how-to-do-static-linking-of-libwinpthread-1-dll-in-mingw
	LDLIBS   = -lbuffer -lm -lpthread -lws2_32 -static
	BINDIR   = $(FIELDTRIP)/realtime/bin/win64
	SUFFIX   = .exe
	# make clean target work on windows
	fixpath  = $(subst /,\,$1)
	RM       = del
endif

ifeq "$(OS)" "Linux"
	LDLIBS += -lpthread -ldl -lrt
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

ifeq "$(OS)" "Darwin"
	LDLIBS += -lpthread
	fixpath = $1
	ifeq "$(MACHINE)" "i386"
		BINDIR = $(FIELDTRIP)/realtime/bin/maci
		CFLAGS   += -m32
		CXXFLAGS += -m32
		LDFLAGS  += -m32
	endif
	ifeq "$(MACHINE)" "x86_64"
		BINDIR = $(FIELDTRIP)/realtime/bin/maci64
		CFLAGS   += -m64
		CXXFLAGS += -m64
		LDFLAGS  += -m64
	endif
	ifeq "$(MACHINE)" "Power Macintosh"
		BINDIR = $(FIELDTRIP)/realtime/bin/mac
	endif
endif

ifndef BINDIR
$(error Unsupported platform: $(PLATFORM) :/.)
endif

TARGETS = $(BINDIR)/buffer$(SUFFIX) $(BINDIR)/buffer$(SUFFIX) $(BINDIR)/buffer_rda$(SUFFIX)

###############################################################################
all: $(TARGETS)

%.o: %.c
	$(CC) $(INCPATH) $(CFLAGS) -c $<

$(BINDIR)/%$(SUFFIX): %.o
	$(CC) -o $@ $^ $(LDFLAGS) $(LIBPATH) $(LDLIBS)

clean:
	$(RM) *.o $(call fixpath, $(TARGETS))
