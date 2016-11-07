# Copyright (C) 2015 Robert Oostenveld
# Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
# Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
# Kapittelweg 29, 6525 EN Nijmegen, The Netherlands

# Use MinGW for compilation on Windows!

PLATFORM = $(shell gcc -dumpmachine)
UNAME = $(shell uname)
MACHINE = $(shell uname -m)

FIELDTRIP = ../../../..
FTBUFFER = $(FIELDTRIP)/realtime/src/buffer
CPPDIR = $(FTBUFFER)/cpp

# better would be to make a proper library of the C++ framework
CPPOBJ = 

# defaults, might be overwritten further down
CFLAGS = -I. -I$(FTBUFFER)/src -I$(FIELDTRIP)/realtime/src/external/inih/ -I$(CPPDIR)
# CFLAGS += -Wunused -Wall -pedantic -O3 -g 
CFLAGS += -w
LDFLAGS = -L$(FTBUFFER)/src 
LDLIBS = -lbuffer $(CPPOBJ)

ifeq "$(PLATFORM)" "mingw32"
	BINDIR = $(FIELDTRIP)/realtime/bin/win32
	SUFFIX = .exe
	CFLAGS += -I$(FIELDTRIP)/realtime/src/external/pthreads-win32/include
	# LDFLAGS += -L$(BINDIR) does not seem to work for DLLs. Hence:
	LDFLAGS += -L$(FIELDTRIP)/realtime/src/external/pthreads-win32/lib
	LDLIBS += -lws2_32 -lwinmm $(BINDIR)/pthreadGC2.dll
	# make clean target work on windows:
	fixpath = $(subst /,\,$1)
endif

ifeq "$(UNAME)" "Linux"
	LDLIBS += -ldl -lstdc++ -lpthread -lrt
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
		LDFLAGS += -m32 
	endif
	ifeq "$(MACHINE)" "x86_64"
		BINDIR = $(FIELDTRIP)/realtime/bin/maci64
		CFLAGS += -m64
		LDFLAGS += -m64
	endif
	ifeq "$(MACHINE)" "Power Macintosh"
		BINDIR = $(FIELDTRIP)/realtime/bin/mac
	endif
endif

ifndef BINDIR
$(error Unsupported platform: $(PLATFORM) :/.)
endif

TARGETS = $(patsubst %,$(BINDIR)/%$(SUFFIX), openbci2ft)

##############################################################################
all: $(TARGETS)

$(TARGETS):

$(BINDIR)/%$(SUFFIX): ini.o serial.o %.o
	$(CC) -o $@ $^ $(LDFLAGS) $(LDLIBS)

ini.o: $(FIELDTRIP)/realtime/src/external/inih/ini.c
	$(CC) $(CFLAGS) -c -o $@ $^ 

clean:
	$(RM) *.o $(call fixpath, $(TARGETS))
