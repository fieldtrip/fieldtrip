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
CPPDIR  = $(FTBUFFER)/cpp

# better would be to make a proper library of the C++ framework
CPPOBJS = $(CPPDIR)/SignalConfiguration.o $(CPPDIR)/GdfWriter.o \
	$(CPPDIR)/FtConnection.o $(CPPDIR)/StringServer.o

# defaults, might be overwritten further down
LDFLAGS = -L$(FTBUFFER)/src
LDLIBS = -lbuffer $(CPPOBJS)
CXXFLAGS = -Wunused -Wall -pedantic -O3 -g
CXXFLAGS += -I$(FTBUFFER)/src -I. -I$(CPPDIR) 

ifeq "$(PLATFORM)" "mingw32"
	BINDIR = $(FIELDTRIP)/realtime/bin/win32
	SUFFIX = .exe
	CXXFLAGS += -I$(FIELDTRIP)/realtime/src/external/pthreads-win32/include
	# LDFLAGS += -L$(BINDIR) does not seem to work for DLLs. Hence:
	LDLIBS += -lws2_32 -lwinmm $(BINDIR)/pthreadGC2.dll
	# make clean target work on windows:
	fixpath = $(subst /,\,$1)
	RM = del
endif

ifeq "$(UNAME)" "Linux"
		LDLIBS += -lbuffer -lpthread -ldl -lrt
		LDFLAGS += -L$(FTBUFFER)/src
    fixpath = $1
    ifeq "$(MACHINE)" "i686"
		    BINDIR = $(FIELDTRIP)/realtime/bin/glnx86
    endif
    ifeq "$(MACHINE)" "x86_64"
		    BINDIR = $(FIELDTRIP)/realtime/bin/glnxa64
    endif
    ifeq "$(MACHINE)" "armv6l"
		    # This is Debian 6.0.5 running on the Raspberry Pi. I have not yet been
		    # able to test whether the binary actually works, that requires having
		    # the amplifier hooked up to the R-Pi.
		    BINDIR = $(FIELDTRIP)/realtime/bin/raspberrypi
    endif
endif

ifeq "$(UNAME)" "Darwin"
    # Currently, only 32 bit OS X is supported since the 64 bit drivers don't seem to
    # work. This means that for correct compilation of biosemi2ft, also libbuffer
    # and ODM need to be (manually) recompiled and linked with '-m32' flag in the
    # CFLAGS and LDFLAGS.
    LDLIBS += -lpthread
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
endif

ifndef BINDIR
$(error Unsupported platform: $(PLATFORM) :/.)
endif

TARGETS = $(BINDIR)/biosemi2ft$(SUFFIX)

##############################################################################
all: $(TARGETS)
	
$(BINDIR)/biosemi2ft$(SUFFIX): biosemi2ft.o BioSemiClient.o 
	$(CXX) -o $@ $^ $(LDFLAGS) $(LDLIBS)

clean:
	$(RM) *.o $(call fixpath, $(TARGETS))
