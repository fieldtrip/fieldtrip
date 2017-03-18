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
FTBUFFER = $(FIELDTRIP)/realtime/src/buffer

# defaults, might be overwritten further down
CXXFLAGS = -Wunused -Wall -pedantic -O3
INCPATH  = -I$(FTBUFFER)/src -I$(FTBUFFER)/cpp -I.
LIBPATH  = -L$(FTBUFFER)/src
LDLIBS   = -lbuffer -lpthread

# better would be to make a proper library of the C++ framework
CPPOBJS = $(patsubst %,$(FTBUFFER)/cpp/%, SignalConfiguration.o GdfWriter.o FtConnection.o StringServer.o)

ifeq "$(PLATFORM)" "mingw32"
	BINDIR 	 = $(FIELDTRIP)/realtime/bin/win32
	INCPATH += -I$(FIELDTRIP)/realtime/src/external/pthreads-win32/include
	# LIBPATH += -L$(BINDIR) does not seem to work for DLLs. Hence:
	LDLIBS  += -lws2_32 -lwinmm $(BINDIR)/pthreadGC2.dll
	SUFFIX   = .exe
	RM       = del
	# make clean target work on windows:
	fixpath  = $(subst /,\,$1)
endif

ifeq "$(OS)" "Linux"
	LDLIBS += -ldl -lrt
	fixpath = $1
	ifeq "$(MACHINE)" "i686"
		BINDIR = $(FIELDTRIP)/realtime/bin/glnx86
	endif
	ifeq "$(MACHINE)" "x86_64"
		BINDIR = $(FIELDTRIP)/realtime/bin/glnxa64
	endif
	# This is Debian 6.0.5 running on the Raspberry Pi. I have not yet been
	# able to test whether the binary actually works, that requires having
	# the amplifier hooked up to the R-Pi.
	ifeq "$(MACHINE)" "armv6l"
		BINDIR = $(FIELDTRIP)/realtime/bin/raspberrypi
	endif
	ifeq "$(MACHINE)" "armv7l"
		BINDIR = $(FIELDTRIP)/realtime/bin/raspberrypi
	endif
endif

ifeq "$(OS)" "Darwin"
    # Currently, only 32 bit OS X is supported since the 64 bit drivers don't seem to
    # work. This means that for correct compilation of biosemi2ft, also libbuffer
    # and ODM need to be (manually) recompiled and linked with '-m32' flag in the
    # CFLAGS and CXXFLAGS.
    fixpath = $1
		ifeq "$(MACHINE)" "i386"
        BINDIR    = $(FIELDTRIP)/realtime/bin/maci
        CFLAGS   += -m32
        CXXFLAGS += -m32
    endif
		ifeq "$(MACHINE)" "x86_64"
        BINDIR    = $(FIELDTRIP)/realtime/bin/maci64
        CFLAGS   += -m64
        CXXFLAGS += -m64
    endif
endif

ifndef BINDIR
$(error Unsupported platform: $(PLATFORM) :/.)
endif

TARGETS = $(BINDIR)/biosemi2ft$(SUFFIX)

##############################################################################
all: $(TARGETS)

%.o: %.cc
	$(CXX) $(INCPATH) $(CXXFLAGS) -c $<

$(BINDIR)/biosemi2ft$(SUFFIX): biosemi2ft.o BioSemiClient.o
	$(CXX) -o $@ $^ $(LIBPATH) $(LDLIBS) $(CPPOBJS)

clean:
	$(RM) *.o $(call fixpath, $(TARGETS))
