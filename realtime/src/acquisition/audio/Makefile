# Copyright (C) 2010 Bart Niehuis, Stefan Klanke
# Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
# Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
# Kapittelweg 29, 6525 EN Nijmegen, The Netherlands

# Use MinGW for compilation on Windows

# Make sure you have PortAudio installed.

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
CXXFLAGS = -Wunused -Wall -pedantic -O3
INCPATH  = -I$(FTBUFFER)/src -I$(FTBUFFER)/cpp
LIBPATH  = -L$(FTBUFFER)/src
LDLIBS   = -lbuffer

# better would be to make a proper library of the C++ framework
CPPOBJS = $(patsubst %,$(FTBUFFER)/cpp/%, SignalConfiguration.o GdfWriter.o FtConnection.o StringServer.o)

ifeq "$(PLATFORM)" "mingw32"
	BINDIR   = $(FIELDTRIP)/realtime/bin/win32
	INCPATH += -I../../external/portaudio
	# LIBPATH += -L$(BINDIR) does not seem to work for DLLs. Hence:
	LDLIBS  += -lws2_32 $(BINDIR)/pthreadGC2.dll $(BINDIR)/PortAudio.dll
	SUFFIX   = .exe
	RM       = del
	# make clean target work on windows:
	fixpath  = $(subst /,\,$1)
endif

ifeq "$(PLATFORM)" "i686-w64-mingw32"
	# prevent dependency on libwinpthread-1.dll by linking statically
	# see http://stackoverflow.com/questions/13768515/how-to-do-static-linking-of-libwinpthread-1-dll-in-mingw
	INCPATH += -I../../external/portaudio
	LDLIBS   = -lbuffer -lm -lpthread -lws2_32 -static
	BINDIR   = $(FIELDTRIP)/realtime/bin/win32
	SUFFIX   = .exe
	RM       = del /Q
	# make clean target work on windows
	fixpath  = $(subst /,\,$1)
endif

ifeq "$(PLATFORM)" "x86_64-w64-mingw32"
	# prevent dependency on libwinpthread-1.dll by linking statically
	# see http://stackoverflow.com/questions/13768515/how-to-do-static-linking-of-libwinpthread-1-dll-in-mingw
	INCPATH += -I../../external/portaudio
	LDLIBS   = -lbuffer -lm -lpthread -lws2_32 -static
	BINDIR   = $(FIELDTRIP)/realtime/bin/win64
	SUFFIX   = .exe
	RM       = del /Q
	# make clean target work on windows
	fixpath  = $(subst /,\,$1)
endif

ifeq "$(OS)" "Linux"
	fixpath = $1
	LDLIBS += -ldl -lpthread -lportaudio
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
	fixpath = $1
	LDLIBS +=-lpthread -lportaudio
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

TARGETS = $(BINDIR)/audio2ft$(SUFFIX)

##############################################################################
all: $(TARGETS)

%.o: %.cc
	$(CXX) $(INCPATH) $(CXXFLAGS) -c $<

$(BINDIR)/%$(SUFFIX): %.o $(CPPOBJS)
	$(CXX) -o $@ $^ $(LIBPATH) $(LDLIBS)

clean:
	$(RM) *.o $(call fixpath, $(TARGETS))
