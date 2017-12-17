# Copyright (C) 2010 Stefan Klanke
# Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
# Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
# Kapittelweg 29, 6525 EN Nijmegen, The Netherlands

# Use MinGW for compilation on Windows

# To recompile, please download the Emotiv EDK, and update the EDKPATH variable below.
# An alternative would be to reimplement using https://github.com/Emotiv/community-sdk

ifneq "$(OS)" "Windows_NT"
	OS      ?= $(shell uname -s)
	MACHINE ?= $(shell uname -m)
endif

PLATFORM ?= $(shell gcc -dumpmachine)

FIELDTRIP = ../../../..
FTBUFFER = $(FIELDTRIP)/realtime/src/buffer
EDKPATH = C:\Program\ Files\ (x86)\Emotiv\ Education\ Edition\ SDK_v1.0.0.5-PREMIUM\doc\examples

# defaults, might be overwritten further down
CC       = gcc
CXX      = g++
CXXFLAGS = -O2 -w # -Wunused -Wall -pedantic
INCPATH  = -I$(FTBUFFER)/src -I$(FTBUFFER)/cpp -I$(FIELDTRIP)/realtime/src/external/emotiv
LIBPATH  = -L$(FTBUFFER)/src
LDLIBS   = -lbuffer -lpthread

# better would be to make a proper library of the C++ framework
CPPOBJS = $(patsubst %,$(FTBUFFER)/cpp/%, SignalConfiguration.o GdfWriter.o FtConnection.o StringServer.o)

ifeq "$(PLATFORM)" "i686-pc-cygwin"
	BINDIR   = $(FIELDTRIP)/realtime/bin/win32
	LDLIBS  += -lm -lws2_32 -L$(BINDIR) -static
	SUFFIX   = .exe
	fixpath  = $1
	RM       = rm -f
endif

ifeq "$(PLATFORM)" "x86_64-pc-cygwin"
	BINDIR   = $(FIELDTRIP)/realtime/bin/win64
	LDLIBS  += -lm -lws2_32 -L$(BINDIR) -static
	SUFFIX   = .exe
	fixpath  = $1
	RM       = rm -f
endif

ifeq "$(PLATFORM)" "i686-w64-mingw32"
	BINDIR   = $(FIELDTRIP)/realtime/bin/win32
	# prevent dependency on libwinpthread-1.dll by linking statically
	# see http://stackoverflow.com/questions/13768515/how-to-do-static-linking-of-libwinpthread-1-dll-in-mingw
	LDLIBS  += -lm -lws2_32 -L$(BINDIR) -static
	SUFFIX   = .exe
	# make clean target work on windows
	fixpath  = $(subst /,\,$1)
	RM       = del
endif

ifeq "$(PLATFORM)" "x86_64-w64-mingw32"
	BINDIR   = $(FIELDTRIP)/realtime/bin/win64
	# prevent dependency on libwinpthread-1.dll by linking statically
	# see http://stackoverflow.com/questions/13768515/how-to-do-static-linking-of-libwinpthread-1-dll-in-mingw
	LDLIBS  += -lm -lws2_32 -L$(BINDIR) -static
	SUFFIX   = .exe
	# make clean target work on windows
	fixpath  = $(subst /,\,$1)
	RM       = del
endif

ifndef BINDIR
$(error Unsupported platform: $(PLATFORM) :/.)
endif

TARGETS = $(BINDIR)/emotiv2ft$(SUFFIX)

##############################################################################
all: $(TARGETS)

%.o: %.cc
	$(CXX) $(INCPATH) $(CXXFLAGS) -c $<

$(BINDIR)/%$(SUFFIX): %.o $(CPPOBJS)
	$(CXX) -o $@ $^ $(LDFLAGS) $(LIBPATH) $(LDLIBS)

clean:
	$(RM) *.o $(call fixpath, $(TARGETS))
