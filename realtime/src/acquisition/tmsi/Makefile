# Copyright (C) 2010 Bart Niehuis, Stefan Klanke
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
CC       = gcc
CXX      = g++
CFLAGS   = -O2 -w # -Wunused -Wall -pedantic
CXXFLAGS = -O2 -w # -Wunused -Wall -pedantic
INCPATH  = -I$(FTBUFFER)/src -I$(FTBUFFER)/cpp -I.
LIBPATH  = -L$(FTBUFFER)/src
LDLIBS   = -lbuffer

# better would be to make a proper library of the C++ framework
CPPOBJS = $(patsubst %,$(FTBUFFER)/cpp/%, SignalConfiguration.o GdfWriter.o FtConnection.o StringServer.o)

ifeq "$(PLATFORM)" "mingw32"
	BINDIR   = $(FIELDTRIP)/realtime/bin/win32
	INCPATH += -I$(FIELDTRIP)/realtime/src/external/pthreads-win32/include
	# LIBPATH += -L$(BINDIR) does not seem to work for DLLs. Hence:
	LDLIBS  += -lws2_32 -lwinmm $(BINDIR)/pthreadGC2.dll
	SUFFIX   = .exe
	RM       = del
	# make clean target work on windows:
	fixpath  = $(subst /,\,$1)
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

ifndef BINDIR
$(error Unsupported platform: $(PLATFORM) :/.)
endif

TARGETS = $(patsubst %,$(BINDIR)/%$(SUFFIX), tmsidriver tmsi2ft)

##############################################################################
all: $(TARGETS)

%.o: %.c
	$(CC) $(INCPATH) $(CFLAGS) -c $<

%.o: %.cpp
	$(CXX) $(INCPATH) $(CXXFLAGS) -c $<

$(BINDIR)/%$(SUFFIX): %.o $(CPPOBJS)
	$(CXX) -o $@ $^ $(LDFLAGS) $(LIBPATH) $(LDLIBS)

$(BINDIR)/tmsidriver$(SUFFIX): Feature.o RTDevice.o

$(BINDIR)/tmsi2ft$(SUFFIX): RTDevice.o

RTDevice.o: RTDevice.cpp RTDevice.h Sadio.h

clean:
	$(RM) *.o $(call fixpath, $(TARGETS))
