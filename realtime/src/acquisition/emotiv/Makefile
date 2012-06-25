# Copyright (C) 2010 Stefan Klanke
# Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
# Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
# Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
#
# Use MinGW for compilation on Windows!
#
# To recompile, please download the Emotiv EDK, and update the EDKPATH variable
# below.

PLATFORM = $(shell gcc -dumpmachine)
FTPATH = ../../../..
FTBUFFER = $(FTPATH)/realtime/src/buffer

EDKPATH = C:\Program\ Files\ (x86)\Emotiv\ Education\ Edition\ SDK_v1.0.0.5-PREMIUM\doc\examples

CXXFLAGS = -I$(FTBUFFER)/src -I$(FTBUFFER)/cpp -Wunused -Wall
LDFLAGS = -L$(FTBUFFER)/src 
LDLIBS = -lbuffer

# TODO: make proper library of C++ framework.
CPPOBJ = $(patsubst %,$(FTBUFFER)/cpp/%,\
	SignalConfiguration.o GdfWriter.o FtConnection.o StringServer.o)

ifeq "$(PLATFORM)" "mingw32"
	BINDIR = $(FTPATH)/realtime/bin/win32
	
	CXXFLAGS += -I$(EDKPATH)/include
	LDFLAGS += -L$(EDKPATH)/lib
	LDLIBS += $(CPPOBJ) -lws2_32 -lwinmm $(BINDIR)/edk.dll $(BINDIR)/pthreadGC2.dll 
	EXTENSION = .exe
	
	# Make clean target work on windows:
	fixpath = $(subst /,\,$1)
	RM = del
	CP = copy /Y
endif

ifndef BINDIR
$(warning Unsupported plaform: $(PLATFORM) :/.)
endif

TARGETS = $(BINDIR)/emotiv2ft$(EXTENSION)

all: $(TARGETS)

$(BINDIR)/%$(EXTENSION): %.cc
	$(CXX) -o $@ $^ $(CPPOBJS) $(CXXFLAGS) $(LDFLAGS) $(LDLIBS)

clean:
	$(RM) *.o $(call fixpath, $(TARGETS))