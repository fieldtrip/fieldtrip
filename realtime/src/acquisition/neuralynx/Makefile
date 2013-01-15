# Copyright (C) 2010 Stefan Klanke
# Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
# Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
# Kapittelweg 29, 6525 EN Nijmegen, The Netherlands

# Use MinGW for compilation on Windows!

PLATFORM = $(shell gcc -dumpmachine)
FIELDTRIP = ../../../..
FTBUFFER = $(FIELDTRIP)/realtime/src/buffer
CPPDIR = $(FTBUFFER)/cpp

CXXFLAGS = -Wunused -Wall -pedantic -O3
CXXFLAGS += -I$(FTBUFFER)/src -I$(CPPDIR) -I.
LDFLAGS = -L$(FTBUFFER)/src
LDLIBS = -lbuffer

ifeq "$(PLATFORM)" "mingw32"
	BINDIR = $(FIELDTRIP)/realtime/bin/win32
	EXTENSION = .exe
	CXXFLAGS += -I$(FIELDTRIP)/realtime/src/external/pthreads-win32/include
	# LDFLAGS += -L$(BINDIR) does not seem to work for DLLs. Hence:
	LDLIBS += -lws2_32 -lwinmm $(BINDIR)/pthreadGC2.dll $(BINDIR)/MatlabNetComClient.dll
	# make clean target work on windows:
	fixpath = $(subst /,\,$1)
	RM = del
endif

ifndef BINDIR
$(error Unsupported platform: $(PLATFORM) :/.)
endif

TARGETS = $(BINDIR)/nlx2ft$(EXTENSION)

##############################################################################
all: $(TARGETS)

$(BINDIR)/%$(EXTENSION): %.o
	$(CXX) -o $@ $^ $(LDFLAGS) $(LDLIBS)

clean:
	$(RM) *.o $(call fixpath, $(TARGETS))
