# Copyright (C) 2010 Stefan Klanke
# Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
# Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
# Kapittelweg 29, 6525 EN Nijmegen, The Netherlands

# Use MinGW for compilation on Windows!

PLATFORM = $(shell gcc -dumpmachine)
FIELDTRIP = ../../../..
FTBUFFER = $(FIELDTRIP)/realtime/src/buffer
CPPDIR  = $(FTBUFFER)/cpp

#defaults, might be over-written further down
LDFLAGS = -L$(FTBUFFER)/src
CXXFLAGS = -Wunused -Wall -pedantic -O3 -g
CXXFLAGS += -I$(FTBUFFER)/src -I. -I$(CPPDIR) 

ifeq "$(PLATFORM)" "mingw32"

	BINDIR = $(FIELDTRIP)/realtime/bin/win32
	CXXFLAGS += -I$(FTBUFFER)/pthreads-win32/include
	LDFLAGS += -LC:/msys/1.0/local/lib
	# LDFLAGS += -L$(BINDIR) does not seem to work for DLLs. Hence:
	LDLIBS += $(CPPOBJ) -lbuffer -lws2_32 -lwinmm $(BINDIR)/pthreadGC2.dll
	EXTENSION = .exe
	
	# Make clean target work on windows:
	fixpath = $(subst /,\,$1)
	RM = del
	CP = copy /Y
endif

ifeq "$(PLATFORM)" "x86_64-redhat-linux"
	LDFLAGS += -L$(FTBUFFER)/src
	LDLIBS += -lbuffer -lpthread -ldl
	BINDIR = $(FIELDTRIP)/realtime/bin/glnxa64
	EXTENSION = 
endif

ifndef BINDIR
$(warning Unsupported plaform: $(PLATFORM) :/.)
endif

# TODO: make proper library of C++ framework.
CPPOBJS = $(CPPDIR)/SignalConfiguration.o $(CPPDIR)/GdfWriter.o $(CPPDIR)/FtConnection.o $(CPPDIR)/StringServer.o

TARGETS = $(BINDIR)/biosemi2ft$(EXTENSION)

all: $(TARGETS)
	
$(BINDIR)/biosemi2ft$(EXTENSION): biosemi2ft.o BioSemiClient.o 
	$(CXX) -o $@ $^ $(CPPOBJS) $(LDFLAGS) $(LDLIBS)

clean:
	$(RM) *.o $(call fixpath, $(TARGETS))