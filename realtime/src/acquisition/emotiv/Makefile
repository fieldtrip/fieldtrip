# Copyright (C) 2010 Stefan Klanke
# Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
# Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
# Kapittelweg 29, 6525 EN Nijmegen, The Netherlands

# Use MinGW for compilation on Windows!
# To recompile, please download the Emotiv EDK, and update the EDKPATH variable below.

PLATFORM = $(shell gcc -dumpmachine)
UNAME = $(shell uname)
MACHINE = $(shell uname -m)

FIELDTRIP = ../../../..
FTBUFFER = $(FIELDTRIP)/realtime/src/buffer

# better would be to make a proper library of the C++ framework
CPPOBJ = $(patsubst %,$(FTBUFFER)/cpp/%,\
	SignalConfiguration.o GdfWriter.o FtConnection.o StringServer.o)

CXXFLAGS = -I$(FTBUFFER)/src -I$(FTBUFFER)/cpp -Wunused -Wall
LDFLAGS = -L$(FTBUFFER)/src 
LDLIBS = -lbuffer

ifeq "$(PLATFORM)" "mingw32"
	EDKPATH = C:\Program\ Files\ (x86)\Emotiv\ Education\ Edition\ SDK_v1.0.0.5-PREMIUM\doc\examples
	BINDIR = $(FIELDTRIP)/realtime/bin/win32
	SUFFIX = .exe
	CXXFLAGS += -I$(EDKPATH)/include
	LDFLAGS += -L$(EDKPATH)/lib
	LDLIBS += $(CPPOBJ) -lws2_32 -lwinmm $(BINDIR)/edk.dll $(BINDIR)/pthreadGC2.dll 
	# Make clean target work on windows:
	fixpath = $(subst /,\,$1)
	RM = del
endif


ifeq "$(UNAME)" "Linux"
	EDKPATH = /home/fraimondo/EmotivResearch_1.0.0.5
	CXX=g++
	LDLIBS += $(CPPOBJ) -ledk -ledk_utils
	LDFLAGS += -L$(EDKPATH)/lib
	CXXFLAGS += -I$(EDKPATH)/doc/examples_Qt/example5/
	ifeq "$(MACHINE)" "i686"
		BINDIR = $(FIELDTRIP)/realtime/bin/glnx86
    endif
    ifeq "$(MACHINE)" "x86_64"
		CXXFLAGS += -m32
		FORCE32 = 1
    	export FORCE32
		BINDIR = $(FIELDTRIP)/realtime/bin/glnxa64
    endif
endif

ifeq "$(UNAME)" "Darwin"
	EDKPATH = /Applications/EmotivResearch
	CXX=g++
	LDLIBS += $(CPPOBJ) -ledk
	LDFLAGS += -L$(EDKPATH)/docs/Examples/EmotivSDKLib/
	CXXFLAGS += -I$(EDKPATH)/docs/Examples/Qt/include/
    ifeq "$(MACHINE)" "x86_64"
		BINDIR = $(FIELDTRIP)/realtime/bin/glnxa64
    endif
endif

ifndef BINDIR
$(error Unsupported platform: $(PLATFORM) :/.)
endif

TARGETS = $(BINDIR)/emotiv2ft$(SUFFIX)

##############################################################################
all: $(TARGETS)

$(BINDIR)/%$(SUFFIX): %.cc
ifdef FORCE32
	cd $(FTBUFFER)/src && $(MAKE) clean
	cd $(FTBUFFER)/src && $(MAKE)
	cd $(FTBUFFER)/cpp && $(MAKE) clean
	cd $(FTBUFFER)/cpp && $(MAKE)
endif
	$(CXX) -o $@ $^ $(CPPOBJS) $(CXXFLAGS) $(LDFLAGS) $(LDLIBS)
ifdef FORCE32
	FORCE32=0
	export FORCE32
	cd $(FTBUFFER)/src && $(MAKE) clean
	cd $(FTBUFFER)/src && $(MAKE)
	cd $(FTBUFFER)/cpp && $(MAKE) clean
	cd $(FTBUFFER)/cpp && $(MAKE)
endif

clean:
	$(RM) $(TARGETS) $(call fixpath, $(TARGETS))
