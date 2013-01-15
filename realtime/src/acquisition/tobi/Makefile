# Copyright (C) 2012 Boris Reudering
# Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
# Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
# Kapittelweg 29, 6525 EN Nijmegen, The Netherlands

# Use MinGW for compilation on Windows!

PLATFORM = $(shell gcc -dumpmachine)
UNAME   = $(shell uname)
MACHINE = $(shell uname -m)

FIELDTRIP = ../../../..
FTBUFFER = $(FIELDTRIP)/realtime/src/buffer

# defaults, might be overwritten further down
CFLAGS = -I$(FTBUFFER)/src -I. -Iinclude
CFLAGS += -Wunused -Wall -pedantic -O3
CXXFLAGS = $(CFLAGS)
LDFLAGS = -L$(FTBUFFER)/src

ifeq "$(PLATFORM)" "x86_64-linux-gnu"
	BINDIR = $(FIELDTRIP)/realtime/bin/glnxa64
	LDFLAGS += -L$(BINDIR) -Llib/ticpp/linux 
	LDLIBS += -lbuffer -lpthread -lboost_thread -lboost_system -lboost_program_options -ltia -lticpp_64 
endif

ifndef BINDIR
$(error Unsupported platform: $(PLATFORM) :/.)
endif

TARGETS = $(patsubst %,$(BINDIR)/%, tia2ft)

##############################################################################
all: $(TARGETS)

$(BINDIR)/%$(SUFFIX): %.o
	$(CXX) $^ $(LDFLAGS) $(LDLIBS) -o $@

$(BINDIR)/tia2ft$(SUFFIX): ft_buffer_helper.o tia2ft.o

clean:
	rm -f *.o $(LOCAL_TARGETS) $(TARGETS)

