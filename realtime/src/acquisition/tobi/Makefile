# Copyright (C) 2012 Boris Reudering
# Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
# Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
# Kapittelweg 29, 6525 EN Nijmegen, The Netherlands

# Use MinGW for compilation on Windows

ifneq "$(OS)" "Windows_NT"
	OS      ?= $(shell uname -s)
	MACHINE ?= $(shell uname -m)
endif

PLATFORM ?= $(shell gcc -dumpmachine)

FIELDTRIP = ../../../..
FTBUFFER = $(FIELDTRIP)/realtime/src/buffer

# defaults, might be overwritten further down
CC       = gcc
CXX      = g++
CFLAGS   = -O2 -w # -Wunused -Wall -pedantic
CXXFLAGS = -O2 -w # -Wunused -Wall -pedantic
INCPATH  = -I$(FTBUFFER)/src -Iinclude -I.
LIBPATH  = -L$(FTBUFFER)/src

ifeq "$(OS)" "Linux"
	ifeq "$(MACHINE)" "x86_64"
		BINDIR   = $(FIELDTRIP)/realtime/bin/glnxa64
		LIBPATH += -L$(BINDIR) -Llib/ticpp/linux
		LDLIBS  += -lbuffer -lpthread -lboost_thread -lboost_system -lboost_program_options -ltia -lticpp_64
	endif
endif

ifndef BINDIR
$(error Unsupported platform: $(PLATFORM) :/.)
endif

TARGETS = $(patsubst %,$(BINDIR)/%, tia2ft)

##############################################################################
all: $(TARGETS)

$(BINDIR)/%$(SUFFIX): %.o
	$(CXX) -o $@ $^ $(LDFLAGS) $(LIBPATH) $(LDLIBS)

$(BINDIR)/tia2ft$(SUFFIX): ft_buffer_helper.o tia2ft.o

clean:
	rm -f *.o $(TARGETS)
