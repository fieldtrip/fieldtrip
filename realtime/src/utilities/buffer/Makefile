# Copyright (C) 2010 Stefan Klanke
# Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
# Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
# Kapittelweg 29, 6525 EN Nijmegen, The Netherlands

# Use MinGW for compilation on Windows!

PLATFORM = $(shell gcc -dumpmachine)
FIELDTRIP = ../../../..
FTBUFFER = $(FIELDTRIP)/realtime/src/buffer

#defaults, might be over-written further down
CFLAGS = -Wunused -Wall -pedantic -O3 -g -I$(FTBUFFER)/src
LDFLAGS = -L$(FTBUFFER)/src
LDLIBS = -lbuffer

ifeq "$(PLATFORM)" "mingw32"
	BINDIR = $(FIELDTRIP)/realtime/bin/win32
	EXTENSION = .exe
	CXXFLAGS += -I$(FIELDTRIP)/realtime/src/external/pthreads-win32/include
	# LDFLAGS += -L$(BINDIR) does not seem to work for DLLs. Hence:
	LDLIBS += -lws2_32 -lwinmm $(BINDIR)/pthreadGC2.dll

	# Make clean target work on windows:
	fixpath = $(subst /,\,$1)
	RM = del
endif

ifeq "$(PLATFORM)" "i686-apple-darwin10"
	# Currently, only 32 bit is supported since the 64 bit drivers don't seem to
	# work. This means that for correct compilation of biosemi2ft, also libbuffer
	# and ODM need to be (manually) recompiled and linked with '-m32' flag in the
	# CFLAGS and LDFLAGS.
	BINDIR = $(FIELDTRIP)/realtime/bin/maci
	CFLAGS += -m32
	LDLIBS += -lpthread -m32
	fixpath = $1
endif

ifeq "$(PLATFORM)" "x86_64-redhat-linux"
	LDLIBS += -lpthread -ldl
	BINDIR = $(FIELDTRIP)/realtime/bin/glnxa64
endif

ifeq "$(PLATFORM)" "x86_64-linux-gnu"
	LDLIBS += -lpthread -ldl
	BINDIR = $(FIELDTRIP)/realtime/bin/glnxa64
endif

ifndef BINDIR
$(warning Unsupported platform: $(PLATFORM) :/.)
endif

TARGETS = $(patsubst %,$(BINDIR)/%$(EXTENSION),\
	demo_buffer demo_buffer_rda)  # demo_buffer_tia does not work.

all: $(TARGETS)

$(BINDIR)/demo_buffer_rda$(EXTENSION): \
	LDLIBS := $(FTBUFFER)/src/rdaserver.o $(LDLIBS)
	
$(BINDIR)/%$(EXTENSION): %.o
	$(CC) -o $@ $^ $(LDFLAGS) $(LDLIBS)

clean:
	$(RM) *.o $(call fixpath, $(TARGETS))
