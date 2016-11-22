# Copyright (C) 2007 Robert Oostenveld
# Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
# Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
# Kapittelweg 29, 6525 EN Nijmegen, The Netherlands

PLATFORM = $(shell gcc -dumpmachine)
UNAME = $(shell uname)
MACHINE = $(shell uname -m)

FIELDTRIP = ../../../..
FTBUFFER = $(FIELDTRIP)/realtime/src/buffer

CFLAGS = -I$(FTBUFFER)/src -Wunused -Wall
LDFLAGS = -L$(FTBUFFER)/src 
LDLIBS = -lbuffer -lncurses -lpthread -lm

ifeq "$(UNAME)" "Linux"
	LDLIBS += -lrt
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

ifeq "$(UNAME)" "Darwin"
	ifeq "$(MACHINE)" "i386"
		BINDIR = $(FIELDTRIP)/realtime/bin/maci
		CFLAGS += -m32
		CXXFLAGS += -m32
		LDFLAGS += -m32
	endif
	ifeq "$(MACHINE)" "x86_64"
		BINDIR = $(FIELDTRIP)/realtime/bin/maci64
		CFLAGS += -m64
		CXXFLAGS += -m64
		LDFLAGS += -m64
	endif
	ifeq "$(MACHINE)" "Power Macintosh"
		BINDIR = $(FIELDTRIP)/realtime/bin/mac
	endif
endif

ifneq "$(PLATFORM)" "i386-redhat-linux"
$(warning For the DCCN, please compile for 32-bit linux.)
endif

ifndef BINDIR
$(error Unsupported platform: $(PLATFORM) :/.)
endif

TARGETS = $(patsubst %,$(BINDIR)/ctf_%$(SUFFIX),\
	AcqBuffer acq2ft acq2ftx fake_meg check_shm acq2ftx_fake acq2ft_fake)

##############################################################################
all: $(TARGETS)

$(BINDIR)/ctf_%$(SUFFIX): %.o
	$(CC) -o $@ $^ $(LDFLAGS) $(LDLIBS)

$(BINDIR)/acq2ftx$(SUFFIX): LDLIBS += -lm

%_fake.o: %.c
	$(CC) $^ $(CFLAGS) -DFAKE_MEG -c -o $@

clean:
	$(RM) *.o $(TARGETS)
