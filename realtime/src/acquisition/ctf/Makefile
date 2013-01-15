# Copyright (C) 2007 Robert Oostenveld
# Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
# Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
# Kapittelweg 29, 6525 EN Nijmegen, The Netherlands

PLATFORM = $(shell gcc -dumpmachine)
FIELDTRIP = ../../../..
FTBUFFER = $(FIELDTRIP)/realtime/src/buffer

CFLAGS = -I$(FTBUFFER)/src -Wunused -Wall
LDFLAGS = -L$(FTBUFFER)/src 
LDLIBS = -lbuffer -lncurses -lpthread -lm

ifneq "$(PLATFORM)" "i386-redhat-linux"
$(warning For the DCCN, please compile for 32-bit linux.)
endif

ifeq "$(PLATFORM)" "x86_64-linux-gnu"
	BINDIR = $(FIELDTRIP)/realtime/bin/glnxa64
endif

ifeq "$(PLATFORM)" "x86_64-redhat-linux"
	BINDIR = $(FIELDTRIP)/realtime/bin/glnxa64
endif

ifeq "$(PLATFORM)" "i386-redhat-linux"
	BINDIR = $(FIELDTRIP)/realtime/bin/glnx86
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
