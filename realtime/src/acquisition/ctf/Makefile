PLATFORM = $(shell gcc -dumpmachine)
FTPATH = ../../../..
FTBUFFER = $(FTPATH)/realtime/src/buffer

CFLAGS = -I$(FTBUFFER)/src -Wunused -Wall
LDFLAGS = -L$(FTBUFFER)/src 
LDLIBS = -lbuffer -lncurses -lpthread

ifneq "$(PLATFORM)" "x86_32-linux-gnu"
$(warning For the DCCN, please compile on 32-bit linux. You are using an untested platform.)
endif

ifeq "$(PLATFORM)" "x86_64-linux-gnu"
	BINDIR = $(FTPATH)/realtime/bin/glnxa64
	EXTENSION = 
endif

TARGETS = $(patsubst %,$(BINDIR)/%$(EXTENSION),\
	AcqBuffer acq2ft acq2ftx fake_meg check_shm acq2ftx_fake acq2ft_fake)

all: $(TARGETS)

$(BINDIR)/%$(EXTENSION): %.c
	$(CC) $< $(CFLAGS) $(LDFLAGS) $(LDLIBS) -o $@

$(BINDIR)/AcqBuffer$(EXTENSION): AcqBuffer.c AcqBuffer.h

$(BINDIR)/acq2ftx$(EXTENSION): LDLIBS += -lm

$(BINDIR)/acq2ft_fake$(EXTENSION): acq2ft.c
	$(CC) $^ $(CFLAGS) -DFAKE_MEG $(LDFLAGS) $(LDLIBS) -o $@
	
$(BINDIR)/acq2ftx_fake$(EXTENSION): acq2ftx.c
	$(CC) $^ $(CFLAGS) -DFAKE_MEG $(LDFLAGS) $(LDLIBS) -lm -o $@

clean:
	rm -f *.o
	rm -f $(TARGETS)
