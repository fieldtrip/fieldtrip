PLATFORM = $(shell gcc -dumpmachine)
FTPATH = ../../../..
FTBUFFER = $(FTPATH)/realtime/src/buffer

CFLAGS = -I$(FTBUFFER)/src -Wunused -Wall
LDFLAGS = -L$(FTBUFFER)/src 
LDLIBS = -lbuffer -lncurses -lpthread

ifneq "$(PLATFORM)" "x86_32-linux-gnu"
$(warning For the DCCN, please compile for 32-bit linux.)
endif

ifeq "$(PLATFORM)" "i686-apple-darwin10"
	BINDIR = $(FTPATH)/realtime/bin/maci
	EXTENSION = .maci
endif

ifeq "$(PLATFORM)" "x86_64-linux-gnu"
	BINDIR = $(FTPATH)/realtime/bin/glnxa64
	EXTENSION = 
endif

ifndef BINDIR
$(error Unsupported platform: $(PLATFORM) :/.)
endif


TARGETS = $(patsubst %,$(BINDIR)/ctf_%$(EXTENSION),\
	AcqBuffer acq2ft acq2ftx fake_meg check_shm acq2ftx_fake acq2ft_fake)

all: $(TARGETS)

$(BINDIR)/ctf_%$(EXTENSION): %.o
	$(CC) $^ $(LDFLAGS) $(LDLIBS) -o $@

$(BINDIR)/acq2ftx$(EXTENSION): LDLIBS += -lm

%_fake.o: %.c
	$(CC) $^ $(CFLAGS) -DFAKE_MEG -c -o $@

clean:
	rm -f *.o $(TARGETS)
