PLATFORM = $(shell gcc -dumpmachine)
FTPATH = ../../../..
FTBUFFER = $(FTPATH)/realtime/src/buffer

CFLAGS = -I$(FTBUFFER)/src -Wunused -Wall
LDFLAGS = -L$(FTBUFFER)/src 
LDLIBS = -lbuffer -lncurses -lpthread -lm

ifneq "$(PLATFORM)" "i386-redhat-linux"
$(warning For the DCCN, please compile for 32-bit linux.)
endif

ifeq "$(PLATFORM)" "x86_64-linux-gnu"
	BINDIR = $(FTPATH)/realtime/bin/glnxa64
endif

ifeq "$(PLATFORM)" "x86_64-redhat-linux"
	BINDIR = $(FTPATH)/realtime/bin/glnxa64
endif

ifeq "$(PLATFORM)" "i386-redhat-linux"
	BINDIR = $(FTPATH)/realtime/bin/glnx86
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
