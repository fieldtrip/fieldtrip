# Use MinGW for compilation on Windows!

PLATFORM = $(shell gcc -dumpmachine)
FIELDTRIP = ../../../..
FTBUFFER = $(FIELDTRIP)/realtime/src/buffer

#defaults, might be over-written further down
INCPATH = -I$(FTBUFFER)/src -I. -Iinclude
LDFLAGS = -L$(FTBUFFER)/src
CFLAGS = -c $(INCPATH) -Wunused -Wall -pedantic -O3
CXXFLAGS = $(CFLAGS)

ifeq "$(PLATFORM)" "x86_64-linux-gnu"
	LDFLAGS += -Llib/tia/amd64 -Llib/ticpp/linux 
	LDLIBS += -lpthread -lboost_thread -lboost_system -lboost_program_options\
		-ltia -lticpp_64 -lbuffer 
	BINDIR = $(FIELDTRIP)/realtime/bin/glnxa64
endif

ifndef BINDIR
$(error Unsupported platform: $(PLATFORM) :/.)
endif

LOCAL_TARGETS = $(patsubst %,%$(EXTENSION), tia2ft)
TARGETS = $(patsubst %,$(BINDIR)/%, $(LOCAL_TARGETS))

all: $(TARGETS)

$(BINDIR)/%$(EXTENSION): %
	cp $< $@

tia2ft: ft_buffer_helper.o tia2ft.o

clean:
	rm -f *.o $(LOCAL_TARGETS) $(TARGETS)

