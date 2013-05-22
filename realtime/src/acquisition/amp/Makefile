# Configure fieldtrip path

FIELDTRIP = ../../../..


# Do not touch!
PLATFORM = $(shell gcc -dumpmachine)
UNAME = $(shell uname)
MACHINE = $(shell uname -m)

FTBUFFER = $(FIELDTRIP)/realtime/src/buffer
CPPDIR  = $(FTBUFFER)/cpp
CPPOBJS = $(CPPDIR)/SignalConfiguration.o $(CPPDIR)/GdfWriter.o \
	$(CPPDIR)/FtConnection.o $(CPPDIR)/StringServer.o
LDFLAGS = -L$(FTBUFFER)/src
LDLIBS = -lbuffer $(CPPOBJS)
CXXFLAGS = -Wunused -Wall -pedantic -O3 -g
CXXFLAGS += -I$(FTBUFFER)/src -I. -I$(CPPDIR) -DDEBUG
CXX=/usr/bin/g++

ifeq "$(UNAME)" "Darwin"
	LDLIBS += -lpthread
    fixpath = $1
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
endif

ifeq "$(UNAME)" "Linux"
	ifeq "$(MACHINE)" "i686"
		LDFLAGS += -L$(FTBUFFER)/src
		LDLIBS += -lbuffer -lpthread -ldl
		BINDIR = $(FIELDTRIP)/realtime/bin/glnx89
    endif
    ifeq "$(MACHINE)" "x86_64"
		LDFLAGS += -L$(FTBUFFER)/src
		LDLIBS += -lbuffer -lpthread -ldl
		BINDIR = $(FIELDTRIP)/realtime/bin/glnxa64
    endif
	LDFLAGS += -L$(FTBUFFER)/src
	LDLIBS += -lbuffer -lpthread -ldl
	BINDIR = $(FIELDTRIP)/realtime/bin/glnxa64
	BINDIR = .
endif

ifndef BINDIR
$(warning Unsupported platform: $(PLAT) :/.)
endif

TARGETS = $(BINDIR)/amp2ft
all: $(TARGETS)
	
$(BINDIR)/amp2ft: amp2ft.o AmpServerClient.o AmpServerClient.h
	$(CXX) -o $@ $^ $(LDFLAGS) $(LDLIBS) $(CXXFLAGS)

clean:
	rm -f *.o $(TARGETS)
