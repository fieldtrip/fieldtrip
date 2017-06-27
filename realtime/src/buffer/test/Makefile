# Copyright (C) 2008, Robert Oostenveld
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
FTBUFFER  = $(FIELDTRIP)/realtime/src/buffer

# defaults, might be overwritten further down
CC       = gcc
CFLAGS   = -O2 -fPIC -w # -Wunused -Wall -pedantic
CXXFLAGS = -O2 -fPIC -w # -Wunused -Wall -pedantic
INCPATH  = -I../src
LIBPATH  = -L../src
LDLIBS   = -lbuffer -lm -lpthread
RM       = rm -f
SUFFIX   =

ifeq "$(PLATFORM)" "i686-pc-cygwin"
	BINDIR   = $(FIELDTRIP)/realtime/bin/win32
	LDLIBS   = -lbuffer -lm -lws2_32
	SUFFIX   = .exe
	RM       = rm -f
	fixpath  = $1
endif

ifeq "$(PLATFORM)" "x86_64-pc-cygwin"
	BINDIR   = $(FIELDTRIP)/realtime/bin/win64
	LDLIBS   = -lbuffer -lm -lws2_32
	SUFFIX   = .exe
	RM       = rm -f
	fixpath  = $1
endif

ifeq "$(PLATFORM)" "i686-w64-mingw32"
	BINDIR   = $(FIELDTRIP)/realtime/bin/win32
	# prevent dependency on libwinpthread-1.dll by linking statically
	# see http://stackoverflow.com/questions/13768515/how-to-do-static-linking-of-libwinpthread-1-dll-in-mingw
	LDLIBS   = -lbuffer -lm -lpthread -lws2_32 -static
	SUFFIX   = .exe
	RM       = del
	# make clean target work on windows:
	fixpath = $(subst /,\,$1)
endif

ifeq "$(PLATFORM)" "x86_64-w64-mingw32"
	BINDIR   = $(FIELDTRIP)/realtime/bin/win64
	# prevent dependency on libwinpthread-1.dll by linking statically
	# see http://stackoverflow.com/questions/13768515/how-to-do-static-linking-of-libwinpthread-1-dll-in-mingw
	LDLIBS   = -lbuffer -lm -lpthread -lws2_32 -static
	SUFFIX   = .exe
	RM       = del
	# make clean target work on windows:
	fixpath = $(subst /,\,$1)
endif

ifeq "$(OS)" "Linux"
	fixpath = $1
	LDLIBS += -ldl -lpthread -lportaudio
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

ifeq "$(OS)" "Darwin"
	fixpath = $1
	ifeq "$(MACHINE)" "i386"
		BINDIR    = $(FIELDTRIP)/realtime/bin/maci
		CFLAGS   += -m32
		CXXFLAGS += -m32
		LDFLAGS  += -m32
	endif
	ifeq "$(MACHINE)" "x86_64"
		BINDIR    = $(FIELDTRIP)/realtime/bin/maci64
		CFLAGS   += -m64
		CXXFLAGS += -m64
		LDFLAGS  += -m64
	endif
endif

ifndef BINDIR
$(error Unsupported platform: $(PLATFORM) :/.)
endif

TARGETS = $(patsubst %, $(BINDIR)/%$(SUFFIX), demo_combined demo_sinewave demo_event test_gethdr test_getdat test_getevt test_flushhdr test_flushdat test_flushevt test_pthread test_benchmark test_nslookup test_waitdat test_connect)

##############################################################################

all: demo test

demo: demo_combined$(SUFFIX) demo_sinewave$(SUFFIX) demo_event$(SUFFIX)

test: test_gethdr$(SUFFIX) test_getdat$(SUFFIX) test_getevt$(SUFFIX) test_flushhdr$(SUFFIX) test_flushdat$(SUFFIX) test_flushevt$(SUFFIX) test_pthread$(SUFFIX) test_benchmark$(SUFFIX) test_nslookup$(SUFFIX) test_waitdat$(SUFFIX) test_connect$(SUFFIX)

demo_combined$(SUFFIX): demo_combined.o sinewave.o ../src/libbuffer.a
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LDFLAGS) $(LIBPATH) $(LDLIBS)

demo_sinewave$(SUFFIX): demo_sinewave.o sinewave.o ../src/libbuffer.a
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LDFLAGS) $(LIBPATH) $(LDLIBS)

demo_event$(SUFFIX): demo_event.o event.o ../src/libbuffer.a
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LDFLAGS) $(LIBPATH) $(LDLIBS)

test_gethdr$(SUFFIX): test_gethdr.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LDFLAGS) $(LIBPATH) $(LDLIBS)

test_getdat$(SUFFIX): test_getdat.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LDFLAGS) $(LIBPATH) $(LDLIBS)

test_getevt$(SUFFIX): test_getevt.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LDFLAGS) $(LIBPATH) $(LDLIBS)

test_flushhdr$(SUFFIX): test_flushhdr.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LDFLAGS) $(LIBPATH) $(LDLIBS)

test_flushdat$(SUFFIX): test_flushdat.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LDFLAGS) $(LIBPATH) $(LDLIBS)

test_flushevt$(SUFFIX): test_flushevt.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LDFLAGS) $(LIBPATH) $(LDLIBS)

test_pthread$(SUFFIX): test_pthread.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LDFLAGS) $(LIBPATH) $(LDLIBS)

test_benchmark$(SUFFIX): test_benchmark.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LDFLAGS) $(LIBPATH) $(LDLIBS)

test_nslookup$(SUFFIX): test_nslookup.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LDFLAGS) $(LIBPATH) $(LDLIBS)

test_waitdat$(SUFFIX): test_waitdat.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LDFLAGS) $(LIBPATH) $(LDLIBS)

test_connect$(SUFFIX): test_connect.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LDFLAGS) $(LIBPATH) $(LDLIBS)

%.o: %.c
	$(CC) $(CFLAGS) $(INCPATH) -c $<

clean:
	$(RM) core *.o *.obj *.a $(call fixpath, $(TARGETS))
