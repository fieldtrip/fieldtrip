##############################################################################
#
# GNU'ish Makefile
#
# $Id$
#
##############################################################################

# locate the MATLAB command line executable, which will help to determine the path and architecture

ifeq ($(MATLABCMD),)
MATLABCMD := $(shell which matlab)
endif

ifeq ($(MATLABCMD),)
MATLABCMD := $(shell which matlab2011a)
endif

ifeq ($(MATLABCMD),)
MATLABCMD := $(shell which matlab2010b)
endif

ifeq ($(MATLABCMD),)
MATLABCMD := $(shell which matlab2010a)
endif

ifeq ($(MATLABCMD),)
MATLABCMD := $(shell which matlab2009b)
endif

ifeq ($(MATLABCMD),)
MATLABCMD := $(shell which matlab2009a)
endif

ifeq ($(MATLABCMD),)
MATLABCMD := $(shell which matlab2008b)
endif

ifeq ($(MATLABCMD),)
MATLABCMD := $(shell which matlab2008a)
endif

ifeq ($(MATLABCMD),)
MATLABCMD := $(shell which matlab79)
endif

ifeq ($(MATLABCMD),)
MATLABCMD := $(shell which matlab78)
endif

ifeq ($(MATLABCMD),)
MATLABCMD := $(shell which matlab77)
endif

ifeq ($(MATLABCMD),)
$(error Cannot locate your MATLAB installation, you should specify MATLABCMD=xxx)
endif

# determine the MATLAB architecture and path
MATLABARCH := $(shell $(MATLABCMD) -e | grep ARCH= | cut -f 2 -d =)
MATLABPATH := $(shell $(MATLABCMD) -e | grep MATLAB= | cut -f 2 -d =)
MATLABLIBS = -L$(MATLABPATH)/bin/$(MATLABARCH) -leng -lmx
MATLABINC  = -I$(MATLABPATH)/extern/include

CC		 = gcc
BINDIR	 = ../bin
INCPATH  = -I. $(MATLABINC)
LIBPATH  = -L. $(MATLABLIBS)
LIBS 	 = -lm -lpthread -lpeer
# CFLAGS   = $(INCPATH) -fPIC -Wall -pedantic -O1 -fpack-struct
# CXXFLAGS = $(INCPATH) -fPIC -Wunused -pedantic -O1
CFLAGS	 = $(INCPATH) -fPIC -Wall -pthread -DSYSLOG=2

# determine the operating system
PLATFORM := $(shell gcc -dumpmachine)

# override the defaults, required for cross compilation (e.g. 32 bit MATLAB on 64 bit OS)
ifeq "$(PLATFORM)" "i686-apple-darwin9"
ifeq "$(MATLABARCH)" "maci"
    CFLAGS := $(CFLAGS) -arch i386
endif
endif

ifeq "$(PLATFORM)" "i686-apple-darwin10"
ifeq "$(MATLABARCH)" "maci"
    CFLAGS := $(CFLAGS) -arch i386
endif
endif

ifeq "$(PLATFORM)" "mingw32"
	# this is only for the peerslave command line executable -> 32-bit Matlab on Windows 7 (DCCN511)
	MATLABARCH = exe
	MATLABPATH = "C:\Program Files (x86)\Matlab\R2009b32bit"
	MATLABLIBS = $(MATLABPATH)/bin/win32/libeng.dll $(MATLABPATH)/bin/win32/libmx.dll
	MATLABINC  = -I$(MATLABPATH)/extern/include
	LIBS 	 = -lm -lpeer ../pthreadGC2.dll -lws2_32 -lwininet
    CFLAGS	 = $(INCPATH) -Wall -I../pthreads-win32/include -DSYSLOG=2 -DLOG_DEBUG=0 -DLOG_ERR=0
endif

ifeq "$(PLATFORM)" "x86_64-w64-mingw32"
	# this is only for the peerslave command line executable -> 32-bit Matlab on Windows 7 (DCCN511)
	MATLABARCH = w64.exe
	MATLABPATH = "C:\Program Files\Matlab\R2009b"
	MATLABLIBS = $(MATLABPATH)/bin/win64/libeng.dll $(MATLABPATH)/bin/win64/libmx.dll
	MATLABINC  = -I$(MATLABPATH)/extern/include
	LIBS 	 = -lm -lpeer ../pthreadGC2-w64.dll -lws2_32 -lwininet
    CFLAGS	 = $(INCPATH) -Wall -I../pthreads-win64/include -DSYSLOG=2 -DLOG_DEBUG=0 -DLOG_ERR=0
endif

LD_LIBRARY_PATH := $(MATLABPATH)/bin/$(MATLABARCH)
export LD_LIBRARY_PATH

TARGET   = peerslave.$(MATLABARCH)

##############################################################################

all: $(TARGET)

libpeer.a: tcpserver.o udsserver.o tcpsocket.o discover.o announce.o expire.o util.o extern.o peerinit.o security.o localhost.o smartshare.o smartmem.o smartcpu.o connect.o parser.o
	ar rv $@ $^

peer: peer.o libpeer.a
	$(CC) $(CFLAGS) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

peerslave.$(MATLABARCH): peerslave.o libpeer.a
	$(CC) $(CFLAGS) $(LIBPATH) $(MATLABLIBS) -o $(BINDIR)/$@ $^ $(LIBS)

test_peerslave: test_peerslave.o
	$(CC) $(CFLAGS) $(LIBPATH) $(MATLABLIBS) -o $@ $^ $(LIBS)

test_engine: test_engine.o
	$(CC) $(CFLAGS) $(LIBPATH) $(MATLABLIBS) -o $@ $^ $(LIBS)

test_send: test_send.o libpeer.a
	$(CC) $(CFLAGS) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

test_recv: test_recv.o libpeer.a
	$(CC) $(CFLAGS) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

test_announce: test_announce.o libpeer.a
	$(CC) $(CFLAGS) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

test_discover: test_discover.o libpeer.a
	$(CC) $(CFLAGS) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

test_sender: test_sender.o
	$(CC) $(CFLAGS) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

test_listener: test_listener.o
	$(CC) $(CFLAGS) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

test_linkedlist: test_linkedlist.o
	$(CC) $(CFLAGS) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

test_localhost: test_localhost.o
	$(CC) $(CFLAGS) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

test_parser: test_parser.o
	$(CC) $(CFLAGS) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

%.o: %.c peer.h extern.h
	$(CC) $(CFLAGS) $(INCPATH) -c $*.c

peerslave.o: peerslave.c peer.h extern.h libpeer.a
	$(CC) $(CFLAGS) $(INCPATH) $(MATLABINC) -c $*.c

test_peerslave.o: test_peerslave.c peer.h extern.h libpeer.a
	$(CC) $(CFLAGS) $(INCPATH) $(MATLABINC) -c $*.c

test_engine.o: test_engine.c
	$(CC) $(CFLAGS) $(INCPATH) $(MATLABINC) -c $*.c

clean:
	rm -f core *.o *.obj *.a *~

distclean: clean
	rm -f ../peerslave.glnx86
	rm -f ../peerslave.glnxa64
	rm -f ../peerslave.maci
	rm -f ../peerslave.maci64

%.o: %.c buffer.h message.h swapbytes.h socket_includes.h unix_includes.h
	$(CC) $(CFLAGS) -c $*.c

