##############################################################################
#
# GNU'ish Makefile
#
# $Log: Makefile,v $
#
##############################################################################

PLATFORM = $(shell gcc -dumpmachine)
PLATFORM = i686-apple-darwin9

CC		 = gcc
BINDIR	 = .
INCPATH  = -I.
LIBPATH  = -L.
LIBS 	 = -lm -lpthread -lpeer
# CFLAGS   = $(INCPATH) -fPIC -Wall -pedantic -O1 -fpack-struct
# CXXFLAGS = $(INCPATH) -fPIC -Wunused -pedantic -O1
CFLAGS	 = $(INCPATH) -fPIC -Wall -pthread -g

ifeq "$(PLATFORM)" "i386-redhat-linux"
	# this is only for the peerslave command line executable -> mentat068
	MATLABARCH = glnx86
	MATLABPATH = /opt/matlab75
	MATLABLIBS = -L$(MATLABPATH)/bin/$(MATLABARCH) -leng -lmx
	MATLABINC  = -I$(MATLABPATH)/extern/include
endif

ifeq "$(PLATFORM)" "x86_64-redhat-linux"
	# this is only for the peerslave command line executable -> mentat001
	MATLABARCH = glnxa64
	MATLABPATH = /opt/matlab77
	MATLABLIBS = -L$(MATLABPATH)/bin/$(MATLABARCH) -leng -lmx
	MATLABINC  = -I$(MATLABPATH)/extern/include
endif

ifeq "$(PLATFORM)" "i686-apple-darwin9"
	# this is only for the peerslave command line executable -> manzana
	MATLABARCH = maci
	# matlab version 2009a is a 32 bit application
	MATLABPATH = /Applications/MATLAB_R2009a.app
	MATLABLIBS = -L$(MATLABPATH)/bin/$(MATLABARCH) -leng -lmx
	MATLABINC  = -I$(MATLABPATH)/extern/include
    # override the architecture defaults, required for 32 bit
    CFLAGS	 = $(INCPATH) -Wall -pthread -arch i386
endif

ifeq "$(PLATFORM)" "i686-apple-darwin10"
	# this is only for the peerslave command line executable -> macbook
	MATLABARCH = maci64
	MATLABPATH = /Applications/MATLAB_R2010a.app
	MATLABLIBS = -L$(MATLABPATH)/bin/$(MATLABARCH) -leng -lmx
	MATLABINC  = -I$(MATLABPATH)/extern/include
endif

ifeq "$(PLATFORM)" "mingw32"
	# this is only for the peerslave command line executable -> 32-bit Matlab on Windows 7 (DCCN511)
	MATLABARCH = exe
	MATLABPATH = "C:\Program Files (x86)\Matlab\R2009b32bit"
	MATLABLIBS = $(MATLABPATH)/bin/win32/libeng.dll $(MATLABPATH)/bin/win32/libmx.dll
	MATLABINC  = -I$(MATLABPATH)/extern/include
	LIBS 	 = -lm -lpeer ../pthreadGC2.dll -lws2_32 -lwininet
    # override the architecture defaults
    CFLAGS	 = $(INCPATH) -Wall -I../pthreads-win32/include -DSYSLOG=2 -DLOG_DEBUG=0 -DLOG_ERR=0
endif

ifeq "$(PLATFORM)" "x86_64-w64-mingw32"
	# this is only for the peerslave command line executable -> 32-bit Matlab on Windows 7 (DCCN511)
	MATLABARCH = w64.exe
	MATLABPATH = "C:\Program Files\Matlab\R2009b"
	MATLABLIBS = $(MATLABPATH)/bin/win64/libeng.dll $(MATLABPATH)/bin/win64/libmx.dll
	MATLABINC  = -I$(MATLABPATH)/extern/include
	LIBS 	 = -lm -lpeer ../pthreadGC2-w64.dll -lws2_32 -lwininet
    # override the architecture defaults
    CFLAGS	 = $(INCPATH) -Wall -I../pthreads-win64/include -DSYSLOG=2 -DLOG_DEBUG=0 -DLOG_ERR=0
endif

## to override the compilation defaults, please specify below where
## your MATLAB is installed and the architecture details
# MATLABARCH = --- string with your matlab platform architecture ---
# MATLABPATH = --- location of your matlab installation ---
# MATLABLIBS = -L$(MATLABPATH)/bin/$(MATLABARCH) -leng -lmx
# MATLABINC  = -I$(MATLABPATH)/extern/include
## to compile the 32-bit version with a 64-bit compiler, you will have to add the following
# CFLAGS += -arch i386

TARGET   = peerslave.$(MATLABARCH)

##############################################################################

all: $(TARGET)

libpeer.a: tcpserver.o udsserver.o tcpsocket.o discover.o announce.o expire.o util.o extern.o peerinit.o security.o localhost.o smartshare.o smartmem.o smartcpu.o connect.o
	ar rv $@ $^

peer: peer.o libpeer.a
	$(CC) $(CFLAGS) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

peerslave.$(MATLABARCH): peerslave.o libpeer.a
	$(CC) $(CFLAGS) $(LIBPATH) $(MATLABLIBS) -o ../$@ $^ $(LIBS)

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

%.o: %.c peer.h extern.h
	$(CC) $(CFLAGS) $(INCPATH) -c $*.c

peerslave.o: peerslave.c peer.h extern.h libpeer.a
	$(CC) $(CFLAGS) $(INCPATH) $(MATLABINC) -c $*.c

test_peerslave.o: test_peerslave.c peer.h extern.h libpeer.a
	$(CC) $(CFLAGS) $(INCPATH) $(MATLABINC) -c $*.c

test_engine.o: test_engine.c
	$(CC) $(CFLAGS) $(INCPATH) $(MATLABINC) -c $*.c

clean:
	rm -f core *.o *.obj *.a ../$(TARGET) *~

distclean:
	rm -f core *.o *.obj *.a ../$(TARGET) *~
	rm -f ../peerslave.glnx86
	rm -f ../peerslave.glnxa64
	rm -f ../peerslave.maci
	rm -f ../peerslave.maci64

%.o: %.c buffer.h message.h swapbytes.h socket_includes.h unix_includes.h
	$(CC) $(CFLAGS) -c $*.c

