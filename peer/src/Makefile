##############################################################################
#
# GNU'ish Makefile
#
# $Log: Makefile,v $
#
##############################################################################

PLATFORM = $(shell gcc -dumpmachine)

CC		 = gcc
BINDIR	 = .
INCPATH  = -I.
LIBPATH  = -L.
LIBS 	 = -lm -lpthread -lpeer
CFLAGS	 = $(INCPATH) -Wunused -pedantic -O1 -fPIC -fpack-struct
# CXXFLAGS = $(INCPATH) -Wunused -pedantic -O1
CFLAGS	 = $(INCPATH) -Wall -pthread

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
	MATLABPATH = /Applications/MATLAB-2009a
	MATLABLIBS = -L$(MATLABPATH)/bin/$(MATLABARCH) -leng -lmx
	MATLABINC  = -I$(MATLABPATH)/extern/include
    # override the architecture defaults
    CFLAGS	 = $(INCPATH) -Wall -pthread -arch i386
endif

ifeq "$(PLATFORM)" "i686-apple-darwin10"
	# this is only for the peerslave command line executable -> macbook
	MATLABARCH = maci
	MATLABPATH = /Applications/MATLAB_R2009a.app
	MATLABLIBS = -L$(MATLABPATH)/bin/$(MATLABARCH) -leng -lmx
	MATLABINC  = -I$(MATLABPATH)/extern/include
    # override the architecture defaults
    CFLAGS	 = $(INCPATH) -Wall -pthread -arch i386
endif

TARGET   = test_announce test_discover test_send test_recv peerslave.$(MATLABARCH)

##############################################################################

all: $(TARGET)

libpeer.a: tcpserver.o udsserver.o tcpsocket.o discover.o announce.o expire.o util.o extern.o peerinit.o fairshare.o security.o localhost.o smartmem.o connect.o
	ar rv $@ $^

peer: peer.o libpeer.a
	$(CC) $(CFLAGS) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

peerslave.$(MATLABARCH): peerslave.o libpeer.a
	$(CC) $(CFLAGS) $(LIBPATH) $(MATLABLIBS) -o ../$@ $^ $(LIBS)

test_send: test_send.o libpeer.a
	$(CC) $(CFLAGS) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

test_recv: test_recv.o libpeer.a
	$(CC) $(CFLAGS) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

test_announce: test_announce.o libpeer.a
	$(CC) $(CFLAGS) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

test_discover: test_discover.o libpeer.a
	$(CC) $(CFLAGS) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

sender: sender.o
	$(CC) $(CFLAGS) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

listener: listener.o
	$(CC) $(CFLAGS) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

linkedlist: linkedlist.o
	$(CC) $(CFLAGS) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

%.o: %.c peer.h extern.h
	$(CC) $(CFLAGS) $(INCPATH) -c $*.c

peerslave.o: peerslave.c peer.h extern.h libpeer.a
	$(CC) $(CFLAGS) $(INCPATH) $(MATLABINC) -c $*.c

clean:
	rm -f core *.o *.obj *.a $(TARGET) *~

distclean:

%.o: %.c buffer.h message.h swapbytes.h socket_includes.h unix_includes.h
	$(CC) $(CFLAGS) -c $*.c

