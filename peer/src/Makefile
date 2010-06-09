##############################################################################
#
# GNU'ish Makefile
#
# $Log: Makefile,v $
#
##############################################################################

CC		 = gcc
BINDIR	 = .
INCPATH  = -I.
LIBPATH  = -L.
LIBS 	 = -lm -lpthread -lpeer
CFLAGS	 = $(INCPATH) -Wunused -pedantic -O1 -fPIC -fpack-struct
# CXXFLAGS = $(INCPATH) -Wunused -pedantic -O1
CFLAGS	 = $(INCPATH) -Wall -pthread
TEST = test_announce test_discover test_send test_recv

##############################################################################

all: $(TEST)

libpeer.a: tcpserver.o tcpsocket.o discover.o announce.o expire.o util.o extern.o peerinit.o fairshare.o security.o localhost.o
	ar rv $@ $^

peer: peer.o libpeer.a
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

test_send: test_send.o libpeer.a
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

test_recv: test_recv.o libpeer.a
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

test_announce: test_announce.o libpeer.a
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

test_discover: test_discover.o libpeer.a
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

sender: sender.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

listener: listener.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

linkedlist: linkedlist.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

%.o: %.c peer.h extern.h
	$(CC) $(CFLAGS) $(INCPATH) -c $*.c

clean:
	rm -f core *.o *.obj *.a $(TEST) *~

distclean:

%.o: %.c buffer.h message.h swapbytes.h socket_includes.h unix_includes.h
	$(CC) $(CFLAGS) -c $*.c

