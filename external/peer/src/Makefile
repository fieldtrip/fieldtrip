##############################################################################
#
# GNU'ish Makefile
#
# $Log: Makefile,v $
# Revision 1.2  2010/01/10 00:01:21  roboos
# Many small changes, added some helper functions in util.c, obsoleted
# cleanup.c.  Added three lists for checking whether a connection (in
# tcpsocket) should be accepted and whether an announcement (in
# discover, see also expire) should be accepted: based on user, group
# and hostname. For each there is an allowed list.  The group can be
# specified by the user, the user and hostname are determined from
# the OS.
#
# The allowed-hostlist checking should be improved, because it now
# doe literal checking of the hostname as a string, whereas it should
# look at IP addresses.
#
# Revision 1.1  2010/01/06 16:23:12  roboos
# moved all source code to subdirectory
# updated compile script
# use random peer in peerfeval instead of the first
# retry in peerfeval in case the peerput fails
#
# Revision 1.10  2010/01/05 16:56:59  roboos
# updated makefile, improved clean
#
# Revision 1.9  2010/01/03 13:08:06  roboos
# fixed various small bugs which were found after a proper test on multiple mentats
#
# Revision 1.8  2010/01/02 11:02:12  roboos
# thread cancelation seems to work, all memory is deallocated
# "clear all" does not retuls in segfault any more
#
# Revision 1.7  2010/01/02 00:12:28  roboos
# first partially functional implementation of thread cancelation
# the mex file cannot be unloaded yet
# the 'bind' state is not properly cleaned upon unloading of the discover thread
#
# Revision 1.6  2009/12/31 08:25:00  roboos
# intermediate commit, many changes
# this version is reasonably functional in matlab, can get peerlist and start all threads
#
# Revision 1.5  2009/12/26 10:28:44  roboos
# some small changes, tryinmg to get detection of host and ip working
#
# Revision 1.4  2009/12/16 22:01:54  roboos
# implemented the expire thread, made seperate tests for announce and discover
#
# Revision 1.3  2009/12/14 21:28:23  roboos
# many small and large changes
# the code compiles, the test application starts three threads
# the announce thread works
# the discover thread works, but does not yet create a linked list
# have not yet looked at the tcpserver thread, nor at the tcpsocket threads that it should spawn
#
# Revision 1.2  2009/12/13 21:53:15  roboos
# many small changes, partially compiles
#
# Revision 1.1.1.1  2009/12/13 20:31:38  roboos
# initial import
#
#
##############################################################################

CC		 = gcc
BINDIR	 = .
INCPATH  = -I.
LIBPATH  = -L.
LIBS 	 = -lm -lpthread -lpeer
# CFLAGS	 = $(INCPATH) -Wunused -pedantic -O1 -fPIC -fpack-struct
# CXXFLAGS = $(INCPATH) -Wunused -pedantic -O1
CFLAGS	 = $(INCPATH) -Wall -pthread
TEST = test_announce test_discover test_send test_recv

##############################################################################

all: $(TEST)

libpeer.a: tcpserver.o tcpsocket.o discover.o announce.o expire.o util.o extern.o peerinit.o
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

