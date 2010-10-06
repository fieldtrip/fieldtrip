##############################################################################
#
# GNU'ish Makefile
#
# $Log: Makefile,v $
# Revision 1.20  2009/01/23 08:26:44  roboos
# fixed a serious bug that caused a lot of memory to leak (in fact all packets that were sent over the socket would eventually leack away), both on the client and server side
#
# Revision 1.19  2008/11/14 15:48:34  roboos
# number of small changes, nothing significant
#
# Revision 1.18  2008/07/09 13:34:57  roboos
# compiling the buffer mex file on matlab75/64bt linux requires -fPIC
#
# Revision 1.17  2008/07/09 10:39:44  roboos
# cleaned up, ensure that all dependencies are matched
#
# Revision 1.16  2008/07/08 18:53:50  release
# *** empty log message ***
#
# Revision 1.15  2008/07/08 18:40:01  release
# cleanup, removed demo stuff, only keep core functions for libbuffer and libclient
#
##############################################################################

CC = gcc

PLATFORM = $(shell gcc -dumpmachine)

# default
INCPATH	 = -I.
LIBPATH  = -L.
LIBS 	 = -lbuffer -lm -lpthread
CFLAGS	 = $(INCPATH) -Wunused -Wall -pedantic -O0 -fPIC -g

ifeq "$(PLATFORM)" "i686-apple-darwin10"
    # override the defaults
	INCPATH	 = -I.
	LIBPATH  = -L.
	LIBS 	 = -lbuffer -lm -lpthread
	CFLAGS	 = $(INCPATH) -Wunused -Wall -pedantic -O0 -fPIC -g
endif

ifeq "$(PLATFORM)" "mingw32"
    # override the defaults
	INCPATH	 = -I. -I../pthreads-win32/include
	LIBPATH  = -L. -L../pthreads-win32/lib
	LIBS 	 = -lbuffer -lm -lpthreadGC2
	CFLAGS	 = $(INCPATH) -Wunused -pedantic -Wall -O0 -g
endif

ifeq "$(PLATFORM)" "x86_64-w64-mingw32"
    # override the defaults
	INCPATH	 = -I. -I../pthreads-win64/include
	LIBPATH  = -L. -L../pthreads-win64/lib
	LIBS 	 = -lbuffer -lm -lpthread
	CFLAGS	 = $(INCPATH) -Wunused -pedantic -Wall -O0 -g 
endif

BINDIR	 = .
CXXFLAGS = $(INCPATH) -Wunused -pedantic -O0

##############################################################################

all: libbuffer.a

libbuffer.a: tcpserver.o tcpsocket.o tcprequest.o clientrequest.o dmarequest.o cleanup.o util.o printstruct.o swapbytes.o extern.o endianutil.o 
	 ar rv $@ $^

libclient.a: tcprequest.o util.o
	 ar rv $@ $^

%.o: %.c buffer.h message.h extern.h swapbytes.h socket_includes.h unix_includes.h
	$(CC) $(CFLAGS) -c $*.c

clean:
ifeq "$(PLATFORM)" "mingw32"
	del *.o *.a
else ifeq "$(PLATFORM)" "x86_64-w64-mingw32"
	del *.o *.a
else
	rm -f core *.o *.obj *.a *~
endif
