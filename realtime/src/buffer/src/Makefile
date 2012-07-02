##############################################################################
#
# GNU'ish Makefile
#
# $Id$
##############################################################################

CC = gcc

PLATFORM = $(shell gcc -dumpmachine)

# default
INCPATH = -I.
LIBPATH = -L.
LIBS = -lbuffer -lm -lpthread
CFLAGS += $(INCPATH) -Wunused -Wall -pedantic -O0 -fPIC -g
CXXFLAGS += $(INCPATH) -Wunused -Wall -pedantic -O0 -fPIC -g

ifeq "$(PLATFORM)" "i686-apple-darwin10"
	RM = rm -f
endif

ifeq "$(PLATFORM)" "mingw32"
	INCPATH	+= -I. -I../../external/pthreads-win32/include
	LIBPATH += -L. -L../../external/pthreads-win32/lib
	LIBS += -lpthreadGC2
	RM = del
endif

ifeq "$(PLATFORM)" "x86_64-w64-mingw32"
	INCPATH	+= -I. -I../../external/pthreads-win64/include
	LIBPATH += -L. -L../../external/pthreads-win64/lib
	RM       = del
endif

BINDIR	 = .

#---

all: libbuffer.a

libbuffer.a: tcpserver.o tcpsocket.o tcprequest.o clientrequest.o dmarequest.o cleanup.o util.o printstruct.o swapbytes.o extern.o endianutil.o socketserver.o
	 ar rv $@ $^

libclient.a: tcprequest.o util.o
	 ar rv $@ $^

%.o: %.c buffer.h message.h extern.h swapbytes.h socket_includes.h unix_includes.h
	$(CC) $(CFLAGS) -c $*.c

clean:
	$(RM) *.o *.a *.obj core
