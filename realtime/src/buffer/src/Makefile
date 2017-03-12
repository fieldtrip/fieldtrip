# Copyright (C) 2008, Robert Oostenveld
# Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
# Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
# Kapittelweg 29, 6525 EN Nijmegen, The Netherlands

ifneq "$(OS)" "Windows_NT"
	OS = $(shell uname -s)
	MACHINE = $(shell uname -m)
endif

PLATFORM = $(shell gcc -dumpmachine)

# defaults, might be overwritten further down
CC = gcc
INCPATH = -I.
LIBPATH = -L.
LIBS = -lbuffer -lm -lpthread

ifeq "$(PLATFORM)" "mingw32"
	INCPATH	+= -I. -I../../external/pthreads-win32/include
	LIBPATH += -L. -L../../external/pthreads-win32/lib
	LIBS += -lpthreadGC2
	RM = rm -f
endif

ifeq "$(PLATFORM)" "i686-w64-mingw32"
	RM = del
endif

ifeq "$(PLATFORM)" "x86_64-w64-mingw32"
	RM = del
endif

ifeq "$(OS)" "Darwin"
  # Snow Leopard 10.6 is the first Mac OS X version to ship with a 64 bit kernel, and obviously your
	# Mac has to have a 64 bit processor to utilize that kernel, but some 64 bit Mac's are not defaulting
	# to the 64 bit kernel. Furthermore, irrespective of the kernel, the gcc compiler is able to (cross)
    # compile objects, libraries and applications for 32 and 64 bit.
    ifeq "$(MACHINE)" "i386"
        # compile for 32 bit
        CFLAGS += -m32
        CXXFLAGS += -m32
        LDFLAGS += -m32
    endif
    ifeq "$(MACHINE)" "x86_64"
        # compile for 64 bit
        CFLAGS += -m64
        CXXFLAGS += -m64
        LDFLAGS += -m64
    endif
endif

#CFLAGS += -O0 -fPIC -g -Wunused -Wall -pedantic
CFLAGS += -O0 -fPIC -g -w
CFLAGS += $(INCPATH)

CXXFLAGS += -O0 -fPIC -g -Wunused -Wall -pedantic
CXXFLAGS += $(INCPATH)

##############################################################################
all: libbuffer.a

libbuffer.a: tcpserver.o tcpsocket.o tcprequest.o clientrequest.o dmarequest.o cleanup.o clock_gettime.o gettimeofday.o timestamp.o util.o printstruct.o swapbytes.o extern.o endianutil.o socketserver.o rdaserver.o
	 ar rv $@ $^

libclient.a: tcprequest.o util.o
	 ar rv $@ $^

%.o: %.c buffer.h message.h extern.h timestamp.h swapbytes.h socket_includes.h unix_includes.h
	$(CC) $(CFLAGS) -c $*.c

clean:
	$(RM) *.o *.a *.obj core
