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
CC		 = gcc
INCPATH	 = -I../src
LIBPATH  = -L../src
LIBS 	 = -lbuffer -lm -lpthread
SUFFIX   =
RM       = rm -f

ifeq "$(PLATFORM)" "mingw32"
	INCPATH	 = -I../src -I../pthreads-win32/include
	LIBPATH  = -L../src -L../pthreads-win32/lib
	LIBS 	 = -lbuffer -lm -lpthreadGC2 -lws2_32
	CFLAGS	 = $(INCPATH) -Wunused -pedantic -O0 -g
	SUFFIX   = .exe
	RM       = del
endif

ifeq "$(PLATFORM)" "i686-w64-mingw32"
	# prevent dependency on libwinpthread-1.dll by linking statically
	# see http://stackoverflow.com/questions/13768515/how-to-do-static-linking-of-libwinpthread-1-dll-in-mingw
	INCPATH	 = -I../src
	LIBPATH  = -L../src
	LIBS 	 = -lbuffer -lm -lpthread -lws2_32 -static
	CFLAGS	 = $(INCPATH) -Wunused -pedantic -Wall -O0 -g
	SUFFIX   = _win32.exe
	RM       = del
endif

ifeq "$(PLATFORM)" "x86_64-w64-mingw32"
	# prevent dependency on libwinpthread-1.dll by linking statically
	# see http://stackoverflow.com/questions/13768515/how-to-do-static-linking-of-libwinpthread-1-dll-in-mingw
	INCPATH	 = -I../src
	LIBPATH  = -L../src
	LIBS 	 = -lbuffer -lm -lpthread -lws2_32 -static
	CFLAGS	 = $(INCPATH) -Wunused -pedantic -Wall -O0 -g
	SUFFIX   = _win64.exe
	RM       = del
endif

ifeq "$(OS)" "Darwin"
    ifeq "$(MACHINE)" "i386"
        CFLAGS += -m32
        CXXFLAGS += -m32
        LDFLAGS += -m32
    endif
    ifeq "$(MACHINE)" "x86_64"
        CFLAGS += -m64
        CXXFLAGS += -m64
        LDFLAGS += -m64
    endif
endif

BINDIR	 = .
CFLAGS	 += $(INCPATH) -Wunused -pedantic -O0 -fPIC -g
CXXFLAGS += $(INCPATH) -Wunused -pedantic -O0

##############################################################################

all: demo test

demo: demo_combined$(SUFFIX) demo_sinewave$(SUFFIX) demo_event$(SUFFIX)

test: test_gethdr$(SUFFIX) test_getdat$(SUFFIX) test_getevt$(SUFFIX) test_flushhdr$(SUFFIX) test_flushdat$(SUFFIX) test_flushevt$(SUFFIX) test_pthread$(SUFFIX) test_benchmark$(SUFFIX) test_nslookup$(SUFFIX) test_waitdat$(SUFFIX) test_connect$(SUFFIX)

demo_combined$(SUFFIX): demo_combined.o sinewave.o ../src/libbuffer.a
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS) $(LDFLAGS)

demo_sinewave$(SUFFIX): demo_sinewave.o sinewave.o ../src/libbuffer.a
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS) $(LDFLAGS)

demo_event$(SUFFIX): demo_event.o event.o ../src/libbuffer.a
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS) $(LDFLAGS)

test_gethdr$(SUFFIX): test_gethdr.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS) $(LDFLAGS)

test_getdat$(SUFFIX): test_getdat.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS) $(LDFLAGS)

test_getevt$(SUFFIX): test_getevt.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS) $(LDFLAGS)

test_flushhdr$(SUFFIX): test_flushhdr.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS) $(LDFLAGS)

test_flushdat$(SUFFIX): test_flushdat.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS) $(LDFLAGS)

test_flushevt$(SUFFIX): test_flushevt.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS) $(LDFLAGS)

test_pthread$(SUFFIX): test_pthread.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS) $(LDFLAGS)

test_benchmark$(SUFFIX): test_benchmark.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS) $(LDFLAGS)

test_nslookup$(SUFFIX): test_nslookup.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS) $(LDFLAGS)

test_waitdat$(SUFFIX): test_waitdat.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS) $(LDFLAGS)

test_connect$(SUFFIX): test_connect.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS) $(LDFLAGS)

%.o: %.c message.h buffer.h
	$(CC) $(CFLAGS) $(INCPATH) -c $*.c

clean:
	$(RM) core *.o *.obj *.a

distclean:
	$(RM) demo_combined$(SUFFIX) demo_sinewave$(SUFFIX) demo_event$(SUFFIX)
	$(RM) test_gethdr$(SUFFIX) test_getdat$(SUFFIX) test_getevt$(SUFFIX) test_flushhdr$(SUFFIX) test_flushdat$(SUFFIX) test_flushevt$(SUFFIX) test_pthread$(SUFFIX) test_benchmark$(SUFFIX) test_nslookup$(SUFFIX) test_waitdat$(SUFFIX) test_connect$(SUFFIX)
