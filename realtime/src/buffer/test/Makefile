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
CC       = gcc
CFLAGS   = -O2 -fPIC -w # -Wunused -Wall -pedantic
CXXFLAGS = -O2 -fPIC -w # -Wunused -Wall -pedantic
INCPATH  = -I../src
LIBPATH  = -L../src
LDLIBS   = -lbuffer -lm -lpthread
RM       = rm -f
BINDIR   = .
SUFFIX   =

ifeq "$(PLATFORM)" "mingw32"
  INCPATH += -I../pthreads-win32/include
  LIBPATH += -L../pthreads-win32/lib
  LDLIBS   = -lbuffer -lm -lpthreadGC2 -lws2_32
  SUFFIX   = .exe
  RM       = del
endif

ifeq "$(PLATFORM)" "i686-w64-mingw32"
  # prevent dependency on libwinpthread-1.dll by linking statically
  # see http://stackoverflow.com/questions/13768515/how-to-do-static-linking-of-libwinpthread-1-dll-in-mingw
  LDLIBS   = -lbuffer -lm -lpthread -lws2_32 -static
  SUFFIX   = _win32.exe
  RM       = del
endif

ifeq "$(PLATFORM)" "x86_64-w64-mingw32"
  # prevent dependency on libwinpthread-1.dll by linking statically
  # see http://stackoverflow.com/questions/13768515/how-to-do-static-linking-of-libwinpthread-1-dll-in-mingw
  LDLIBS   = -lbuffer -lm -lpthread -lws2_32 -static
  SUFFIX   = _win64.exe
  RM       = del
endif

ifeq "$(OS)" "Darwin"
    ifeq "$(MACHINE)" "i386"
        CFLAGS   += -m32
        CXXFLAGS += -m32
    endif
    ifeq "$(MACHINE)" "x86_64"
        CFLAGS   += -m64
        CXXFLAGS += -m64
    endif
endif

##############################################################################

all: demo test

demo: demo_combined$(SUFFIX) demo_sinewave$(SUFFIX) demo_event$(SUFFIX)

test: test_gethdr$(SUFFIX) test_getdat$(SUFFIX) test_getevt$(SUFFIX) test_flushhdr$(SUFFIX) test_flushdat$(SUFFIX) test_flushevt$(SUFFIX) test_pthread$(SUFFIX) test_benchmark$(SUFFIX) test_nslookup$(SUFFIX) test_waitdat$(SUFFIX) test_connect$(SUFFIX)

demo_combined$(SUFFIX): demo_combined.o sinewave.o ../src/libbuffer.a
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBPATH) $(LDLIBS)

demo_sinewave$(SUFFIX): demo_sinewave.o sinewave.o ../src/libbuffer.a
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBPATH) $(LDLIBS)

demo_event$(SUFFIX): demo_event.o event.o ../src/libbuffer.a
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBPATH) $(LDLIBS)

test_gethdr$(SUFFIX): test_gethdr.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBPATH) $(LDLIBS)

test_getdat$(SUFFIX): test_getdat.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBPATH) $(LDLIBS)

test_getevt$(SUFFIX): test_getevt.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBPATH) $(LDLIBS)

test_flushhdr$(SUFFIX): test_flushhdr.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBPATH) $(LDLIBS)

test_flushdat$(SUFFIX): test_flushdat.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBPATH) $(LDLIBS)

test_flushevt$(SUFFIX): test_flushevt.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBPATH) $(LDLIBS)

test_pthread$(SUFFIX): test_pthread.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBPATH) $(LDLIBS)

test_benchmark$(SUFFIX): test_benchmark.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBPATH) $(LDLIBS)

test_nslookup$(SUFFIX): test_nslookup.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBPATH) $(LDLIBS)

test_waitdat$(SUFFIX): test_waitdat.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBPATH) $(LDLIBS)

test_connect$(SUFFIX): test_connect.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBPATH) $(LDLIBS)

%.o: %.c
	$(CC) $(CFLAGS) $(INCPATH) -c $<

clean:
	$(RM) core *.o *.obj *.a

distclean:
	$(RM) demo_combined$(SUFFIX) demo_sinewave$(SUFFIX) demo_event$(SUFFIX)
	$(RM) test_gethdr$(SUFFIX) test_getdat$(SUFFIX) test_getevt$(SUFFIX) test_flushhdr$(SUFFIX) test_flushdat$(SUFFIX) test_flushevt$(SUFFIX) test_pthread$(SUFFIX) test_benchmark$(SUFFIX) test_nslookup$(SUFFIX) test_waitdat$(SUFFIX) test_connect$(SUFFIX)
