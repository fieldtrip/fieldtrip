# Copyright (C) 2008, Robert Oostenveld
# Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
# Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
# Kapittelweg 29, 6525 EN Nijmegen, The Netherlands

PLATFORM = $(shell gcc -dumpmachine)
UNAME = $(shell uname)
MACHINE = $(shell uname -m)

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

ifeq "$(PLATFORM)" "x86_64-w64-mingw32"
	INCPATH	 = -I../src -I../pthreads-win64/include
	LIBPATH  = -L../src -L../pthreads-win64/lib
	LIBS 	 = -lbuffer -lm -lpthread -lws2_32
	CFLAGS	 = $(INCPATH) -Wunused -pedantic -Wall -O0 -g
	SUFFIX   = _64.exe
	RM       = del
endif

ifeq "$(UNAME)" "Darwin"
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

all: demo

test: test_gethdr test_getdat test_getevt test_flushhdr test_flushdat test_flushevt test_connect test_nslookup test_benchmark test_waitdat

demo: demo_sinewave$(SUFFIX) demo_event$(SUFFIX) demo_combined$(SUFFIX)

demo_combined$(SUFFIX): demo_combined.o sinewave.o ../src/libbuffer.a
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS) $(LDFLAGS)

demo_sinewave$(SUFFIX): demo_sinewave.o sinewave.o ../src/libbuffer.a
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS) $(LDFLAGS)

demo_event$(SUFFIX): demo_event.o event.o ../src/libbuffer.a
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS) $(LDFLAGS)

demo_buffer$(SUFFIX): demo_buffer.o ../src/libbuffer.a
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS) $(LDFLAGS)

demo_buffer_rda$(SUFFIX): demo_buffer_rda.o ../src/libbuffer.a ../src/rdaserver.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS) $(LDFLAGS)	
	
demo_buffer_unix$(SUFFIX): demo_buffer_unix.o ../src/libbuffer.a 
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS) $(LDFLAGS)	

test_gethdr: test_gethdr.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS) $(LDFLAGS)

test_getdat: test_getdat.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS) $(LDFLAGS)

test_getevt: test_getevt.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS) $(LDFLAGS)

test_flushhdr: test_flushhdr.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS) $(LDFLAGS)

test_flushdat: test_flushdat.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS) $(LDFLAGS)

test_flushevt: test_flushevt.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS) $(LDFLAGS)

test_append: test_append.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS) $(LDFLAGS)

test_pthread: test_pthread.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS) $(LDFLAGS)

test_benchmark: test_benchmark.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS) $(LDFLAGS)

test_nslookup: test_nslookup.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS) $(LDFLAGS)
	
test_waitdat: test_waitdat.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS) $(LDFLAGS)

test_connect: test_connect.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS) $(LDFLAGS)

%.o: %.c message.h buffer.h
	$(CC) $(CFLAGS) $(INCPATH) -c $*.c

clean:
	$(RM) core *.o *.obj *.a 

distclean:
	$(RM) test_gethdr$(SUFFIX) test_getdat$(SUFFIX) test_getevt$(SUFFIX) test_flushhdr$(SUFFIX) test_flushdat$(SUFFIX) test_flushevt$(SUFFIX)
	$(RM) rm demo_sinewave$(SUFFIX) demo_event$(SUFFIX) demo_buffer$(SUFFIX) demo_buffer_unix$(SUFFIX) demo_combined$(SUFFIX) 
