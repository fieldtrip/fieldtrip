##############################################################################
#
# Makefile for gcc   on Linux
# Makefile for gcc   on OSX
# Makefile for mingw on Windows
#
##############################################################################

CC		 = gcc

PLATFORM = $(shell gcc -dumpmachine)

# Default options for POSIX systems
INCPATH	 = -I../src
LIBPATH  = -L../src
LIBS 	 = -lbuffer -lm -lpthread
CFLAGS	 = $(INCPATH) -Wunused -pedantic -O0 -fPIC -g
SUFFIX   =

ifeq "$(PLATFORM)" "mingw32"
	INCPATH	 = -I../src -I../pthreads-win32/include
	LIBPATH  = -L../src -L../pthreads-win32/lib
	LIBS 	 = -lbuffer -lm -lpthreadGC2 -lws2_32
	CFLAGS	 = $(INCPATH) -Wunused -pedantic -O0 -g
	SUFFIX   = .exe
	
endif

ifeq "$(PLATFORM)" "x86_64-w64-mingw32"
    # override the defaults
	INCPATH	 = -I../src -I../pthreads-win64/include
	LIBPATH  = -L../src -L../pthreads-win64/lib
	LIBS 	 = -lbuffer -lm -lpthread -lws2_32
	CFLAGS	 = $(INCPATH) -Wunused -pedantic -Wall -O0 -g
	SUFFIX   = _64.exe
endif

BINDIR	 = .
CXXFLAGS = $(INCPATH) -Wunused -pedantic -O0

##############################################################################

all: demo

test: test_gethdr test_getdat test_getevt test_flushhdr test_flushdat test_flushevt test_connect test_nslookup test_benchmark test_waitdat

demo: demo_sinewave$(SUFFIX) demo_event$(SUFFIX) demo_buffer$(SUFFIX) demo_combined$(SUFFIX) demo_buffer_unix$(SUFFIX)

demo_combined$(SUFFIX): demo_combined.o sinewave.o ../src/libbuffer.a
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

demo_sinewave$(SUFFIX): demo_sinewave.o sinewave.o ../src/libbuffer.a
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

demo_event$(SUFFIX): demo_event.o event.o ../src/libbuffer.a
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

demo_buffer$(SUFFIX): demo_buffer.o ../src/libbuffer.a
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

demo_buffer_rda$(SUFFIX): demo_buffer_rda.o ../src/libbuffer.a ../src/rdaserver.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)	
	
demo_buffer_unix$(SUFFIX): demo_buffer_unix.o ../src/libbuffer.a ../src/socketserver.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)	

test_gethdr: test_gethdr.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

test_getdat: test_getdat.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

test_getevt: test_getevt.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

test_flushhdr: test_flushhdr.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

test_flushdat: test_flushdat.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

test_flushevt: test_flushevt.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

test_append: test_append.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

test_pthread: test_pthread.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

test_benchmark: test_benchmark.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

test_nslookup: test_nslookup.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)
	
test_waitdat: test_waitdat.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

test_connect: test_connect.o
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

%.o: %.c message.h buffer.h
	$(CC) $(CFLAGS) $(INCPATH) -c $*.c

clean:
ifeq "$(PLATFORM)" "mingw32"
	del *.o *.a
else ifeq "$(PLATFORM)" "x86_64-w64-mingw32"
	del *.o *.a
else
	rm -f core *.o *.obj *.a *~
endif

distclean:
ifeq "$(PLATFORM)" "mingw32"
	del test_gethdr.exe test_getdat.exe test_getevt.exe test_flushhdr.exe test_flushdat.exe test_flushevt.exe
	del demo_sinewave.exe demo_event.exe demo_buffer.exe demo_combined.exe
else ifeq "$(PLATFORM)" "x86_64-w64-mingw32"
	del test_gethdr.exe test_getdat.exe test_getevt.exe test_flushhdr.exe test_flushdat.exe test_flushevt.exe
	del demo_sinewave_64.exe demo_event_64.exe demo_buffer_64.exe demo_combined_64.exe demo_buffer_unix_64.exe
else
	rm test_gethdr test_getdat test_getevt test_flushhdr test_flushdat test_flushevt
	rm demo_sinewave demo_event demo_buffer demo_combined
endif
