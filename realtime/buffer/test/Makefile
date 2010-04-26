##############################################################################
#
# Makefile for gcc   on Linux
# Makefile for gcc   on OSX
# Makefile for mingw on Windows
#
##############################################################################

CC		 = gcc

PLATFORM = $(shell gcc -dumpmachine)

ifeq "$(PLATFORM)" "mingw32"
	INCPATH	 = -I../src -I../pthreads-win32/include
	LIBPATH  = -L../src -L../pthreads-win32/lib
	LIBS 	 = -lbuffer -lm -lpthreadGC2 -lws2_32
	CFLAGS	 = $(INCPATH) -Wunused -pedantic -O0 -fpack-struct -g 
else
	INCPATH	 = -I../src
	LIBPATH  = -L../src
	LIBS 	 = -lbuffer -lm -lpthread
	CFLAGS	 = $(INCPATH) -Wunused -pedantic -O0 -fPIC -fpack-struct -g
endif

BINDIR	 = .
CXXFLAGS = $(INCPATH) -Wunused -pedantic -O0

##############################################################################

all: demo

test: test_gethdr test_getdat test_getevt test_flushhdr test_flushdat test_flushevt test_connect test_nslookup test_benchmark test_waitdat

demo: demo_sinewave demo_event demo_buffer demo_combined

demo_combined: demo_combined.o sinewave.o ../src/libbuffer.a
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

demo_sinewave: demo_sinewave.o sinewave.o ../src/libbuffer.a
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

demo_event: demo_event.o event.o ../src/libbuffer.a
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

demo_buffer: demo_buffer.o ../src/libbuffer.a
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

demo_buffer_rda: demo_buffer_rda.o ../src/libbuffer.a ../src/rdaserver.o
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
else
	rm -f core *.o *.obj *.a *~
endif

distclean:
ifeq "$(PLATFORM)" "mingw32"
	del test_gethdr.exe test_getdat.exe test_getevt.exe test_flushhdr.exe test_flushdat.exe test_flushevt.exe
	del demo_sinewave.exe demo_event.exe demo_buffer.exe demo_combined.exe
else
	rm test_gethdr test_getdat test_getevt test_flushhdr test_flushdat test_flushevt
	rm demo_sinewave demo_event demo_buffer demo_combined
endif
