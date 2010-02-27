##############################################################################
#
# GNU'ish Makefile
#
# $Log: Makefile,v $
# Revision 1.7  2009/05/13 20:46:07  roboos
# added test application that writes as fast as it can
#
# Revision 1.6  2009/01/23 17:48:45  roboos
# added some test applications, changed dependencies to improve compile behaviour
#
# Revision 1.5  2008/11/14 15:48:34  roboos
# number of small changes, nothing significant
#
# Revision 1.4  2008/07/09 10:39:44  roboos
# cleaned up, ensure that all dependencies are matched
#
# Revision 1.3  2008/07/09 10:09:40  roboos
# renamed the various demo applications into demo_xxx
#
# Revision 1.2  2008/07/08 20:24:43  roboos
# added event example
#
# Revision 1.1  2008/07/08 18:52:30  roboos
# initial version, helps in compiling test and demo applications
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

all: test demo

test: test_gethdr test_getdat test_getevt test_flushhdr test_flushdat test_flushevt

demo: demo_sinewave demo_event demo_buffer demo_combined

demo_combined: demo_combined.o sinewave.o ../src/libbuffer.a
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

demo_sinewave: demo_sinewave.o sinewave.o ../src/libbuffer.a
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

demo_event: demo_event.o event.o ../src/libbuffer.a
	$(CC) $(LIBPATH) -o $(BINDIR)/$@ $^ $(LIBS)

demo_buffer: demo_buffer.o ../src/libbuffer.a
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
