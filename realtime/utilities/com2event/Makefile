# Copyright (C) 2010 Stefan Klanke
# Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
# Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
# Kapittelweg 29, 6525 EN Nijmegen, The Netherlands

# Use MinGW for compilation on Windows!

PLATFORM = $(shell gcc -dumpmachine)

ifeq "$(PLATFORM)" "mingw32"
   LDFLAGS = pthreadGC2.dll ../buffer/src/libbuffer.a -LC:/msys/1.0/local/lib -lws2_32 -lwinmm 
   ADDINC = -I../buffer/pthreads-win32/include
   SUFFIX = .exe
else
   LDFLAGS = -lpthread ../buffer/src/libbuffer.a
   ADDINC = 
   SUFFIX = 
endif

CFLAGS = -I../buffer/src -I. $(ADDINC) -Wunused -Wall -pedantic

all: serial_event$(SUFFIX)

###############################################################################
# Platform-independent C libraries / utilities
###############################################################################
serial.o: serial.c serial.h
	gcc -c serial.c -I.
		
###############################################################################
# Platform-independent applications
###############################################################################
serial_event.o: serial_event.c serial.h
	gcc -c serial_event.c $(CFLAGS)
	
serial_event$(SUFFIX): serial_event.o serial.o
	gcc -o serial_event$(SUFFIX) serial_event.o serial.o $(LDFLAGS)
