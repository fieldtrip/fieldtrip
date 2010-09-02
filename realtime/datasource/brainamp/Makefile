# Copyright (C) 2010 Stefan Klanke
# Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
# Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
# Kapittelweg 29, 6525 EN Nijmegen, The Netherlands

# Use MinGW for compilation on Windows!

PLATFORM = $(shell gcc -dumpmachine)

ifeq "$(PLATFORM)" "mingw32"
   LDFLAGS = pthreadGC2.dll ../../buffer/src/libbuffer.a -LC:/msys/1.0/local/lib -lws2_32 -lwinmm 
   INCPATH = -I../../buffer/src -I. -I../../buffer/pthreads-win32/include
   SUFFIX = .exe
else
   LDFLAGS = -lpthread ../../buffer/src/libbuffer.a
   INCPATH = -I../../buffer/src -I. 
   SUFFIX = 
endif

CFLAGS = $(INCPATH) -Wunused -Wall -pedantic

all: rda2ft$(SUFFIX)
		
rda2ft.o: rda2ft.c
	gcc -c rda2ft.c $(CFLAGS)
	
rda2ft$(SUFFIX): rda2ft.o ../../buffer/src/socketserver.o
	gcc -o $@ $^ $(LDFLAGS)
