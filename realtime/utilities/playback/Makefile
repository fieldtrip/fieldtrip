# Copyright (C) 2010 Stefan Klanke
# Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
# Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
# Kapittelweg 29, 6525 EN Nijmegen, The Netherlands

# Use MinGW for compilation on Windows!

PLATFORM = $(shell gcc -dumpmachine)

LDFLAGS = -lpthread ../../buffer/src/libbuffer.a
INCPATH = -I../../buffer/src -I. 
SUFFIX = 

ifeq "$(PLATFORM)" "mingw32"
   LDFLAGS = pthreadGC2.dll ../../buffer/src/libbuffer.a -LC:/msys/1.0/local/lib -lws2_32 -lwinmm 
   INCPATH = -I../../buffer/src -I. -I../../buffer/pthreads-win32/include
   SUFFIX = .exe
endif

CFLAGS = $(INCPATH) -Wunused -Wall -pedantic -g

all: playback$(SUFFIX) recording$(SUFFIX)

playback.o: playback.c
	gcc -c playback.c $(CFLAGS)
	
recording.o: recording.c
	gcc -c $^ $(CFLAGS)
	
ft_offline.o: ft_offline.c
	gcc -c $^ $(CFLAGS)
			
playback$(SUFFIX): playback.o
	gcc -o $@ $^ $(LDFLAGS)
		
recording$(SUFFIX): recording.o ft_offline.o ../../buffer/src/socketserver.o
	gcc -o $@ $^ $(LDFLAGS)	

clean:
	rm *.o
