# Copyright (C) 2010 Bart Niehuis, Stefan Klanke
# Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
# Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
# Kapittelweg 29, 6525 EN Nijmegen, The Netherlands

# Use MinGW for compilation on Windows!

PLATFORM = $(shell gcc -dumpmachine)

ifeq "$(PLATFORM)" "mingw32"
   LDFLAGS = pthreadGC2.dll ../../buffer/src/libbuffer.a -lws2_32
   INCPATH = -I. -I../../buffer/src -I../../buffer/pthreads-win32/include
   SUFFIX = .exe
else 
   LDFLAGS = -lpthread ../../buffer/src/libbuffer.a -ldl 
   INCPATH =  -I.  -I../../buffer/src
   SUFFIX = 
endif

CFLAGS = $(INCPATH) -Wunused -Wall -pedantic -O3 -g

all: tmsidriver$(SUFFIX)

tmsidriver.o: tmsidriver.cpp 
	g++ -c $^ $(CFLAGS)

	
tmsidriver$(SUFFIX): tmsidriver.o RTDevice.o Feature.o
	g++ -o $@ $^ $(LDFLAGS) 

clean:
	del *.o