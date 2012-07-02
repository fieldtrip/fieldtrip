# Copyright (C) 2010 Stefan Klanke
# Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
# Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
# Kapittelweg 29, 6525 EN Nijmegen, The Netherlands

# Use MinGW for compilation on Windows!

PLATFORM = $(shell gcc -dumpmachine)

#defaults, might be over-written further down
LDFLAGS += -L../src
LDFLAGS += -lpthread -lbuffer
INCPATH = -I../src -I. 
RM = rm -f

ifeq "$(PLATFORM)" "mingw32"
	LDFLAGS += ../../external/pthreads-win32/lib/pthreadGC2.dll \
		-LC:/msys/1.0/local/lib -lws2_32 -lwinmm 
	INCPATH += -I../../external/pthreads-win32/include
	SUFFIX = .exe
	RM = del
endif

ifeq "$(PLATFORM)" "x86_64-w64-mingw32"
	LDFLAGS += ../../external/pthreads-win64/lib/pthreadGC2-w64.dll \
		LC:/msys/1.0/local/lib -lws2_32 -lwinmm 
	INCPATH += -I../../external/pthreads-win64/include
	SUFFIX = .exe
	RM = del
endif

ifeq "$(PLATFORM)" "i686-apple-darwin10"
	LDFLAGS += -ldl
endif

ifeq "$(PLATFORM)" "x86_64-redhat-linux"
	LDFLAGS += -ldl
endif

CFLAGS += $(INCPATH) -Wunused -Wall -pedantic -O3

all: odmTest$(SUFFIX)

FtConnection.o: FtConnection.cc FtBuffer.h
	g++ -c FtConnection.cc $(CFLAGS)	
	
SignalConfiguration.o: SignalConfiguration.cc SignalConfiguration.h
	g++ -c SignalConfiguration.cc $(CFLAGS)	
	
StringServer.o: StringServer.cc StringServer.h
	g++ -c StringServer.cc $(CFLAGS)	
	
GdfWriter.o: GdfWriter.cc GdfWriter.h 
	g++ -c GdfWriter.cc $(CFLAGS)	
	
odmTest.o: odmTest.cc OnlineDataManager.h StringServer.o GDF_BackgroundWriter.h
	g++ -c odmTest.cc $(CFLAGS)

odmTest$(SUFFIX): odmTest.o SignalConfiguration.o GdfWriter.o FtConnection.o StringServer.o
	g++ -o $@ $^ $(LDFLAGS) 

clean:
	$(RM) *.o *.exe
