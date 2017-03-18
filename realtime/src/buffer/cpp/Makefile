# Copyright (C) 2010 Stefan Klanke
# Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
# Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
# Kapittelweg 29, 6525 EN Nijmegen, The Netherlands

# Use MinGW for compilation on Windows!

ifneq "$(OS)" "Windows_NT"
	OS = $(shell uname -s)
	MACHINE = $(shell uname -m)
endif

PLATFORM = $(shell gcc -dumpmachine)

# defaults, might be overwritten further down
CC       = gcc
CXX      = g++
CFLAGS   = -Wunused -Wall -pedantic -O3
CXXFLAGS = -Wunused -Wall -pedantic -O3
LDFLAGS  = -L../src
LDLIBS   = -lpthread -lbuffer
INCPATH  = -I. -I../src
SUFFIX   =

ifeq "$(PLATFORM)" "mingw32"
  INCPATH += -I../../external/pthreads-win32/include
	LDFLAGS += ../../external/pthreads-win32/lib/pthreadGC2.dll -LC:/msys/1.0/local/lib -lws2_32 -lwinmm
	SUFFIX   = .exe
	RM       = del
endif

ifeq "$(PLATFORM)" "i686-w64-mingw32"
	# prevent dependency on libwinpthread-1.dll by linking statically
	# see http://stackoverflow.com/questions/13768515/how-to-do-static-linking-of-libwinpthread-1-dll-in-mingw
	LDFLAGS += -lws2_32 -lwinmm -static
	SUFFIX   = _win32.exe
	RM       = del
endif

ifeq "$(PLATFORM)" "x86_64-w64-mingw32"
	# prevent dependency on libwinpthread-1.dll by linking statically
	# see http://stackoverflow.com/questions/13768515/how-to-do-static-linking-of-libwinpthread-1-dll-in-mingw
	LDFLAGS += -lws2_32 -lwinmm -static
	SUFFIX   = _win64.exe
	RM       = del
endif

ifeq "$(OS)" "Linux"
	ifeq "$(MACHINE)" "i686"
		LDLIBS += -ldl -lrt
	endif
	ifeq "$(MACHINE)" "x86_64"
		LDLIBS += -ldl -lrt
	endif
	ifeq "$(MACHINE)" "armv6l"
		LDLIBS += -ldl -lrt
	endif
	ifeq "$(MACHINE)" "armv7l"
		LDLIBS += -ldl -lrt
	endif
endif

ifeq "$(OS)" "Darwin"
	LDFLAGS += -ldl
	ifeq "$(MACHINE)" "i386"
		CFLAGS   += -m32
		CXXFLAGS += -m32
		LDFLAGS  += -m32
	endif
	ifeq "$(MACHINE)" "x86_64"
		CFLAGS   += -m64
		CXXFLAGS += -m64
		LDFLAGS  += -m64
	endif
endif

includes = $(wildcard *.h)

##############################################################################
all: odmTest$(SUFFIX)

%.o: %.cc ${includes}
	$(CXX) $(CXXFLAGS) $(INCPATH) -c $<

odmTest$(SUFFIX): odmTest.o SignalConfiguration.o GdfWriter.o FtConnection.o StringServer.o
	$(CXX) -o $@ $^ $(LDFLAGS) $(LDLIBS)

clean:
	$(RM) *.o odmTest$(SUFFIX)
