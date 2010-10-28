# Copyright (C) 2010 Stefan Klanke
# Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
# Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
# Kapittelweg 29, 6525 EN Nijmegen, The Netherlands

# Use MinGW for compilation on Windows!
# This Makefile assumes that FLTK has been installed to a place where the compiler can find it.
# On Windows, this location is C:/msys/1.0/local -- change this if need be
CC = gcc
CXX = g++

PLATFORM = $(shell gcc -dumpmachine)

ifeq "$(PLATFORM)" "mingw32"
	FLPREFIX = C:/msys/1.0/local
	LDFLAGS  = pthreadGC2.dll ../buffer/src/libbuffer.a	
	FLFLAGS  = -I$(FLPREFIX)/include 
	LDFL     = -mno-cygwin -lfltk -lole32 -luuid -lcomctl32 -lgdi32 -lwsock32 -L$(FLPREFIX)/lib
	SUFFIX   = .exe
else
	LDFLAGS = -lpthread ../buffer/src/libbuffer.a
	FLFLAGS = `fltk-config --cxxflags`
	LDFL    = `fltk-config --ldflags`
	SUFFIX  = 
endif

CXXFLAGS = $(FLFLAGS) -I../buffer/src -I../buffer/cpp -I. -Wunused -Wall -pedantic

all: bufferViewer$(SUFFIX)
	
bufferViewer.o: bufferViewer.cc
	g++ -c bufferViewer.cc $(CXXFLAGS)
	
bufferViewer$(SUFFIX): bufferViewer.o ../buffer/cpp/FtConnection.o
	g++ -o bufferViewer$(SUFFIX) $^ $(LDFLAGS) $(LDFL)
