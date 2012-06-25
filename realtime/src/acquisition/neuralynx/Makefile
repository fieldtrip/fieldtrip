# Copyright (C) 2010 Stefan Klanke
# Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
# Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
# Kapittelweg 29, 6525 EN Nijmegen, The Netherlands

# Use MinGW for compilation on Windows!
CC = gcc
CXX = g++

LDFLAGS = pthreadGC2.dll ../../buffer/src/libbuffer.a -lws2_32 
CXXFLAGS = -I../../buffer/src -I../../buffer/cpp -I. -Wunused -Wall -pedantic -O3 -I../../buffer/pthreads-win32/include

all: nlx2ft.exe
	
nlx2ft.o: nlx2ft.cc
	g++ -c nlx2ft.cc $(CXXFLAGS)
	
nlx2ft.exe: nlx2ft.o
	g++ -o nlx2ft.exe nlx2ft.o $(LDFLAGS) MatlabNetComClient.dll
