# Copyright (C) 2010 Stefan Klanke
# Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
# Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
# Kapittelweg 29, 6525 EN Nijmegen, The Netherlands

# Use MinGW for compilation on Windows!

PLATFORM = $(shell gcc -dumpmachine)
CPPDIR   = ../../buffer/cpp

LDFLAGS = -lpthread ../../buffer/src/socketserver.o ../../buffer/src/libbuffer.a -ldl 
INCPATH = -I../../buffer/src -I. -I$(CPPDIR)
SUFFIX = 

ifeq "$(PLATFORM)" "mingw32"
   LDFLAGS = pthreadGC2.dll ../../buffer/src/socketserver.o ../../buffer/src/libbuffer.a -LC:/msys/1.0/local/lib -lws2_32 -lwinmm 
   INCPATH = -I../../buffer/src -I. -I$(CPPDIR) -I../../buffer/pthreads-win32/include
   SUFFIX = .exe
endif

ifeq "$(PLATFORM)" "i686-apple-darwin10"
   LDFLAGS = -lpthread ../../buffer/src/socketserver.o ../../buffer/src/libbuffer.a -ldl 
   INCPATH = -I../../buffer/src -I. -I$(CPPDIR)
   SUFFIX = .maci64
endif

CFLAGS = $(INCPATH) -Wunused -Wall -pedantic -O3 -g

all: modeeg2ft$(SUFFIX)  modeeg2ft_2chn$(SUFFIX) setupbluesmirf$(SUFFIX)

serial.o: serial.c
	gcc -c serial.c -I.
	
modeeg2ft.o: modeeg2ft.cc
	g++ -c modeeg2ft.cc $(CFLAGS)
	
modeeg2ft$(SUFFIX): modeeg2ft.o serial.o
	g++ -o $@ $^ $(CPPDIR)/SignalConfiguration.o $(CPPDIR)/GdfWriter.o $(CPPDIR)/FtConnection.o $(CPPDIR)/StringServer.o $(LDFLAGS) 

modeeg2ft_2chn.o: modeeg2ft_2chn.cc
	g++ -c modeeg2ft_2chn.cc $(CFLAGS)
	
modeeg2ft_2chn$(SUFFIX): modeeg2ft_2chn.o serial.o
	g++ -o $@ $^ $(LDFLAGS) $(CFLAGS) $(CPPDIR)/GdfWriter.cc

setupbluesmirf$(SUFFIX): setupbluesmirf.c serial.o
	gcc -o $@ $^ -I.

clean:
	rm -f *.o
	rm -f modeeg2ft$(SUFFIX)
	rm -f modeeg2ft_2chn$(SUFFIX)
	rm -f setupbluesmirf$(SUFFIX)

