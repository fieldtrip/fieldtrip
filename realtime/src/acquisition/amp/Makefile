 # Copyright (C) 2013 Federico Raimondo
 # Applied Artificial Intelligence Lab
 # Computer Sciences Department
 # University of Buenos Aires, Argentina
 
 # This file is part of Amp2ft
 #
 # Amp2ft is free software: you can redistribute it and/or modify
 # it under the terms of the GNU General Public License as published by
 # the Free Software Foundation, either version 3 of the License, or
 # (at your option) any later version.
 #
 # Amp2ft is distributed in the hope that it will be useful,
 # but WITHOUT ANY WARRANTY; without even the implied warranty of
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 # GNU General Public License for more details.
 #
 # You should have received a copy of the GNU General Public License
 # along with Amp2ft.  If not, see <http://www.gnu.org/licenses/>.

ifneq "$(OS)" "Windows_NT"
	OS = $(shell uname -s)
	MACHINE = $(shell uname -m)
endif

PLATFORM = $(shell gcc -dumpmachine)

FIELDTRIP = ../../../..
FTBUFFER = $(FIELDTRIP)/realtime/src/buffer

CXXFLAGS = -Wunused -Wall -pedantic -O3 -DDEBUG
INCPATH  = -I$(FTBUFFER)/src -I$(FTBUFFER)/cpp -I.
LIBPATH  = -L$(FTBUFFER)/src
LDLIBS   = -lbuffer -lpthread

# better would be to make a proper library of the C++ framework
CPPOBJS = $(patsubst %,$(FTBUFFER)/cpp/%, SignalConfiguration.o GdfWriter.o FtConnection.o StringServer.o)

ifeq "$(OS)" "Linux"
	fixpath = $1
	LDLIBS += -ldl -lrt
	ifeq "$(MACHINE)" "i686"
		BINDIR = $(FIELDTRIP)/realtime/bin/glnx86
	endif
	ifeq "$(MACHINE)" "x86_64"
		BINDIR = $(FIELDTRIP)/realtime/bin/glnxa64
	endif
endif

ifeq "$(OS)" "Darwin"
	fixpath = $1
	ifeq "$(MACHINE)" "i386"
		BINDIR    = $(FIELDTRIP)/realtime/bin/maci
		CFLAGS   += -m32
		CXXFLAGS += -m32
	endif
	ifeq "$(MACHINE)" "x86_64"
		BINDIR    = $(FIELDTRIP)/realtime/bin/maci64
		CFLAGS   += -m64
		CXXFLAGS += -m64
	endif
endif

ifndef BINDIR
$(error Unsupported platform: $(PLATFORM) :/.)
endif

TARGETS = $(BINDIR)/amp2ft

##############################################################################
all: $(TARGETS)

%.o: %.cpp
	$(CXX) $(INCPATH) $(CXXFLAGS) -c $<

AmpServerClient.o: AmpServerClient.h

$(BINDIR)/amp2ft: amp2ft.o AmpServerClient.o
	$(CXX) -o $@ $^ $(LIBPATH) $(LDLIBS) $(CPPOBJS)

clean:
	rm -f *.o $(TARGETS)
