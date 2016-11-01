 # Copyright (C) 2013 Federico Raimondo
 # Applied Artificial Intelligence Lab
 # Computer Sciences Department 
 # University of Buenos Aires, Argentina
 #
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

FIELDTRIP = ../../../..


# Do not touch!
PLATFORM = $(shell gcc -dumpmachine)
UNAME = $(shell uname)
MACHINE = $(shell uname -m)

FTBUFFER = $(FIELDTRIP)/realtime/src/buffer
CPPDIR  = $(FTBUFFER)/cpp
CPPOBJS = $(CPPDIR)/SignalConfiguration.o $(CPPDIR)/GdfWriter.o \
	$(CPPDIR)/FtConnection.o $(CPPDIR)/StringServer.o
LDFLAGS = -L$(FTBUFFER)/src
LDLIBS = -lbuffer $(CPPOBJS)
CXXFLAGS = -Wunused -Wall -pedantic -O3 -g
CXXFLAGS += -I$(FTBUFFER)/src -I. -I$(CPPDIR) -DDEBUG
CXX = /usr/bin/g++

ifeq "$(UNAME)" "Darwin"
	LDLIBS += -lpthread
    fixpath = $1
	ifeq "$(MACHINE)" "i386"
        BINDIR = $(FIELDTRIP)/realtime/bin/maci
        CFLAGS += -m32
        CXXFLAGS += -m32
        LDFLAGS += -m32
    endif
    ifeq "$(MACHINE)" "x86_64"
        BINDIR = $(FIELDTRIP)/realtime/bin/maci64
        CFLAGS += -m64
        CXXFLAGS += -m64
        LDFLAGS += -m64
    endif
endif

ifeq "$(UNAME)" "Linux"
    LDLIBS += -lbuffer -lpthread -ldl -lrt
    LDFLAGS += -L$(FTBUFFER)/src
    ifeq "$(MACHINE)" "i686"
        BINDIR = $(FIELDTRIP)/realtime/bin/glnx86
    endif
    ifeq "$(MACHINE)" "x86_64"
        BINDIR = $(FIELDTRIP)/realtime/bin/glnxa64
    endif
endif

ifndef BINDIR
$(warning Unsupported platform: $(PLAT) :/.)
endif

TARGETS = $(BINDIR)/amp2ft
all: $(TARGETS)

AmpServerClient.o: AmpServerClient.h

$(BINDIR)/amp2ft: amp2ft.o AmpServerClient.o
	$(CXX) -o $@ $^ $(LDFLAGS) $(LDLIBS) $(CXXFLAGS)

clean:
	rm -f *.o $(TARGETS)
