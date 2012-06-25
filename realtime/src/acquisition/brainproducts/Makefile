# Copyright (C) 2010 Stefan Klanke
# Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
# Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
# Kapittelweg 29, 6525 EN Nijmegen, The Netherlands

# Use MinGW for compilation on Windows!

PLATFORM = $(shell gcc -dumpmachine)
FIELDTRIP = ../../../..
FTBUFFER = $(FIELDTRIP)/realtime/src/buffer

#defaults, might be over-written further down
CFLAGS = -Wunused -Wall -pedantic
CFLAGS += -I$(FTBUFFER)/src -I. 
LDFLAGS = -L$(FTBUFFER)/src
LDLIBS = -lbuffer -lpthread
TARGET = $(FTPATH)/realtime/bin/glnxa64/rda2ft


ifeq "$(PLATFORM)" "mingw32"
	fixpath = $(subst /,\,$1)
	BINDIR = $(FIELDTRIP)/realtime/bin/win32
	CFLAGS += -I$(FTBUFFER)/pthreads-win32/include
	LDFLAGS += -LC:/msys/1.0/local/lib
	# LDFLAGS += -L$(BINDIR) does not seem to work for DLLs. Hence:
	LDLIBS += -lws2_32 -lwinmm $(BINDIR)/pthreadGC2.dll
	EXTENSION = .exe
	RM = del
	CP = copy /Y
endif

ifeq "$(PLATFORM)" "i686-apple-darwin10"
 TARGET = $(FTPATH)/realtime/bin/maca64/rda2ft
endif

ifndef BINDIR
$(warning Unsupported plaform: $(PLATFORM) :/.)
endif

TARGETS = $(BINDIR)/rda2ft$(EXTENSION)

all: $(TARGETS)
	
$(BINDIR)/rda2ft$(EXTENSION): rda2ft.c
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(LDLIBS)

clean:
	$(RM) $(call fixpath, *.o $(TARGETS))
