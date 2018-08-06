# Makefile for compiling and linking the proxy for sending data
# from Elekta Neuromag MEG system to the FieldTrip real-time buffer
#
# (C) Lauri Parkkonen and Gustavo Sudre, 2010

ifneq "$(OS)" "Windows_NT"
	OS      ?= $(shell uname -s)
	MACHINE ?= $(shell uname -m)
endif

PLATFORM ?= $(shell gcc -dumpmachine)

ifeq ($(OS), Linux)
	# The Elekta Linux acquisition system is always a 64-bit system
	ARCH            = x86_64-pc-linux-gnu
	CFLAGS_PLATFORM = -I ./include/dacq_rel5up
	LIBS_PLATFORM   = -ldacq_rel5up
else ifeq ($(OS), HP-UX)
	# The Elekta HP-UX acquisition system shall be HP-UX 11.11
	ARCH            = hppa2.0-hp-hpux11.11
	CFLAGS_PLATFORM = -D_HPUX_SOURCE -DDACQ_OLD_CONNECTION_SCHEME -I ./include/dacq_rel5up
	LIBS_PLATFORM   = -ldacq_rel4
else
	$(error Unsupported platform. Compilation and linking only on Linux or HP-UX)
endif

CC       = gcc
PROGRAM  = bin/$(ARCH)/neuromag2ft
OBJECTS  = neuromag2ft.o process_tag.o process_data.o
HEADER   = neuromag2ft.h

CFLAGS   = -Wall -I ../../buffer/src -I ./include $(CFLAGS_PLATFORM)
LIBPATH  = -L ../../buffer/src -L ./lib/$(ARCH)
USERLIBS = -lbuffer $(LIBS_PLATFORM) -ldacqcomm -lfiff -lerr-2.1
SYSLIBS  = -lpthread -lm -lrt

##############################################################################
all : $(OBJECTS) ../../buffer/src/libbuffer.a
	$(CC) $(CFLAGS) $(LIBPATH) $(OBJECTS) $(USERLIBS) $(SYSLIBS) -o $(PROGRAM)

*.o: *.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $^

clean:
	rm -f $(OBJECTS) $(PROGRAM)
