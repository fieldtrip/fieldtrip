# Makefile for compiling and linking the proxy for sending data
# from Elekta Neuromag MEG system to the FieldTrip real-time buffer
#
# (C) Lauri Parkkonen and Gustavo Sudre, 2010

PLATFORM = x86_64-pc-linux-gnu     # Use with -ldacq_rel5up in USERLIBS
#PLATFORM = hppa2.0-hp-hpux11.11  # Use with -ldacq_rel4 in USERLIBS and -D_HPUX_SOURCE and -DDACQ_OLD_CONNECTION_SCHEME in CFLAGS

CC       = gcc
PROGRAM  = neuromag2ft
OBJECTS  = neuromag2ft.o process_tag.o process_data.o
HEADER   = neuromag2ft.h

CFLAGS   = -g -Wall -I ../../buffer/src -I ./include
#CFLAGS   = -g -Wall -I ../../buffer/src -I ./include -D_HPUX_SOURCE -DDACQ_OLD_CONNECTION_SCHEME
LDFLAGS  = -L ../../buffer/src -L ./lib/$(PLATFORM)
USERLIBS = -lbuffer -ldacq_rel5up -ldacqcomm -lfiff -lerr-2.1
#USERLIBS = -lbuffer -ldacq_rel4 -ldacqcomm -lfiff -lerr-2.1
SYSLIBS  = -lpthread -lm

all : $(OBJECTS) ../../buffer/src/libbuffer.a
	$(CC) $(CFLAGS) $(LDFLAGS) $(OBJECTS) $(USERLIBS) $(SYSLIBS) -o $(PROGRAM)

*.o : *.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $^ 

clean :
	rm -f $(OBJECTS) $(PROGRAM)
