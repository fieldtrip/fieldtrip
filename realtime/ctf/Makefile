CC=gcc
CFLAGS=-I. -Wall
LFLAGS=-lncurses
DEPS = AcqBuffer.h

all: AcqBuffer

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

AcqBuffer: AcqBuffer.o
	gcc -o $@ $^ $(CFLAGS) $(LFLAGS)

