CC=gcc
CFLAGS=-I. -Wall
LFLAGS=-lncurses
DEPS = AcqBuffer.h

all: AcqBuffer acq2ft fake_meg

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

AcqBuffer: AcqBuffer.o
	gcc -o $@ $^ $(CFLAGS) $(LFLAGS)

acq2ft: acq2ft.c 
	gcc -o acq2ft acq2ft.c  -I../buffer/src -Wall ../buffer/src/libbuffer.a  -lpthread
	
fake_meg: fake_meg.c
	gcc -o fake_meg fake_meg.c -I. -Wall
