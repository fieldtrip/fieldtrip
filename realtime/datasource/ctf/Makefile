CC=gcc
CFLAGS=-I. -Wall
LFLAGS=-lncurses
DEPS = AcqBuffer.h

all: AcqBuffer acq2ft acq2ftx fake_meg check_shm acq2ftx_fake

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

AcqBuffer: AcqBuffer.o
	gcc -o $@ $^ $(CFLAGS) $(LFLAGS)

acq2ft: acq2ft.c 
	gcc -o acq2ft acq2ft.c  -I../../buffer/src -Wall ../../buffer/src/libbuffer.a  -lpthread

acq2ftx: acq2ftx.c 
	gcc -o acq2ftx acq2ftx.c  -I../../buffer/src -Wall ../../buffer/src/libbuffer.a  -lpthread -lm
	
acq2ftx_fake: acq2ftx.c 
	gcc -o acq2ftx_fake acq2ftx.c  -DFAKE_MEG -I../../buffer/src -Wall ../../buffer/src/libbuffer.a  -lpthread -lm
	
fake_meg: fake_meg.c
	gcc -o fake_meg fake_meg.c -I. -Wall

check_shm: check_shm.c
	gcc -o check_shm check_shm.c -I. -Wall
