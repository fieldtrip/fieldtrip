

AR=ar -r
RM=rm -f
IFLAGS=-I../include 

SRC=attr.c  basic.c  binarize.c  convert.c  file.c  graph.c  image.c  option.c
OBJ=$(SRC:%.c=%.o)

.SUFFIXES: .c .o
.c.o:
	$(CC) -c $(IFLAGS) $(CFLAGS) $<

all: lib

lib: $(OBJ)
	echo $(SRC)
	$(RM) libvista.a
	$(AR) libvista.a $(OBJ)

clean:
	rm -f *.o *~ *.a