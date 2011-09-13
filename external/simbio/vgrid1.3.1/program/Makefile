

AR=ar -r
RM=rm -f
IFLAGS=-I../include -I../vista/include
LDFLAGS=-L../src -L../vista/src
LIBS=-lvgrid -lvista

SRC=vgrid.C
OBJ=$(SRC:%.C=%.o)

.SUFFIXES: .C .o
.C.o:
	$(CXX) -c $(IFLAGS) $(CXXFLAGS) $<

all: vgrid

vgrid: $(OBJ)
	$(CXX) -o vgrid vgrid.o $(LDFLAGS) $(LIBS)

clean:
	rm -f *.o *~ *.a vgrid