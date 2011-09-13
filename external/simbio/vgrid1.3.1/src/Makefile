

AR=ar -r
RM=rm -f
IFLAGS=-I../include -I../vista/include

SRC=generateGrid.C  marching-tets.C  materialDB.C  mesh.C  trilinear.C

OBJ=$(SRC:%.C=%.o)

.SUFFIXES: .C .o
.C.o:
	$(CXX) -c $(IFLAGS) $(CXXFLAGS) $<

all: lib

lib: $(OBJ)
	$(RM) libvgrid.a
	$(AR) libvgrid.a $(OBJ)

clean:
	rm -f *.o *~ *.a