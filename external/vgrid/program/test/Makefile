
.SUFFIXES: .C .o

.C.o:
	$(CXX) -o $*.o -c $(CXXFLAGS) $(IFLAGS) $<


IFLAGS=-I../../vista/include
LDFLAGS=-L../../vista/src
LIBS=-lvista

all: genimg
	./run_test_suite

genimg.o: genimg.C

genimg: genimg.o
	$(CXX) -o genimg genimg.o $(LDFLAGS) $(LDPATH) $(LIBS)

alltests: simple check4x4x4-min2  check-mt

simple: genimg
	./test2x2x2.sh 
	./test4x4x4.sh


# 2x2x2 image with voxel (0,0,0) different material
check2x2x2: genimg
	./test2x2x2.sh 

# 4x4x4 image with voxel (0,0,0) different material
check4x4x4: genimg
	./test4x4x4.sh 

# 4x4x4 image with voxel (0,0,0) different material
check4x4x4-min2: genimg
	./test4x4x4-min2.sh 

check-mt: check-mt-two-halves check-mt-oblique2d

# marching tets with axis-pararallel material interface
check-mt-two-halves: genimg
	./test-mt-two-halves.sh 

# marching tets with interfaces parallel to plane x+y = const.
check-mt-oblique2d: genimg
	./test-mt-oblique2d.sh 



clean:
	rm -f genimg *.o *~ *.log  *.v *.ascii *.out *.gmv *.complex2d *.cells *.vertices 
	rm -f core core.*
	rm -rf out-mt
