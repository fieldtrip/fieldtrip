
SUBDIRS= vista src program

CFLAGS   = -O2 -Wall # C flags
CXXFLAGS = -O2 -Wall # C++ flags
CC = cc   # C compiler
CXX= c++  # C++ compiler

all:
	@set -e; \
        for i in $(SUBDIRS); \
        do cd $$i; $(MAKE) all;\
        cd ..; done

clean:
	@set -e; \
        for i in $(SUBDIRS); \
        do cd $$i; $(MAKE) clean;\
        cd ..; done


