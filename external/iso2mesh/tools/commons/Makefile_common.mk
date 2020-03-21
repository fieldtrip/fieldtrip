########################################################
#
#  Makefile for Compiling External Tools needed by Iso2Mesh
#  Copyright (C) 2018 Qianqian Fang <q.fang at neu.edu>
#
########################################################

########################################################
# Base Makefile for all example/tests and main program
#
# This file specifies the compiler options for compiling
# and linking
########################################################

ifndef ROOTDIR
ROOTDIR := .
endif

ifndef I2MDIR
I2MDIR := $(ROOTDIR)
endif

BXDSRC :=$(I2MDIR)/src

CXX        := g++
CC         := gcc
AR         := g++
BIN        := ../bin
BUILT      := built
BINDIR     := $(BIN)
OBJDIR 	   := $(BUILT)
CCFLAGS    := -c -Wall -O3  #-g
INCLUDEDIR := $(I2MDIR)/src
ARFLAGS    :=
AROUTPUT   := -o
MAKE       :=make
CMAKE      :=cmake

COPY	 := cp
ECHO	 := echo
MKDIR    := mkdir

OBJSUFFIX        := .o
BINSUFFIX        := 

OBJS      := $(addprefix $(OBJDIR)/, $(FILES))
OBJS      := $(subst $(OBJDIR)/$(BXDSRC)/,$(BXDSRC)/,$(OBJS))
OBJS      := $(addsuffix $(OBJSUFFIX), $(OBJS))

TARGETSUFFIX:=$(suffix $(BINARY))

ifeq ($(TARGETSUFFIX),.so)
	CCFLAGS+= -fPIC 
	ARFLAGS+= -shared -Wl,-soname,$(BINARY).1 
endif

ifeq ($(TARGETSUFFIX),.a)
        CCFLAGS+=
	AR         := ar
        ARFLAGS    := r
	AROUTPUT   :=
endif

all: $(SUBDIRS) makedirs copybin

$(SUBDIRS):
	if test -f $@/CMakeLists.txt; then cd $@ && $(CMAKE) . && cd ..; fi
	$(MAKE) -C $@ --no-print-directory

copybin:
	@$(COPY) cgalmesh/mesh_3D_image           ../bin/cgalmesh
	@$(COPY) cgalmesh/mesh_polyhedral_domain  ../bin/cgalpoly
	@$(COPY) cgalsurf/mesh_a_3d_gray_image    ../bin/cgalsurf
	@$(COPY) cgalsimp2/edge_collapse_enriched_polyhedron  ../bin/cgalsimp2
	@$(COPY) cork/bin/cork  ../bin/cork
	@$(COPY) meshfix/meshfix  ../bin/meshfix
	@$(COPY) meshfix/contrib/JMeshLib/test/jmeshlib  ../bin/jmeshlib
	@$(COPY) tetgen/tetgen  ../bin/tetgen1.5

makedirs:
	@if test ! -d $(BINDIR); then $(MKDIR) $(BINDIR); fi
#	@if test ! -d $(OBJDIR); then $(MKDIR) $(OBJDIR); fi

.SUFFIXES : $(OBJSUFFIX) .cpp

##  Compile .cpp files ##
$(OBJDIR)/%$(OBJSUFFIX): %.cpp
	@$(ECHO) Building $@
	$(CXX) $(CCFLAGS) $(USERCCFLAGS) -I$(INCLUDEDIR) -o $@  $<

##  Compile .cpp files ##
%$(OBJSUFFIX): %.cpp
	@$(ECHO) Building $@
	$(CXX) $(CCFLAGS) $(USERCCFLAGS) -I$(INCLUDEDIR) -o $@  $<

##  Compile .c files  ##
$(OBJDIR)/%$(OBJSUFFIX): %.c
	@$(ECHO) Building $@
	$(CC) $(CCFLAGS) $(USERCCFLAGS) -I$(INCLUDEDIR) -o $@  $<

##  Link  ##
$(BINDIR)/$(BINARY): $(OBJS)
	@$(ECHO) Building $@
	$(AR)  $(ARFLAGS) $(AROUTPUT) $@ $(OBJS) $(USERARFLAGS)

## Clean
clean:
	rm -rf $(OBJS) $(OBJDIR) #$(BINDIR)
ifdef SUBDIRS
	for i in $(SUBDIRS); do $(MAKE) --no-print-directory -C $$i clean; done
endif

.PHONY: regression clean arch makedirs dep $(SUBDIRS)

