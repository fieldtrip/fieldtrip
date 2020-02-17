#!/usr/bin/env make -f
#
# Work in Progress: Makefile for FieldTrip
#
# Currently only supports generating .mex files,
# and removing existing ones.
#
# Usage:
#
# 	make mex MATLAB=`which matlab`
#   make mex OCTAVE=`which octave`
# 	make clean MATLAB=`which matlab`
#   make clean OCTAVE=`which octave`
#
# "make mex" builds the .mex files in src, and then copies the resulting
# mex files to other directories in FieldTrip.
# "make clean" deletes .mex files in src and other directories in FieldTrip
#
# Depending on the Matlab or Octave binary provided, this Makefile
# selects the appropriate extension for .mex files.
#
# Note that Octave uses a single extension (".mex") for different platforms
# (e.g. Linux, OSX, Windows; 32 and 64 bit) so it seems a Bad Idea
# to include compiled .mex files for Octave. 
# As of Feb 2016, FieldTrip comes with pre-compiled .mex* files for 
# a variety of platforms
#
# Exclusions: 
# 	src/mxDeserialize*
#   src/mxSerialize*
#   src/rfbevent*
#   compat/matlablt2010b/@uint64/*
#
#
# 2017 Nikolaas N. Oosterhof

help:
	@echo "Build or remove FieldTrip .mex* files"
	@echo ""
	@echo "Usage: make <pathspec> <target>"
	@echo ""
	@echo "where <pathspec> can be:"
	@echo "------------------------------------------------------------------"
	@echo "OCTAVE=\`which octave\`  # or any other path to Octave"
	@echo "MALTAB=\`which matlab\`  # or any other path to Matlab"
	@echo "------------------------------------------------------------------"
	@echo "and  <target> is one of:"
	@echo "------------------------------------------------------------------"
	@echo "  mex                compile FieldTrip mex* files"
	@echo "  clean              remove FieldTrip mex* files"
	

# FieldTrip source directory
SRC=$(shell pwd)/src

# Some helpers for running mex commands
MEXRUNPRE=succes=false;try,succes=~mex(
MEXRUNPOST=);catch,end;exit(~succes);

MEXPREFIX="[pth,nm,ext]=fileparts('
MEXSUFFIX=');cd(pth);$(MEXRUNPRE)[nm,ext]$(MEXRUNPOST)"

# determine mex extension from Matlab's or Octave's "mexext" function
octave_or_matlab:
	$(if $(MATLAB),,$(if $(OCTAVE),,$(error Set either the \
			MATLAB and OCTAVE variable to the \
			appropriate binary, for example \
			OCTAVE=`which octave` or MATLAB=`which matlab`)))

ifdef MATLAB
  BIN=$(MATLAB)
  FLAGS:=-nojvm -nodisplay -nosplash  -r
else ifdef OCTAVE
  BIN=$(OCTAVE)
  FLAGS=--eval
endif

ifdef BIN
  EXT:=$(shell $(BIN) $(FLAGS) "disp(mexext());exit" | tail -2 | head -1)
endif

.PHONY: all

# set which files are not included 
EXCLUDE_FIND_PAT= -not -iname 'mxDeserialize*' -not -iname 'mxSerialize*' -not -iname 'rfbevent*'

# The command below was used to generate
# the list of files to be compiled. It can be 
# copy/pasted in the terminal:
# EXCLUDE_FIND_PAT="! -iname 'mxDeserialize*' ! -iname 'mxSerialize*' ! -iname 'rfbevent*'"
#  
#
# CORE="./@config ./bin ./compat ./connectivity  \
  ./engine ./fileio ./forward ./inverse ./multivariate \
  ./plotting ./preproc ./private ./specest \
  ./statfun ./trialfun ./utilities"; \
  for i in `find src -iname "*.mexmaci64" ${EXCLUDE_FIND_PAT} `; do \
	j=`basename $i`; \
    echo $i $j; continue; \
    for c in ${CORE}; do \
            find ${c} -iname ${j} | sed s/\^\.\\///g; \
    done | sed s/\.mexmaci64/\ \\\\/g; done

# targets with .mex* extension outside $SRC
OTHER_TARGETS_PF=connectivity/private/det2x2 \
		connectivity/private/ft_getopt \
		connectivity/private/inv2x2 \
		connectivity/private/mtimes2x2 \
		connectivity/private/sandwich2x2 \
		engine/private/ft_getopt \
		fileio/private/ft_getopt \
		fileio/private/read_16bit \
		fileio/private/read_24bit \
		forward/private/ft_getopt \
		forward/private/lmoutr \
		forward/private/meg_leadfield1 \
		forward/private/plgndr \
		forward/private/ptriproj \
		forward/private/routlm \
		forward/private/solid_angle \
		inverse/private/ft_getopt \
		inverse/private/solid_angle \
		plotting/private/ft_getopt \
		plotting/private/ltrisect \
		private/det2x2 \
		private/inv2x2 \
		private/lmoutr \
		private/ltrisect \
		private/mtimes2x2 \
		private/plgndr \
		private/ptriproj \
		private/routlm \
		private/sandwich2x2 \
		private/solid_angle \
		specest/private/ft_getopt \
		utilities/ft_getopt \
		utilities/private/lmoutr \
		utilities/private/ptriproj \
    external/stats/nanmean \
    external/stats/nanstd \
    external/stats/nansum \
    external/stats/nanvar	

# targets with .mex* extension outside $SRC
OTHER_TARGETS=$(patsubst %,%.$(EXT),$(OTHER_TARGETS_PF))

# targets with .mex* extension inside $SRC
SRC_TARGETS=$(shell find $(SRC) -iname "*.$(EXT)" $(EXCLUDE_FIND_PAT))

comma:=,

# The files here depend on $(SRC)/geometry.c
$(SRC)/ptriproj.$(EXT) \
    $(SRC)/lmoutr.$(EXT) \
	$(SRC)/routlm.$(EXT) \
	$(SRC)/solid_angle.$(EXT) \
	$(SRC)/ltrisect.$(EXT): 			$(SRC)/geometry.c
	$(BIN) $(FLAGS) "cd('$(@D)');$(MEXRUNPRE)'$(basename $(@F)).c'$(addprefix $(comma)',$(addsuffix .c',$(basename $(?F))))$(MEXRUNPOST)"

# All other files do not have dependencies
$(SRC)/%.$(EXT): $(SRC)/%.c
	$(BIN) $(FLAGS) $(MEXPREFIX)$(<)$(MEXSUFFIX)

# Helper function to define how to get .mex* files not in $(SRC).
define addrule
$(1) : src/$(notdir $(1))
	cp src/$(notdir $(1)) $(1)
endef

# Apply helper function to each element in $(OTHER_TARGETS)
$(foreach a, $(OTHER_TARGETS),$(eval $(call addrule, $(a))))

# Helper rule, see https://stackoverflow.com/questions/16467718/how-to-print-out-a-variable-in-makefile
print-%: ; @echo $* = $($*)

#######################################
#
# Main rules

all: mex

clean: mex_clean

# For mex files
mex: octave_or_matlab $(SRC_TARGETS) $(OTHER_TARGETS)

mex_clean: octave_or_matlab
	rm -f $(SRC_TARGETS) $(OTHER_TARGETS) 
