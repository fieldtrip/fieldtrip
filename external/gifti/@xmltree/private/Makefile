#!/usr/bin/env make -f
# XMLTree Makefile called by {SPM}/src/Makefile
#
# Copyright (C) 2008-2013 Wellcome Trust Centre for Neuroimaging
#
# $Id$

include ../../src/Makefile.var

SPMMEX = xml_findstr.$(SUF)

all: $(SPMMEX)
	@:

clean:
	@:
        
distclean: clean
	$(DEL) $(SPMMEX)

install:
	@:

tarball: all
	$(TAR) cf spm_mex.tar $(SPMMEX)

%.$(SUF) : %.c
	$(MEX) $< $(MEXEND)
