#!/usr/bin/env make -f
# GIfTI Makefile
#
# Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging
#
# $Id: Makefile 6402 2015-04-09 18:11:06Z guillaume $

include Makefile.var

SPMMEX = zstream.$(SUF)

ifeq (mex,$(SUF))
	export CFLAGS = $(shell $(MEX) -p CFLAGS) -std=c99
else
	MEXOPTS += CFLAGS='$$CFLAGS -std=c99'
endif

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
