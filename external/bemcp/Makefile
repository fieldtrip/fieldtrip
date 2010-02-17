#!/usr/bin/env make -f
# EEG-BEM Makefile called by {SPM}/external/Makefile
#
# Christophe Phillips
# $Id: Makefile 2990 2009-03-28 17:33:52Z guillaume $

include ../../src/Makefile.var

BEMCPMEX = \
	bem_Cii_cog.$(SUF) bem_Cii_cst.$(SUF) bem_Cii_lin.$(SUF) \
	bem_Cij_cog.$(SUF) bem_Cij_cst.$(SUF) bem_Cij_lin.$(SUF) \
	bem_Gi_cog.$(SUF) bem_Gi_vert.$(SUF)

all: $(BEMCPMEX)

clean:

distclean: clean
	$(DEL) $(BEMCPMEX)

install:

tarball: all
	$(TAR) -cf spm_mex.tar --transform='s,^,bemcp/,g' $(BEMCPMEX)

%.$(SUF) : %.c
	$(MEX) $< $(MEXEND)
