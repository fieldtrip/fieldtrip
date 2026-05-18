ifeq (, $(shell which matlab))
$(error "The MATLAB executable cannot be found, please add it to your PATH.")
endif

MATLAB = $(shell which matlab)
MATLABBIN = $(dir $(realpath $(MATLAB)))
MEX = $(MATLABBIN)mex
MEXEXT = $(shell $(MATLABBIN)mexext)

$(info MATLAB = $(MATLAB))
$(info MEX    = $(MEX))
$(info MEXEXT = $(MEXEXT))

FIELDTRIPSRC = $(dir $(realpath $(lastword $(MAKEFILE_LIST))))
FIELDTRIPROOT = $(dir $(realpath $(FIELDTRIPSRC)))

MEXFILES = \
	det2x2.$(MEXEXT) \
	det3x3.$(MEXEXT) \
	ft_getopt.$(MEXEXT) \
	ft_spike_sub_crossx.$(MEXEXT) \
	getpid.$(MEXEXT) \
	inv2x2.$(MEXEXT) \
	inv3x3.$(MEXEXT) \
	meg_leadfield1.$(MEXEXT) \
	mtimes2x2.$(MEXEXT) \
	mtimes3x3.$(MEXEXT) \
	nanmean.$(MEXEXT) \
	nanstd.$(MEXEXT) \
	nansum.$(MEXEXT) \
	nanvar.$(MEXEXT) \
	plgndr.$(MEXEXT) \
	read_16bit.$(MEXEXT) \
	read_24bit.$(MEXEXT) \
	rename.$(MEXEXT) \
	sandwich2x2.$(MEXEXT) \
	sandwich3x3.$(MEXEXT) \
	splint_gh.$(MEXEXT)

GEOMETRYMEXFILES = \
	lmoutr.$(MEXEXT) \
	ltrisect.$(MEXEXT) \
	plinproj.$(MEXEXT) \
	ptriproj.$(MEXEXT) \
	routlm.$(MEXEXT) \
	solid_angle.$(MEXEXT)

CTFMEXFILES = \
	read_ctf_shm.$(MEXEXT) \
	write_ctf_shm.$(MEXEXT)

D3DESMEXFILES = \
	rfbevent.$(MEXEXT)

CPPMEXFILES = \
	combineClusters.$(MEXEXT)

OBSOLETEMEXFILES = \
	mxDeserialize_cpp.$(MEXEXT) \
	mxSerialize_cpp.$(MEXEXT) \
	mxDeserialize.$(MEXEXT) \
	mxSerialize.$(MEXEXT)

ALLMEXFILES = $(MEXFILES) $(GEOMETRYMEXFILES) $(CPPMEXFILES)

ifeq ($(MEXEXT),mexa64)
	ALLMEXFILES += $(D3DESMEXFILES)
	ALLMEXFILES += $(CTFMEXFILES)
else ifeq ($(MEXEXT), mexmaca64)
	ALLMEXFILES += $(D3DESMEXFILES)
endif

all: $(ALLMEXFILES)

$(MEXFILES) $(CTFMEXFILES): %.$(MEXEXT): %.c
	$(MEX) $<

$(CPPMEXFILES): %.$(MEXEXT): %.cpp
	$(MEX) $<

$(GEOMETRYMEXFILES): %.$(MEXEXT): %.c geometry.c geometry.h
	$(MEX) $< geometry.c 

$(D3DESMEXFILES): %.$(MEXEXT): %.c d3des.c d3des.h
	$(MEX) $< d3des.c

install: $(ALLMEXFILES)
	@echo ""
	@echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
	@echo "📦 Installing compiled MEX files ..."
	@echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
	@for file in $^; do \
		TARGETS=$$(find $(FIELDTRIPROOT) -type f -name "$$file" -not -path "$(FIELDTRIPSRC)*" ); \
		if [ -z "$$TARGETS" ]; then \
			echo "⚠️ No target found for $$file, skipping..."; \
			continue; \
		fi; \
		for target in $$TARGETS; do \
			cp -f "$(FIELDTRIPSRC)/$$file" "$$target" && \
			echo "✅ Installed $$file to $$target" || \
			echo "❌ Failed to install $$file to $$target"; \
		done; \
	done
	@echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
	@echo "✔️ Installation complete"
	@echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

clean:
	@echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
	@echo "🧹 Removing compiled MEX files ..."
	@echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
	rm -f $(ALLMEXFILES)
	@echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
	@echo "✅ Clean complete"
	@echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
