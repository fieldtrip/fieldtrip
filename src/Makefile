MEXCMD = mex
MEXEXT = $(shell mexext)

FIELDTRIPSRC = $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
FIELDTRIPROOT = $(shell dirname $(FIELDTRIPSRC))

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

CPPMEXFILES = \
	combineClusters.$(MEXEXT)

SPECIALMEXFILES = \
	mxDeserialize_c.$(MEXEXT) \
	mxDeserialize_cpp.$(MEXEXT) \
	mxSerialize_c.$(MEXEXT) \
	mxSerialize_cpp.$(MEXEXT) \
	plgndr.$(MEXEXT) \
	rfbevent.$(MEXEXT)

ALLMEXFILES = $(MEXFILES) $(GEOMETRYMEXFILES)

all: $(MEXFILES) $(GEOMETRYMEXFILES)

$(MEXFILES): %.$(MEXEXT): %.c
	$(MEXCMD) $<

$(GEOMETRYMEXFILES): %.$(MEXEXT): %.c geometry.c geometry.h
	$(MEXCMD) $< geometry.c 

$(CPPMEXFILES): %.$(MEXEXT): %.cpp
	$(MEXCMD) $<

install:
	@echo ""
	@echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
	@echo "📦 Installing MEX files"
	@echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
	@for file in $(ALLMEXFILES); do \
		if [ -f "$$file" ]; then \
			count=$$(find $(FIELDTRIPROOT) -type f -name "$$file" -not -path "$(FIELDTRIPROOT)/src/*" | wc -l | tr -d ' '); \
			if [ $$count -gt 0 ]; then \
				echo "  📄 $$file → updating $$count location(s)"; \
				find $(FIELDTRIPROOT) -type f -name "$$file" -not -path "$(FIELDTRIPROOT)/src/*" -exec cp -f "$$file" {} \; ; \
			else \
				echo "  📄 $$file → no matching files found in tree"; \
			fi \
		else \
			echo "  ⚠️  Warning: $$file not found, skipping..."; \
		fi \
	done
	@echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
	@echo "✅ Deployment complete"
	@echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"


clean:
	@echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
	@echo "🧹 Removing compiled MEX files here..."
	@echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
	rm -f $(ALLMEXFILES)
	@echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
	@echo "✅ Clean complete"
	@echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
