# clear all the implicit rules
.SUFFIXES :


# determine the correct matlab root
MATLABROOT:=$(wildcard /usr/local/matlab/bin /opt/matlab/bin /local/share/matlab/bin /Applications/MATLAB*/bin)
MATLABROOT:=$(firstword $(MATLABROOT))
MEX:=$(MATLABROOT)/mex

#MEXEXT:=mexa64
MEXEXT:=mexglx  # N.B. override in the call with MEXEXT=`mexext` for macs etc.
# alternative which finds the extension automatically, when mexext exists
# MEXEXT:=$(shell $(MATLABROOT)/mexext)

FLAGS:=-g -ansi -pedantic -Wall
#FLAGS:=-O -g  

OPS=plus minus times rdivide ldivide power eq ne lt gt le ge

.PHONY: all
all: tprod repop 

.PHONY: tprod
tprod : tprod.$(MEXEXT)

.PHONY: repop
repop : repop.$(MEXEXT)

.PHONY: repops
repops : $(OPS:%=r%.$(MEXEXT))


# sompe dependency information
repop.$(MEXEXT) : repop_mex.c repop_util.c ddrepop.c dsrepop.c sdrepop.c ssrepop.c repop.def repop.h mxInfo.c mxInfo_mex.c
tprod.$(MEXEXT) : tprod_mex.c tprod_util.c ddtprod.c dstprod.c sdtprod.c sstprod.c tprod.h tprod.def mxInfo.c mxInfo_mex.c

tprod_testcases : tprod_testcases.c mxInfo.c mxInfo.h mxUtils.c mxUtils.h tprod_util.c ddtprod.c dstprod.c sdtprod.c sstprod.c tprod.h tprod.def
	$(CC) $(filter %.c,$^) $(FLAGS) -o $@

# general rule to make utility object codes
%.o : %.c %.h 
	$(MEX) -c $< $(FLAGS)

# general rule to make mex files
# ordered in decreasing specificity so the first one which matches fires
%.$(MEXEXT) : mxInfo.c %.def %.h
	$(MEX) $(filter %.c,$^) $(FLAGS) -output $* $(filter %.o,$^)

%.$(MEXEXT) : mxInfo.c %.def
	$(MEX) $(filter %.c,$^) $(FLAGS) -output $* $(filter %.o,$^)

%.$(MEXEXT) : %.c mxInfo.c 
	$(MEX) $(filter %.c,$^) $(FLAGS) $(filter %.o,$^)

%.$(MEXEXT) : %.c
	$(MEX) $(filter %.c,$^) $(FLAGS) $(filter %.o,$^)

# repops are special cases like genops
r%.$(MEXEXT) : repop_ind.c repop.c repop.def repop.h mxInfo.c
	$(MEX) $(filter %.c,$^) $(FLAGS) -D_`echo $* | tr a-z A-Z`_ -output r$* $(filter %.o,$^)
# 	$(MEX) $< $(FLAGS) -D_`echo $* | sed -e 's/^.//' | tr a-z A-Z`_ -output $*

clean:
	rm -f *.$(MEXEXT) *.o
