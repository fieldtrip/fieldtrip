# clear all the implicit rules
.SUFFIXES :


# determine the correct matlab root
MATLABROOT:=$(wildcard /usr/local/matlab/bin /opt/matlab/bin /local/share/matlab/bin)
MATLABROOT:=$(firstword $(MATLABROOT))
MEX:=$(MATLABROOT)/mex

#MEXEXT:=mexa64
MEXEXT:=mexglx  # N.B. override in the call with MEXEXT=`mexext` for macs etc.
# alternative which finds the extension automatically, when mexext exists
# MEXEXT:=$(shell $(MATLABROOT)/mexext)

FLAGS:=-g3 
#FLAGS:=-O -g  

OPS=plus minus times rdivide ldivide power eq ne lt gt le ge

.PHONY: all
all: repop 

.PHONY: repop
repop : repop.$(MEXEXT)

.PHONY: repops
repops : $(OPS:%=r%.$(MEXEXT))


# sompe dependency information
repop.$(MEXEXT) : repop_mex.c repop_util.c ddrepop.c dsrepop.c sdrepop.c ssrepop.c repop.def repop.h mxInfo.c mxInfo_mex.c

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
