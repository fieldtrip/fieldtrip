FFLAGS=$(F77FLAGS)
CFLAGS=$(CCFLAGS)
AR=$(AR_PROG)
FC=$(FC_PROG)
CC=$(CC_PROG)
RANLIB=$(RANLIB_PROG)
MAKE=$(MAKE_PROG)

LIB    =$(LIBDIR)/libsolver$(SUFFIX).a
CLEAN  =$(LIB)

OBJ=	$(LIB)(axunsu.o) $(LIB)(axunsy.o) $(LIB)(cgstab.o) $(LIB)(esv.o) \
	$(LIB)(freund.o) $(LIB)(illcon.o) $(LIB)(kograd.o) $(LIB)(matpro.o) \
	$(LIB)(partch.o) $(LIB)(pccgsy.o) $(LIB)(prcgst.o) $(LIB)(precnd.o) \
	$(LIB)(ptfqmr.o) $(LIB)(ransol.o) $(LIB)(scalen.o) $(LIB)(scalun.o) \
	$(LIB)(solchs.o) $(LIB)(solvel.o) $(LIB)(solveu.o) $(LIB)(solvlu.o) \
	$(LIB)(unscal.o) $(LIB)(unscun.o)  

.f.o:
	$(FC) -c $(FFLAGS) $<

.f.a:
	$(FC) -c $(FFLAGS) $<
	$(AR) -rv $(LIB) $*.o 
	-rm -f $*.o $*.map

.c.a:
	$(CC) -c $(CFLAGS) $<
	$(AR) -rv $(LIB) $*.o 
	-rm -f $*.o $*.map

all: $(LIB) #$(LIBDIR)/bdata/$(BDATA)$(SUFFIX).o

$(LIB): $(OBJ) $(INC)
	$(RANLIB) $(LIB)
	-for n in *.map "XX";do if test -w $$n; then rm $$n;fi;done	

#$(LIBDIR)/bdata/$(BDATA)$(SUFFIX).o: $(BDATA).f $(INC) 
#	$(FC) -c $(FFLAGS) $(BDATA).f 
#	mv $(BDATA).o $(LIBDIR)/bdata/$(BDATA)$(SUFFIX).o

clean:  
	-for n in $(CLEAN) "XX";do if test -w $$n; then rm $$n;fi;done

$(LIB)(axunsu.o): axunsu.f
$(LIB)(axunsy.o): axunsy.f
$(LIB)(cgstab.o): cgstab.f
$(LIB)(esv.o): esv.f
$(LIB)(freund.o): freund.f
$(LIB)(illcon.o): illcon.f
$(LIB)(kograd.o): kograd.f
$(LIB)(matpro.o): matpro.f
$(LIB)(partch.o): partch.f
$(LIB)(pccgsy.o): pccgsy.f
$(LIB)(prcgst.o): prcgst.f
$(LIB)(precnd.o): precnd.f
$(LIB)(ptfqmr.o): ptfqmr.f
$(LIB)(ransol.o): ransol.f
$(LIB)(scalen.o): scalen.f
$(LIB)(scalun.o): scalun.f
$(LIB)(solchs.o): solchs.f
$(LIB)(solvel.o): solvel.f
$(LIB)(solveu.o): solveu.f
$(LIB)(solvlu.o): solvlu.f
$(LIB)(unscal.o): unscal.f
$(LIB)(unscun.o): unscun.f

