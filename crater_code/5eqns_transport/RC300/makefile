#
# Makefile for multiphase code
#
FFLAGS = -O  
LFLAGS =               
FC     = gfortran     
#
OBJECTS = \
driver.o \
setprob.o \
qinit.o \
../lib/tclaw2ez.o \
../lib/tclaw2.o \
../lib/tclaw2_maloc.o \
../lib/step2.o \
../lib/flux2.o \
../lib/ts2_erk1.o \
../lib/ts2_erk2_ssp22.o \
../lib/ts2_erk3_ssp33.o \
../lib/q2qlqr_interp_tvd2_qp.o \
../lib/q2qlqr_thinc_m.o \
../lib/rpn2.o \
../lib/hllc.o \
../lib/tfluct2.o \
../lib/prmtoc.o \
../lib/ctomfd.o \
../lib/ctoqe.o \
../lib/rhoe_prho_ideal.o \
../lib/rhoe_prho_sg.o \
../lib/tsrc2_gaxis.o \
../lib/setaux2.o \
../lib/setaux2_grav.o \
../lib/mapc2p.o \
../lib/out2.o \
../lib/out1.o \
../lib/bc2.o \
../lib/bc2_lr.o \
../lib/estdt2.o \
../lib/cellave2.o \
../lib/fss.o \
../lib/a4step2.o \
../lib/b4step2.o \
../lib/copyaux2.o \
../lib/copyq2.o \
../lib/restart2.o \
../lib/bc2_aux.o \
../lib/auxbc2.o \
../lib/auxbc2_lr.o \
../lib/dg2.o \
../lib/outgrd.o \
../lib/thinc_LR.o \
../lib/philim.o \
../lib/zero1d_aux.o \
../lib/zeroin.o

SOURCES = \

xclaw: $(OBJECTS)
	$(FC) $(LFLAGS) $(OBJECTS) -o xclaw 

clear: $(OBJECTS)
	rm $(OBJECTS) xclaw

program : $(SOURCES)
	cat $(SOURCES) > claw2program.f

### DO NOT remove this line - make depends on it ###
