
#gfortran
LDF	= gfortran
LDflags   = -O3 -ffree-form -ffree-line-length-none -cpp -Warray-bounds -Wall -Wunderflow -Warray-temporaries -Wcharacter-truncation -Wtabs -Wintrinsic-shadow -Walign-commons -frange-check -Wconversion -fbounds-check -ffpe-trap=invalid,zero,overflow -fbacktrace -pedantic -fimplicit-none

#ifort
#LDF     = ifort
#LDflags   = -O3 -fpe0 -check all -traceback -warn

LIBS = -llapack -lblas

OBJSF0  = rbm_fam.f90


.SUFFIXES : .f90 .o .c++

FRUN	 = RBM_FAM

all:	$(FRUN)

$(FRUN): $(OBJSF0) $(OBJCPP)
	@echo "$(LDF) $(OBJCPP) $(OBJSF0) $(LDflags) -o $(FRUN) $(LIBS)"
	$(LDF) $(OBJCPP) $(OBJSF0) $(LDflags) -o $(FRUN) $(LIBS)
	@echo "done"

.c++.o:
	@echo compilation of $<
	$(CPCPP) $(CPCPPflags) $< -c -o $@

clean:
	rm *.o *.mod RBM_FAM
