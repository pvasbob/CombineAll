## in .profile
#PGI=/opt/pgi; export PGI;
#PATH=/opt/pgi/osx86-64/10.1/bin:$PATH; export PATH;
#MANPATH=$MANPATH:/opt/pgi/osx86-64/10.1/man; export MANPATH;
#LM_LICENSE_FILE=$LM_LICENSE_FILE:/opt/pgi/license.dat; export LM_LICENSE_FILE; 

## in a case
#PGI:=/opt/pgi
#export PGI
#PATH:=$(PATH):/opt/pgi/osx86-64/10.1/bin:/opt/local/bin:/opt/local/sbin
#export PATH
#MANPATH:=/opt/pgi/osx86-64/10.1/man
#export MANPATH
#LM_LICENSE_FILE:=/opt/pgi/license.dat
#export LM_LICENSE_FILE
	
#gfortran
CPF	= f95
LDF	= f95

#PGI
#CPF	= /opt/pgi/osx86-64/10.1/bin/pgfortran
#LDF	= /opt/pgi/osx86-64/10.1/bin/pgfortran  

## for yrast
#LDflags  = -O4 -traceback -C -CB -CA -lmkl_ia32 -lguide -lpthread  
#LDflags  = -parallel -openmp -ipo -O4 
#LDflags  = -O4 -traceback -lmkl_ia32 -lguide -lpthread  -Wunused -g -F
#LDflags   = -O0 -ffree-line-length-none -Dhide_dme -Dhide_tho -fbounds-check
LDflags   = -O3 -ffree-line-length-none -Dhide_dme -Dhide_tho -Dhide_mpi_qrpa -Dhide_array_qrpa
 

# for apple vectorized lapack and blas lib
#LIBS = /System/Library/Frameworks/veclib.framework/vecLib
#LIBS = /Users/stoitsov/Downloads/lapack-3.3.0/blas_pghpf_SUNMP.a /Users/stoitsov/Downloads/lapack-3.3.0/lapack_pghpf_SUNMP.a
#LIBS = -llapack -lblas -L . -lgfortran -L /usr/lib/gcc/i386-redhat-linux/4.1.1/
LIBS = -llapack -lblas 

OBJSF0  = HFBTHO137noQRPA.F90 bessik.f90 FAMQRPA18.F90

.SUFFIXES : .f90 .o .c++

FRUN	 = fam.out


all:	$(FRUN)

$(FRUN): $(OBJSF0) $(OBJCPP)
	@echo "$(LDF) $(OBJCPP) $(OBJSF0) $(LDflags) -o $(FRUN) $(LIBS)"
	$(LDF) $(OBJCPP) $(OBJSF0) $(LDflags) -o $(FRUN) $(LIBS)
	@echo "done"

.c++.o:
	@echo compilation of $<
	$(CPCPP) $(CPCPPflags) $< -c -o $@

clean:
	rm ../RUN/*.o ../RUN/*.mod ../RUN/*.out  ../RUN/*.f
