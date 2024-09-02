FC=gfortran
FFLAGS=
FFLAGSDEBUG=
FDEBUG=-g -Wall -Wextra -fcheck=all -fbacktrace
FOPTIM=-O3


MAINPRG=Program
MAINSRC=$(MAINPRG).f90

PRECPRG=PrecisionMod
PRECSRC=$(PRECPRG).f90

TYPEPRG=TypeCustom
TYPESRC=$(TYPEPRG).f90

SUBPRG=SubRoutines
SUBSRC=$(SUBPRG).f90

FFLAGSDEBUG+=$(FDEBUG)
FFLAGS+=$(FOPTIM)

LIB=usr/lib/x86_64-linux-gnu/hdf5/serial
LIB2=usr/local/lib
INC1=usr/include/hdf5/serial
INC2=usr/local/include
INCLOCAL=./



build: $(OBJ)
	$(FC) $(FFLAGS) -I/$(INC1) -I/$(INC2) -L/$(LIB) -L/$(LIB2) -L/$(LIB3) -c -o $(PRECPRG).o $(PRECSRC) -lhdf5 -lhdf5_fortran -fopenmp
	$(FC) $(FFLAGS) -I/$(INC1) -I/$(INC2) -I/$(INCLOCAL) -L/$(LIB) -L/$(LIB2) -L/$(LIB3) -o $(TYPEPRG).o $(TYPESRC) $(PRECPRG).o -lhdf5 -lhdf5_fortran -fopenmp
#	$(FC) $(FFLAGS) -I/$(INC1) -I/$(INC2) -L/$(LIB) -L/$(LIB2) -L/$(LIB3) -c $(SUBSRC) $(PRECPRG).o $(TYPEPRG).o -lhdf5 -lhdf5_fortran -fopenmp
#	$(FC) $(FFLAGS) -I/$(INC1) -I/$(INC2) -L/$(LIB) -L/$(LIB2) -L/$(LIB3) -o $(MAINPRG) $(MAINSRC) $(SUBPRG).o $(PRECPRG).o $(TYPEPRG).o -lhdf5 -lhdf5_fortran -fopenmp

#optim: $(OBJ)


#debug: $(OBJ)
#	$(FC) $(FFLAGSDEBUG) -I/$(INC1) -I/$(INC2) -L/$(LIB) -L/$(LIB2) -o $(APP) $(SRC) -lhdf5 -lhdf5_fortran -lfftw3 -fopenmp -llapack -lblas

run:
	make build; ./$(MAINPRG)

debugrun:
	make debug; ./$(MAINPRG)

cleanrun:
	make cleanall; make run

clean:
	rm -rf $(OBJ)
	rm -rf *.mod
	rm -rf $(MAINPRG).x
	rm -rf *.o
	rm -rf ./$(MAINPRG)
