OBJECTS = sparse_matrices.o system_generation.o AIRMs.o phantomgallery.o
BLAS_PATH = #/Users/bavli/lib
LFLAGS =# -L$(BLAS_PATH)
LIB = #-llapack -lblas


FFLAGS = -g -C -O3 -fbackslash -fopenmp

.PHONY: clean help

main.exe: main.f90 $(OBJECTS)
#	make modules
	gfortran $(OBJECTS) main.f90 $(FFLAGS) $(LFLAGS) $(LIB) -o main.exe

run: main.exe
	./main.exe -t -v

modules: sparse_matrices.f90 phantomgallery.f90 system_generation.f90 AIRMS.f90
	gfortran -c sparse_matrices.f90
	gfortran -c system_generation.f90
	gfortran -c phantomgallery.f90
	gfortran -c AIRMS.f90

%.o : %.f90
	gfortran $(FFLAGS) $(LFLAGS) $(LIB) -c  $< 

%.mod: %.f90
	gfortran $(FFLAGS) -c $<

clean:
	rm -f *.o *.exe *.mod *.txt

help:
	@echo "main.exe:	build main program"
	@echo "run:		build and run main program"
