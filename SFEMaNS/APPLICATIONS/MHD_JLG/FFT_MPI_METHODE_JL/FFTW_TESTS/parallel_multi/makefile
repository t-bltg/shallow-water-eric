# LNK=mpxlf_r  -qsuffix=f=f90
##FF=mpxlf_r $$FFTW3_I4S4D8 -C -qfree=f90
FF=gfortran -C -l LD_LIBRARY_PATH $(LIBRARY_PATH)

SRC = test_data.f90 multi_transform.f90

OBJ = test_data.o multi_transform.o

AOUT = exe

$(AOUT) : $(SRC) $(OBJ)
	$(FF) $(OBJ) -o $(AOUT) 

test_data.o : test_data.f90
	$(FF) -c test_data.f90

multi_transform.o : multi_transform.f90
	$(FF) -c multi_transform.f90

clean : 
	rm -rf $(OBJ) $(AOUT) *.mod
