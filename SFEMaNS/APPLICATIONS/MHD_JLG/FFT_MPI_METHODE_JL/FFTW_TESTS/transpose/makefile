# turing:
CFLAGS =-I/opt/fftw/pgi/include
LDFLAGS = /opt/fftw/pgi/lib/libfftw3.a 
CMP = mpif77
LNK= mpif77

# zahir: 
CMP = mpxlf_r
LNK= mpxlf_r
CFLAGS =  $$FFTW3_I4S4D8 -C -qfree=f90
LDFLAGS =  $$FFTW3_I4S4D8 -C -qfree=f90

#############

SRC = transpose.f90

OBJ = transpose.o

AOUT = exe

$(AOUT) : $(SRC) $(OBJ)
	$(LNK) $(OBJ) -o $(AOUT) $(LDFLAGS)

transpose.o : transpose.f90
	$(CMP)  $(CFLAGS) -c $<

clean : 
	rm -rf $(OBJ) $(AOUT)
