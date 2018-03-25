# zahir: 
# FF=xlf_r  $$FFTW3_I4S4D8 -C

# turing:

CFLAGS =-I/opt/fftw/pgi/include
LDFLAGS = /opt/fftw/pgi/lib/libfftw3.a 

CMP = mpif77
FF= $(CMP) $(CFLAGS)
LNK= $(CMP) 

SRC = test_data.f90 transform1D.f90

OBJ = test_data.o transform1D.o

AOUT = exe

$(AOUT) : $(SRC) $(OBJ)
	$(LNK) $(OBJ) -o $(AOUT) $(LDFLAGS)

test_data.o : test_data.f90
	$(FF) -c test_data.f90

transform1D.o : transform1D.f90
	$(FF) -c transform1D.f90

clean : 
	rm -rf $(OBJ) $(AOUT)
