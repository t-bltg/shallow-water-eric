
# zahir: 
FF=xlf_r  
CFLAGS = $$FFTW3_I4S4D8 -C

# turing:
FF=pgf95
CFLAGS = -C -lfftw3 -lfftw3f -L/opt/fftw/pgi/lib -I/opt/fftw/pgi/include

exe: generate_test_data.f90 
	$(FF) generate_test_data.f90 -o exe $(CFLAGS)

clean : 
	rm -rf $(OBJ) $(AOUT)
