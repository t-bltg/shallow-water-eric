####################TURING
#CMP = mpif77 -c -r8    $(FFTW_PGI_INCLUDE)
#LNK = mpif77
#RACINE =/home/coro/guermond/SFEMaNS
#INC_FFTW = $(FFTW_PGI_INCLUDE)
#LIB_FFTW = $(FFTW_PGI_64)
####################TURING

#####################zahir
#CMP = mpxlf_r -c -qfree=f90  -qrealsize=8
#LNK = mpxlf_r -qfree=f90  -qrealsize=8
#RACINE =/homegpfs/rech/nor/rnor522/KB/SFEMaNS
#INC_FFTW = $$FFTW3_I4S4D8
#LIB_FFTW = $$FFTW3_I4S4D8
#####################zahir

#####################Charm (mono)
#CMP	 = gfortran -c -fdefault-real-8
#LNK	 = gfortran
#####################Charm (parallel)
CMP = mpif90 -c -fdefault-real-8
LNK = mpif90
RACINE =/u/guermond/GUERMOND/FEM90/SFEMaNS
INC_FFTW = -I/usr/include/
LIB_FFTW = -l fftw3
#Utilisation parallel
# mpirun -np 8 ./b.exe -host charm down (different hosts)
# mpirun -np 4 ./b.exe (sur le meme host)
#####################Charm 

OPT= -O

####################################
# Version mono proceceur
####################################
mono.exe:  mpif.h mpi.o fft_parallel.o test_fft_par.o 
	$(LNK) -o mono.exe fft_parallel.o  mpi.o test_fft_par.o $(LIB_FFTW) 

mpi.o: $(RACINE)/FEMSUB_OBJET/mpi.f90
	(cp $(RACINE)/FEMSUB_OBJET/mpi.f90 mpi.f90; \
	$(CMP) $(OPT) mpi.f90; \
	rm mpi.f90)

mpif.h: $(RACINE)/FEMSUB_OBJET/mpif.h
	(cp $(RACINE)/FEMSUB_OBJET/mpif.h mpif.h)
####################################
# Version mono proceceur
####################################

b.exe : fft_parallel.o test_fft_par.o
	$(LNK) -o b.exe test_fft_par.o fft_parallel.o $(LIB_FFTW)

exec=b.exe

run1: $(exec)
	mpiexec -boot -machinefile hostfile1 -n 1 ./$(exec)
run2: $(exec)
	mpiexec -boot -machinefile hostfile2 -n 2 ./$(exec)
run3: $(exec)
	mpiexec -boot -machinefile hostfile3 -n 3 ./$(exec)
run4: $(exec)
	mpiexec -boot -machinefile hostfile4 -n 4 ./$(exec)
run8: $(exec)
	mpiexec -boot -machinefile hostfile8 -n 8 ./$(exec)


fft_parallel.o : $(RACINE)/FEMSUB_OBJET/fft_parallel.f90
	$(CMP) $(OPT) $(INC_FFTW) $(RACINE)/FEMSUB_OBJET/fft_parallel.f90 

test_fft_par.o : test_fft_par.f90
		$(CMP) $(OPT) $(INC_FFTW) $< -o $(*).o


clean :
	rm -rf *.mod *.o


