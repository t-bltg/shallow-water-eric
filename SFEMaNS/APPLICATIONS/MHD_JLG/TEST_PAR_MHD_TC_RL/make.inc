#####################Charm (mono)
#CMP	 = gfortran -c -fdefault-real-8
#LNK	 = gfortran
#####################Charm (parallel)
CMP = mpif90 -c -fdefault-real-8
LNK = mpif90
RACINE =/u/guermond/GUERMOND/FEM90/SFEMaNS
INC_FFTW = -I/usr/include/
LIB_FFTW = -l fftw3
SPLIB_DIR = $(RACINE)/../SPLIB
BLAS_LIB =  $(SPLIB_DIR)/blas/localblas.a 
#Utilisation parallel
# mpirun -np 8 ./b.exe -host charm down (different hosts)
# mpirun -np 4 ./b.exe (sur le meme host)
#####################Charm 


OPT     =   -O
LIB_TIME =





