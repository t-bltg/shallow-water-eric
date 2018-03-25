#####################TURING
CMP=mpif90 -c #Parallele et monoproc
LNK=mpif90
INC_FFTW=$(FFTW_PGI_INCLUDE) $(FFTW_PGI_64)
RACINE =/u/guermond/GUERMOND/FEM90/SFEMaNS
SPLIB_DIR = /u/guermond/GUERMOND/FEM90/SPLIB
#SPLIB_DIR = /home/nore/nore00/SPLIB
BLAS_LIB =  $(ACML_PGI_64)
LIB_PARDISO = $(MKL_INTEL_64)
INC_PARDISO = /opt/intel/Compiler/11.0/083/mkl/include 
#####################TURING

OPT     =   -O
LIB_TIME =


