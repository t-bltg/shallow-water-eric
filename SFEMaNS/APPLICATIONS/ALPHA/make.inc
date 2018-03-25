#####################ZAHIR
CMP=mpxlf_r -qfree=f90 -c #Paralle et monoproc
LNK=mpxlf_r 
INC_FFTW=$$FFTW3_I4S4D8
RACINE=../..
SPLIB_DIR = /homegpfs/rech/nor/rnor352/femcode/FEM_CODE/SPLIB/
BLAS_LIB =  -lessl
#####################ZAHIR

OPT      = -O3
#OPT     =   -C #-g
LIB_TIME =


