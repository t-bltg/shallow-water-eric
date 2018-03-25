ROOT  = /u/guermond/GUERMOND/FEM90/SFEMaNS
OPT   = -O #-Wall -fimplicit-none -fbounds-check 
##OPT   = -march=native -mtune=native -Ofast
F90   = gfortran $(OPT) 
LNK   = gfortran $(OPT)

RACINE_INTEL = /w/guermond
LIB_PARDISO = $(RACINE_INTEL)/intel/composerxe-2011.0.084/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
INC_PARDISO = $(RACINE_INTEL)/intel/composerxe-2011.0.084/mkl/include

SPLIB_DIR = /u/guermond/GUERMOND/FEM90/SPLIB
SPLIB_OBJ = $(SPLIB_DIR)/splib.a \
	$(SPLIB_DIR)/blas/localblas.a \
	$(SPLIB_DIR)/orderings/reorder.o
