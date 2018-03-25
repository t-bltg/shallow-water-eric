# on zahir, link this file to as make_local.inc
# with ln make_zahir.inc make_local.ind

RACINE = ../../..
SPLIB_DIR = /homegpfs/rech/nor/rnor522/SPLIB

RACINE_MHD = $(RACINE)/MHD_POSTPROC

MPICMP=mpxlf_r  -c 
MPILNK=mpxlf_r 
MONOCMP=xlf_r -c -I$(RACINE)/FEMSUB_OBJET
MONOLNK=xlf_r 


COMMONOPT = -I$(SPLIB_DIR)  -qfree=f90 $$FFTW3_I4S4D8
RELEASEOPT = -O $(COMMONOPT)
DEBUGOPT = -qnooptimize -g -qfullpath $(COMMONOPT)

LIB_TIME =  

TA_OPT = -L/usr/local/prod/traceanalyzer/lib64 -lVT
