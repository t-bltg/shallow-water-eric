include make.inc

FEM_OBJ =  dyn_line.o \
           def_type_mesh.o \
	   associate_gauss.o \
	   gauss_points_2d_p1.o \
	   gauss_points_2d_p2.o \
	   chaine_caractere.o \
	   dir_nodes.o \
	   prep_mesh.o \
	   prep_mesh_interface.o \
	   prep_mesh_periodic.o \
	   fem_s_direct_axi.o \
	   fem_tn_axi.o \
	   fem_tn_navier_mhd.o \
	   tn_par.o \
	   fem_sparsekit_axi.o \
	   st_sparsekit.o \
	   sub_plot.o \
	   bessel.o \
	   condlim_test_conv.o\
           condlim.o \
	   restart.o \
	   fft_parallel.o \
	   sub_ns.o \
	   maxwell_update_time.o \
	   initialisation.o \
	   user_time.o \
       post_proc.o

LOC_OBJ =  

GAUSS_OBJ = gauss_points_2d_p1.o\
	   gauss_points_2d_p2.o

SPLIB_OBJ = $(SPLIB_DIR)/orderings/reorder.o \
	$(SPLIB_DIR)/solve_splib.o \
	$(SPLIB_DIR)/splib.a 

####################################
# Version multiprocesseur 
####################################
test.exe: $(FEM_OBJ) $(SPLIB_OBJ) $(LOC_OBJ) main_mhd.o
	($(LNK) $(OPT) $(LIB_TIME) -o test.exe main_mhd.o $(FEM_OBJ) $(SPLIB_OBJ) $(LOC_OBJ) \
	$(BLAS_LIB) $(LIB_FFTW)) 

main_mhd.o:  $(FEM_OBJ)  main_mhd.f90 
	($(CMP) $(OPT) $(LIB_TIME) $(INC_FFTW) main_mhd.f90)

###################################
# Execution paralelle sur Zahir
###################################
r1: test.exe
	mpirun -np 1 test.exe
r2: test.exe
	mpirun -np 2 test.exe

###################################
# EXECUTION paralelle sur Turing
###################################
run1: test.exe
	mpiexec -boot -machinefile hostfile1 -n 1 ./test.exe
run2: test.exe
	mpiexec -boot -machinefile hostfile2 -n 2 ./test.exe
run3: test.exe
	mpiexec -boot -machinefile hostfile3 -n 3 ./test.exe
run4: test.exe
	mpiexec -boot -machinefile hostfile4 -n 4 ./test.exe

####################################
# Version mono proceceur
####################################
test_mono.exe: mpi.o  mpif.h $(FEM_OBJ) $(SPLIB_OBJ) $(LOC_OBJ) main_mhd.o
	($(LNK) $(OPT) $(LIB_TIME) -o test_mono.exe main_mhd.o mpi.o $(FEM_OBJ) $(SPLIB_OBJ) \
	$(LOC_OBJ) $(BLAS_LIB) $(LIB_FFTW); rm mpif.h)

mpi.o: $(RACINE)/FEMSUB_OBJET/mpi.f90 
	$(CMP) $< -o $(*).o

mpif.h: $(RACINE)/FEMSUB_OBJET/mpif.h
	cp $(RACINE)/FEMSUB_OBJET/mpif.h .

####################################
# Menage
####################################
clean: 
	rm *.o *.mod suite* *.plt


#####################################
# FEM OBJECTS       
#####################################

gauss_points_2d_p1.o: $(RACINE)/GAUSS_POINTS_OBJET/gauss_points_2d_p1.f90
	$(CMP) $(OPT) $< -o $(*).o

gauss_points_2d_p2.o: $(RACINE)/GAUSS_POINTS_OBJET/gauss_points_2d_p2.f90
	$(CMP) $(OPT) $< -o $(*).o

associate_gauss.o: $(RACINE)/FEMSUB_OBJET/associate_gauss.f90
	$(CMP) $(OPT) $< -o $(*).o

user_time.o: $(RACINE)/FEMSUB/user_time.f90
	$(CMP) $(OPT) $< -o $(*).o

chaine_caractere.o: $(RACINE)/FEMSUB/chaine_caractere.f90
	$(CMP) $(OPT) $< -o $(*).o

dyn_line.o: $(RACINE)/FEMSUB/dyn_line.f90
	$(CMP) $(OPT) $< -o $(*).o

dir_nodes.o: $(RACINE)/FEMSUB/dir_nodes.f90 dyn_line.o
	$(CMP) $(OPT) $< -o $(*).o

def_type_mesh.o: $(RACINE)/FEMSUB_OBJET/def_type_mesh.f90
	$(CMP) $(OPT) $< -o $(*).o

prep_mesh.o: $(RACINE)/FEMSUB_OBJET/prep_mesh.f90 $(GAUSS_OBJ)
	$(CMP) $(OPT) $< -o $(*).o

prep_mesh_periodic.o: $(RACINE)/FEMSUB_OBJET/prep_mesh_periodic.f90 $(GAUSS_OBJ)
	cp $(SPLIB_DIR)/solve_sp.mod .; $(CMP) $(OPT) $< -o $(*).o 

prep_mesh_interface.o: $(RACINE)/FEMSUB_OBJET/prep_mesh_interface.f90
	$(CMP) $(OPT) $< -o $(*).o

fem_s_direct_axi.o: $(RACINE)/FEMSUB_OBJET/fem_s_direct_axi.f90 $(GAUSS_OBJ)
	$(CMP) $(OPT) $< -o $(*).o

fem_tn_axi.o: $(RACINE)/FEMSUB_OBJET/fem_tn_axi.f90 $(GAUSS_OBJ)
	$(CMP) $(OPT) $< -o $(*).o

fem_sparsekit_axi.o: $(RACINE)/FEMSUB_OBJET/fem_sparsekit_axi.f90 $(GAUSS_OBJ)
	$(CMP) $(OPT) $< -o $(*).o

st_sparsekit.o: $(RACINE)/FEMSUB/st_sparsekit.f90
	$(CMP) $(OPT) $< -o $(*).o

sub_plot.o: $(RACINE)/FEMSUB/sub_plot.f90
	$(CMP) $(OPT) $< -o $(*).o

condlim.o: ./condlim.f90
	$(CMP) $(OPT) $< -o $(*).o

condlim_test_conv.o: $(RACINE)/MHD_DATA_TEST_CONV/condlim_test_conv.f90
	$(CMP) $(OPT) $< -o $(*).o

fft_parallel.o: $(RACINE)/FEMSUB_OBJET/fft_parallel.f90
	$(CMP) $(OPT) $< -o $(*).o $(INC_FFTW)

fem_tn_navier_mhd.o: $(RACINE)/FEMSUB_OBJET/fem_tn_navier_mhd.f90
	$(CMP) $(OPT) $< -o $(*).o

tn_par.o: $(RACINE)/FEMSUB_OBJET/tn_par.f90 fem_tn_navier_mhd.o
	$(CMP) $(OPT) $< -o $(*).o

restart.o: $(RACINE)/FEMSUB_OBJET/restart.f90
	$(CMP) $(OPT) $< -o $(*).o

sub_ns.o: $(RACINE)/FEMSUB_OBJET/sub_ns.f90
	$(CMP) $(OPT) $< -o $(*).o

maxwell_update_time.o: $(RACINE)/FEMSUB_OBJET/maxwell_update_time.f90
	$(CMP) $(OPT) $< -o $(*).o

initialisation.o: $(RACINE)/FEMSUB_OBJET/initialisation.f90
	$(CMP) $(OPT) $< -o $(*).o

bessel.o: $(RACINE)/FEMSUB/bessel.f90
	$(CMP) $(OPT) $< -o $(*).o

post_proc.o: $(RACINE)/MHD_POSTPROC/post_proc.f90
	$(CMP) $(OPT) $< -o $(*).o

