include make.inc

a.exe: main.o space_dim.o read_data.o tools.o mesh.o def_type_mesh.o dyn_line.o dir_nodes.o prep_mesh.o chaine_caractere.o gauss_points_2d_p1_april_2009.o gauss_points_2d_p2_april_2009.o toto.o type_matrix.o st_sparsekit.o fem_sparsekit.o associate_gauss.o condlim.o solver_pardiso.o sub_plot.o fem_tn.o
	$(LNK) -o a.exe  main.o space_dim.o read_data.o tools.o mesh.o def_type_mesh.o dyn_line.o dir_nodes.o prep_mesh.o chaine_caractere.o gauss_points_2d_p1_april_2009.o gauss_points_2d_p2_april_2009.o toto.o type_matrix.o st_sparsekit.o fem_sparsekit.o associate_gauss.o condlim.o solver_pardiso.o sub_plot.o -L$(LIB_PARDISO) fem_tn.o
# Here are the compile steps

main.o:main.f90 space_dim.o read_data.o mesh.o toto.o condlim.o sub_plot.o fem_tn.o associate_gauss.o
	$(F90) -c  main.f90

space_dim.o:./space_dim.f90
	$(F90) -c  ./space_dim.f90

read_data.o:./read_data.f90 tools.o
	$(F90) -c  ./read_data.f90

tools.o:./tools.f90
	$(F90) -c  ./tools.f90

mesh.o:./mesh.f90 def_type_mesh.o dir_nodes.o space_dim.o read_data.o prep_mesh.o
	$(F90) -c  ./mesh.f90

def_type_mesh.o:$(ROOT)/FEMSUB_OBJET/def_type_mesh.f90 dyn_line.o
	$(F90) -c  $(ROOT)/FEMSUB_OBJET/def_type_mesh.f90

dyn_line.o:$(ROOT)/FEMSUB/dyn_line.f90
	$(F90) -c  $(ROOT)/FEMSUB/dyn_line.f90

dir_nodes.o:$(ROOT)/FEMSUB/dir_nodes.f90
	$(F90) -c  $(ROOT)/FEMSUB/dir_nodes.f90

prep_mesh.o:$(ROOT)/FEMSUB_OBJET/prep_mesh.f90 def_type_mesh.o chaine_caractere.o gauss_points_2d_p1_april_2009.o gauss_points_2d_p2_april_2009.o dir_nodes.o
	$(F90) -c  $(ROOT)/FEMSUB_OBJET/prep_mesh.f90

chaine_caractere.o:$(ROOT)/FEMSUB/chaine_caractere.f90
	$(F90) -c  $(ROOT)/FEMSUB/chaine_caractere.f90

gauss_points_2d_p1_april_2009.o:$(ROOT)/GAUSS_POINTS_OBJET/gauss_points_2d_p1_april_2009.f90 def_type_mesh.o
	$(F90) -c  $(ROOT)/GAUSS_POINTS_OBJET/gauss_points_2d_p1_april_2009.f90

gauss_points_2d_p2_april_2009.o:$(ROOT)/GAUSS_POINTS_OBJET/gauss_points_2d_p2_april_2009.f90 def_type_mesh.o
	$(F90) -c  $(ROOT)/GAUSS_POINTS_OBJET/gauss_points_2d_p2_april_2009.f90

toto.o:./toto.f90 type_matrix.o space_dim.o read_data.o st_sparsekit.o mesh.o fem_sparsekit.o condlim.o solver_pardiso.o sub_plot.o
	$(F90) -c  ./toto.f90

type_matrix.o:$(ROOT)/FEMSUB_OBJET/type_matrix.f90
	$(F90) -c  $(ROOT)/FEMSUB_OBJET/type_matrix.f90

st_sparsekit.o:$(ROOT)/FEMSUB/st_sparsekit.f90 def_type_mesh.o
	$(F90) -c  $(ROOT)/FEMSUB/st_sparsekit.f90

fem_sparsekit.o:$(ROOT)/FEMSUB_OBJET/fem_sparsekit.f90 associate_gauss.o def_type_mesh.o
	$(F90) -c  $(ROOT)/FEMSUB_OBJET/fem_sparsekit.f90

associate_gauss.o:$(ROOT)/FEMSUB_OBJET/associate_gauss.f90 def_type_mesh.o
	$(F90) -c  $(ROOT)/FEMSUB_OBJET/associate_gauss.f90

condlim.o:./condlim.f90 read_data.o mesh.o space_dim.o
	$(F90) -c  ./condlim.f90

solver_pardiso.o:$(ROOT)/FEMSUB_OBJET/solver_pardiso.f90
	$(F90) -c  $(ROOT)/FEMSUB_OBJET/solver_pardiso.f90 -I$(INC_PARDISO)

sub_plot.o:$(ROOT)/FEMSUB/sub_plot.f90 def_type_mesh.o chaine_caractere.o
	$(F90) -c  $(ROOT)/FEMSUB/sub_plot.f90

fem_tn.o:$(ROOT)/FEMSUB_OBJET/fem_tn.f90 associate_gauss.o
	$(F90) -c  $(ROOT)/FEMSUB_OBJET/fem_tn.f90

# This entry allows you to type " make cl " to get rid of
# all object and module files
cl:
	 rm -f -r f_{files,modd}* *.o *.mod *.M *.d V*.inc *.vo  V*.f *.dbg album F.err
