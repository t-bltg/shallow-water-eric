!
!Authors: Jean-Luc Guermond, Copyright 2000
!
MODULE def_type_mesh
  USE dyn_line
  IMPLICIT NONE

  !------------------------------------------------------------------------------
  !  REAL(KIND=8), DIMENSION(n_w,  l_G),  PUBLIC :: ww
  !  REAL(KIND=8), DIMENSION(n_ws, l_Gs), PUBLIC :: wws
  !  REAL(KIND=8), DIMENSION(k_d,  n_w,  l_G,   me),   PUBLIC :: dw
  !  REAL(KIND=8), DIMENSION(n_w,l_G,1:2,me),          PUBLIC :: dwni !d/dn, interface (JLG, April 2009)
  !  REAL(KIND=8), DIMENSION(k_d,        l_Gs,  mes),  PUBLIC :: rnorms
  !  REAL(KIND=8), DIMENSION(l_G,   me),               PUBLIC :: rj
  !  REAL(KIND=8), DIMENSION(l_Gs,  mes),              PUBLIC :: rjs
  !------------------------------------------------------------------------------

  TYPE gauss_type
     INTEGER  :: k_d, n_w, l_G, n_ws, l_Gs
     REAL(KIND=8), DIMENSION(:, :),       POINTER :: ww
     REAL(KIND=8), DIMENSION(:, :),       POINTER :: wws
     REAL(KIND=8), DIMENSION(:, :),       POINTER :: wwsi !Interface shape function (JLG, June 4 2012)
     REAL(KIND=8), DIMENSION(:, :, :, :), POINTER :: dw
     REAL(KIND=8), DIMENSION(:, :, :, :), POINTER :: dwni !Interface gradient (JLG, April 2009)
     REAL(KIND=8), DIMENSION(:,    :, :), POINTER :: rnorms
     REAL(KIND=8), DIMENSION(:, :, :, :), POINTER :: rnormsi !Interface normal (JLG, June 4 2012)
     REAL(KIND=8), DIMENSION(:, :),       POINTER :: rj   !Interface weight (JLG, April 2009)
     REAL(KIND=8), DIMENSION(:, :),       POINTER :: rji 
     REAL(KIND=8), DIMENSION(:, :),       POINTER :: rjs
     REAL(KIND=8), DIMENSION(:, :, :, :), POINTER :: dw_s !gradient at the boundary
     REAL(KIND=8), DIMENSION(:, :, :, :), POINTER :: dwps !special!
     REAL(KIND=8), DIMENSION(:, :, :, :), POINTER :: dws  !SPECIAL!
  END TYPE gauss_type

  TYPE periodic_type
     TYPE(dyn_int_line),  DIMENSION(20)    :: list
     TYPE(dyn_int_line),  DIMENSION(20)    :: perlist
     INTEGER, POINTER,    DIMENSION(: )    :: pnt
     INTEGER                               :: n_bord
  END TYPE periodic_type

  !------------------------------------------------------------------------------
  !  jj(n_w,   me)     nodes of the  volume_elements
  !  jji(n_w, 1:2, mi) edge to node conectivity array --> volume numbering (JLG April 2009)
  !  neighi(1:2, mi)   interfaces to volume elements --> cell 1 has lowest cell number
  !  jjsi(n_ws, mi)    nodes of the interface elements --> volume numbering (JLG April 2009)
  !  jjs(n_ws, mes)    nodes of the surface_elements --> volume numbering
  !  iis(n_ws, mes)    nodes of the surface_elements --> surface numbering 
  !  mm(me)           (possibly sub) set of elements for quadrature
  ! mms(mes)          (possibly sub) set of surface_elements for surf_quadrature
  !------------------------------------------------------------------------------

  TYPE mesh_type
     INTEGER,      POINTER, DIMENSION(:,:) :: jj, jjs, iis
     INTEGER,      POINTER,DIMENSION(:,:,:):: jji  ! (JLG April 2009)
     INTEGER,      POINTER, DIMENSION(:,:) :: jjsi ! (JLG April 2009)
     INTEGER,      POINTER, DIMENSION(:)   :: j_s  ! boundary nodes --> volume numbering
     REAL(KIND=8), POINTER, DIMENSION(:,:) :: rr
     INTEGER,      POINTER, DIMENSION(:,:) :: neigh
     INTEGER,      POINTER, DIMENSION(:,:) :: neighi ! (JLG April 2009)
     INTEGER,      POINTER, DIMENSION(:)   :: sides, neighs
     INTEGER,      POINTER, DIMENSION(:)   :: i_d
     INTEGER                               :: me, mes, np, nps, mi ! (JLG April 2009)
     LOGICAL                               :: edge_stab ! edge stab, yes/no, (JLG April 2009) 
     TYPE(gauss_type)                      :: gauss
     TYPE(periodic_type)                   :: periodic
  END TYPE mesh_type

  TYPE mesh_type_interface
     INTEGER,      POINTER, DIMENSION(:)   :: slave_elem ! list slave elemt in interface 
     INTEGER,      POINTER, DIMENSION(:)   :: list_slave_node ! list of slave nodes on interface 
     INTEGER,      POINTER, DIMENSION(:,:) :: master_node ! local --> global numbering; master nodes
     INTEGER,      POINTER, DIMENSION(:,:) :: slave_node  ! local --> global numbering; slave nodes
     INTEGER                               :: me ! nb of slave elemt in interface  
  END TYPE mesh_type_interface

  TYPE mesh_type_boundary
     INTEGER,      POINTER, DIMENSION(:)   :: master ! list master boundary elemts
     INTEGER,      POINTER, DIMENSION(:)   :: slave ! list slave boundary elemts not in interface
     INTEGER,      POINTER, DIMENSION(:)   :: INTERFACE ! list slave boundary elemts in the interface
     INTEGER,      POINTER, DIMENSION(:,:) :: master_node ! local --> global numbering; master nodes
  END TYPE mesh_type_boundary

  TYPE interface_type
     INTEGER                               :: mes ! number of interface elements 
     INTEGER,      POINTER, DIMENSION(:)   :: mesh1 ! list slave interface elements
     INTEGER,      POINTER, DIMENSION(:)   :: mesh2 ! list master interface elements
     INTEGER,      POINTER, DIMENSION(:,:) :: jjs1 ! list of slave node on interface elements 
     INTEGER,      POINTER, DIMENSION(:,:) :: jjs2 ! list of master nodes on interface elements 
  END TYPE interface_type


END MODULE def_type_mesh
