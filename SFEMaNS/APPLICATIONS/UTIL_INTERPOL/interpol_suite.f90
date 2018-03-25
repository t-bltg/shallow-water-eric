PROGRAM  interpol
  USE def_type_mesh
  USE prep_maill
  USE chaine_caractere
  USE restart
  USE sub_plot
  USE mesh_interpolation  
  IMPLICIT NONE
  include 'mpif.h'

  !Champs pour NS-------------------------------------------------------------
  TYPE(mesh_type)                                 :: pp_mesh, vv_mesh  
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:)     :: un, un_m1
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:)     :: pn, pn_m1
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:)     :: incpn, incpn_m1
  TYPE(mesh_type)                                 :: pp_mesh_new, vv_mesh_new 
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:)     :: un_new, un_m1_new
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:)     :: pn_new, pn_m1_new
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:)     :: incpn_new, incpn_m1_new
  !---------------------------------------------------------------------------

  !Champs pour Maxwell--------------------------------------------------------
  TYPE(mesh_type), TARGET                         :: H_mesh, phi_mesh
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:)     :: Hn, Hn1, phin, phin1
  TYPE(mesh_type), TARGET                         :: H_mesh_new, phi_mesh_new
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:)     :: Hn_new, Hn1_new, phin_new, phin1_new
  !---------------------------------------------------------------------------

  !Nom du fichier pour restart-----------------------------------------------
  CHARACTER(len=64)                               :: directory, file_name
  CHARACTER(len=64)                               :: directory_new, file_name_new
  !---------------------------------------------------------------------------

  INTEGER, DIMENSION(:), ALLOCATABLE      :: list_dom_ns, list_dom_H, list_dom_H_in, list_dom_phi
  INTEGER                                 :: nb_dom_ns, nb_dom_H, nb_dom_phi
  INTEGER                                 :: type_fe_H, type_fe_phi
  LOGICAL                                 :: iformatted
  REAL(KIND=8)                            :: time_u, time_h
  CHARACTER(len=3)                        :: type_pb
  INTEGER, DIMENSION(:), ALLOCATABLE      :: list_mode
  INTEGER                                 :: k, m, npr_vv, npr_pp, npr_h, npr_phi, m_max, m_max_c, m_max_cr, nb_procs_r
  INTEGER                                 :: code, rang, nb_procs

  REAL(KIND=8), DIMENSION(2) :: r_0, r_1

  CALL MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,code)

  !-------------CONTROL OF DATA INPUTS-------------------------------------------
  OPEN(UNIT = 21, FILE = 'data_interpol', FORM = 'formatted', STATUS = 'unknown')
  !------------------------------------------------------------------------------

  !-------------TYPE OF THE PROBLEM----------------------------------------------
  CALL read_until(21, 'problem_type')
  READ (21, *) type_pb
  !------------------------------------------------------------------------------

  !-------------LE VIEUX MAILLAGE------------------------------------------------
  CALL read_until(21, 'data_old_mesh')
  READ (21, *) iformatted
  READ (21, *) directory, file_name
  !------------------------------------------------------------------------------

  !-------------LE NOUVEAU MAILLAGE------------------------------------------------
  CALL read_until(21, 'data_new_mesh')
  READ (21, *) iformatted
  READ (21, *) directory_new, file_name_new
  !------------------------------------------------------------------------------

  !-------------ORGANISATION DU MAILLAGE NAVIER_STOKES---------------------------
  IF (type_pb=='nst' .OR. type_pb=='mhd' .OR. type_pb=='mxx') THEN
     CALL read_until(21, 'mesh_navier_stokes')
     READ (21, *)  nb_dom_ns 
     ALLOCATE(list_dom_ns(nb_dom_ns))
     READ (21, *)  list_dom_ns
     CALL load_mesh_free_format_ordered(directory, file_name, list_dom_ns, 1, pp_mesh, iformatted)
     CALL load_mesh_free_format_ordered(directory, file_name, list_dom_ns, 2, vv_mesh, iformatted)
     CALL load_mesh_free_format_ordered(directory_new, file_name_new, list_dom_ns, 1, pp_mesh_new, iformatted)
     CALL load_mesh_free_format_ordered(directory_new, file_name_new, list_dom_ns, 2, vv_mesh_new, iformatted)
  END IF
  !------------------------------------------------------------------------------

  !-------------ORGANISATION DU MAILLAGE MAXWELL---------------------------------
  IF (type_pb=='mxw' .OR. type_pb=='mhd' .OR. type_pb=='mxx') THEN
     CALL read_until(21, 'mesh_maxwell')
     READ (21, *) type_fe_H, type_fe_phi
     READ (21, *) nb_dom_H  ! number of sub_domains for H
     ALLOCATE(list_dom_H_in(nb_dom_H), list_dom_H(nb_dom_H))
     READ (21, *) list_dom_H_in
     READ (21, *) nb_dom_phi  ! number of sub_domains for phi
     ALLOCATE(list_dom_phi(nb_dom_phi))
     READ (21, *) list_dom_phi  

     IF (type_pb=='mhd' .OR.  type_pb=='mxx') THEN
        IF (SIZE(list_dom_H) < SIZE(list_dom_ns)) THEN
           WRITE(*,*) ' BUG: NS must be a subset of Maxwell ' 
           STOP
        END IF
        DO k = 1, nb_dom_ns 
           IF (MINVAL(ABS(list_dom_H_in - list_dom_ns(k))) /= 0) THEN
              WRITE(*,*) ' BUG : NS must be a subset of Maxwell '
              STOP
           END IF
           list_dom_H(k) = list_dom_ns(k)
        END DO
        m = nb_dom_ns 
        DO k = 1, nb_dom_H
           IF (MINVAL(ABS(list_dom_H_in(k) - list_dom_ns)) == 0) CYCLE 
           m = m + 1
           list_dom_H(m) = list_dom_H_in(k)
        END DO
        IF (m/=nb_dom_H) THEN
           WRITE(*,*) ' BUG : m/=nb_dom_H ' 
           STOP
        END IF
     ELSE
        list_dom_H = list_dom_H_in
     END IF

     CALL load_mesh_free_format_ordered(directory, file_name, list_dom_H, type_fe_H, H_mesh, iformatted)  
     CALL load_mesh_free_format_ordered(directory, file_name, list_dom_phi, type_fe_phi, phi_mesh, iformatted)
     CALL load_mesh_free_format_ordered(directory_new, file_name_new, list_dom_H, type_fe_H, H_mesh_new, iformatted)  
     CALL load_mesh_free_format_ordered(directory_new, file_name_new, list_dom_phi, type_fe_phi, phi_mesh_new, iformatted)  
  END IF
  !------------------------------------------------------------------------------

  !------------ARRAY ALLOCATION FOR NAVIER STOKES--------------------------------
  IF (type_pb=='nst' .OR. type_pb=='mhd' .OR. type_pb=='mxx') THEN
     OPEN(UNIT = 10, FILE = 'suite_ns.'//file_name, FORM = 'unformatted', STATUS = 'unknown')
     READ(10) time_u, npr_vv , npr_pp, nb_procs_r, m_max_cr
     CLOSE(10)
     IF (npr_vv/=vv_mesh%np .or. npr_pp/=pp_mesh%np) THEN
        WRITE(*,*) ' Fichier suite_ns incompatible avec le maillage'
        STOP
     END IF

     m_max = nb_procs_r*m_max_cr 
     m_max_c = m_max/nb_procs
     IF (MOD(m_max,nb_procs)/=0 .OR. m_max_c==0) THEN
        WRITE(*,*) ' BUG '
        STOP
     END IF
     IF (.NOT.ALLOCATED(list_mode)) ALLOCATE(list_mode(m_max_c))
     list_mode = -1
     ALLOCATE(un_m1   (vv_mesh%np, 6, m_max_c))
     ALLOCATE(un      (vv_mesh%np, 6, m_max_c))
     ALLOCATE(pn_m1   (pp_mesh%np, 2, m_max_c))
     ALLOCATE(pn      (pp_mesh%np, 2, m_max_c))
     ALLOCATE(incpn_m1(pp_mesh%np, 2, m_max_c))
     ALLOCATE(incpn   (pp_mesh%np, 2, m_max_c))
     ALLOCATE(un_m1_new   (vv_mesh_new%np, 6, m_max_c))
     ALLOCATE(un_new      (vv_mesh_new%np, 6, m_max_c))
     ALLOCATE(pn_m1_new   (pp_mesh_new%np, 2, m_max_c))
     ALLOCATE(pn_new      (pp_mesh_new%np, 2, m_max_c))
     ALLOCATE(incpn_m1_new(pp_mesh_new%np, 2, m_max_c))
     ALLOCATE(incpn_new   (pp_mesh_new%np, 2, m_max_c))
  END IF
  !------------------------------------------------------------------------------


  !------------ARRAY ALLOCATION FOR MAXWELL---------------------------------------
  IF (type_pb=='mxw' .OR. type_pb=='mhd' .OR. type_pb=='mxx' ) THEN
     OPEN(UNIT = 10, FILE = 'suite_maxwell.'//file_name, FORM = 'unformatted', STATUS = 'unknown')
     READ(10) time_h, npr_h , npr_phi, nb_procs_r, m_max_cr
     CLOSE(10)
     IF (npr_h/=H_mesh%np .or. npr_phi/=phi_mesh%np) THEN
        WRITE(*,*) ' Fichier suite_maxwell incompatible avec le maillage'
        STOP
     END IF

     m_max = nb_procs_r*m_max_cr 
     m_max_c = m_max/nb_procs
     IF (MOD(m_max,nb_procs)/=0 .OR. m_max_c==0) THEN
        WRITE(*,*) ' BUG '
        STOP
     END IF
     IF (.NOT.ALLOCATED(list_mode)) ALLOCATE(list_mode(m_max_c))
     list_mode = -1
     ALLOCATE(Hn1  (H_mesh%np,  6,  m_max_c))
     ALLOCATE(Hn   (H_mesh%np,  6,  m_max_c))
     ALLOCATE(phin1(phi_mesh%np,2,  m_max_c))   
     ALLOCATE(phin (phi_mesh%np,2,  m_max_c))
     ALLOCATE(Hn1_new   (H_mesh_new%np,  6,  m_max_c))
     ALLOCATE(Hn_new    (H_mesh_new%np,  6,  m_max_c))
     ALLOCATE(phin1_new (phi_mesh_new%np,2,  m_max_c))   
     ALLOCATE(phin_new  (phi_mesh_new%np,2,  m_max_c))
  END IF
  !------------------------------------------------------------------------------

  !-----------INTERPOLATION FOR NAVIER-STOKES------------------------------------
  IF (type_pb=='nst' .OR. type_pb=='mhd' .OR. type_pb=='mxx') THEN
     CALL read_restart_ns(vv_mesh, pp_mesh, time_u, list_mode, un, un_m1, pn, pn_m1, &
          incpn, incpn_m1, file_name, interpol=.true.)

     CALL mesh_interp(pp_mesh,pp_mesh_new,pn(:,1,1),pn_new(:,1,1),.true.)
     DO m = 1, m_max_c
        DO k = 1, 2
           CALL mesh_interp(pp_mesh,pp_mesh_new,      pn(:,k,m),      pn_new(:,k,m),.false.)
           CALL mesh_interp(pp_mesh,pp_mesh_new,   pn_m1(:,k,m),   pn_m1_new(:,k,m),.false.)
           CALL mesh_interp(pp_mesh,pp_mesh_new,   incpn(:,k,m),   incpn_new(:,k,m),.false.)
           CALL mesh_interp(pp_mesh,pp_mesh_new,incpn_m1(:,k,m),incpn_m1_new(:,k,m),.false.)
        END DO
     END DO

     CALL mesh_interp(vv_mesh,vv_mesh_new,un(:,1,1),un_new(:,1,1),.true.)
     DO m = 1, m_max_c
        DO k =1, 6
           CALL mesh_interp(vv_mesh,vv_mesh_new,   un(:,k,m),   un_new(:,k,m),.false.)
           CALL mesh_interp(vv_mesh,vv_mesh_new,un_m1(:,k,m),un_m1_new(:,k,m),.false.)
        END DO
     END DO
     CALL write_restart_ns(vv_mesh_new, pp_mesh_new, time_u, list_mode, &
          un_new, un_m1_new, pn_new, pn_m1_new, incpn_new, incpn_m1_new, file_name_new, 999, 1)
  END IF
  !------------------------------------------------------------------------------

  !-----------INTERPOLATION FOR MAXWELL------------------------------------------
  IF (type_pb=='mxw' .OR. type_pb=='mhd' .OR. type_pb=='mxx') THEN
     CALL read_restart_maxwell(H_mesh, phi_mesh, time_h, list_mode, Hn, Hn1, &
          phin, phin1, file_name)

     CALL mesh_interp(phi_mesh,phi_mesh_new,phin(:,1,1),phin_new(:,1,1),.true.)
     DO m = 1, m_max_c
        DO k = 1, 2
           CALL mesh_interp(phi_mesh,phi_mesh_new,    phin(:,k,m),    phin_new(:,k,m),.false.)
           CALL mesh_interp(phi_mesh,phi_mesh_new,   phin1(:,k,m),   phin1_new(:,k,m),.false.)
        END DO
     END DO

     CALL mesh_interp(H_mesh,H_mesh_new,Hn(:,1,1),Hn_new(:,1,1),.true.)
     DO m = 1, m_max_c
        DO k =1, 6
           CALL mesh_interp(H_mesh,H_mesh_new, Hn(:,k,m), Hn_new(:,k,m),.false.)
           CALL mesh_interp(H_mesh,H_mesh_new,Hn1(:,k,m),Hn1_new(:,k,m),.false.)
        END DO
     END DO
  END IF 
  CALL write_restart_maxwell(H_mesh_new, phi_mesh_new, time_h, list_mode, &
       Hn_new, Hn1_new, phin_new, phin1_new, file_name_new, 999, 1)
  !------------------------------------------------------------------------------


END PROGRAM interpol
