PROGRAM interpol_neg
!
! CN faire (H,phi) -> -(H,phi)
!
  USE def_type_mesh
  USE prep_maill
  USE chaine_caractere
  USE restart
  USE sub_plot
  USE mesh_interpolation  
  IMPLICIT NONE
  include 'mpif.h'

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

  INTEGER, DIMENSION(:), ALLOCATABLE      :: list_dom_H, list_dom_phi
  INTEGER                                 :: nb_dom_H, nb_dom_phi
  INTEGER                                 :: type_fe_H, type_fe_phi
  LOGICAL                                 :: iformatted
  REAL(KIND=8)                            :: time_h
  CHARACTER(len=3)                        :: type_pb
  INTEGER, DIMENSION(:), ALLOCATABLE      :: list_mode_H
  INTEGER                                 :: k, m, npr_vv, npr_pp, npr_h, npr_phi, m_max, m_max_c, m_max_cr, nb_procs_r
  INTEGER                                 :: code, rang, nb_procs


  !-------------DEBUT PARALLELISATION-----------------------------------------
  CALL MPI_INIT(code)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,code)

  !-------------CONTROL OF DATA INPUTS-------------------------------------------
  OPEN(UNIT = 21, FILE = 'data_interpol_neg', FORM = 'formatted', STATUS = 'unknown')
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

  !-------------ORGANISATION DU MAILLAGE MAXWELL---------------------------------
  IF (type_pb=='mxw' .OR. type_pb=='mhd' .OR. type_pb=='mxx') THEN
     CALL read_until(21, 'mesh_maxwell')
     READ (21, *) type_fe_H, type_fe_phi
     READ (21, *) nb_dom_H  ! number of sub_domains for H
     ALLOCATE(list_dom_H(nb_dom_H))
     READ (21, *) list_dom_H
     READ (21, *) nb_dom_phi  ! number of sub_domains for phi
     ALLOCATE(list_dom_phi(nb_dom_phi))
     READ (21, *) list_dom_phi  


     CALL load_mesh_free_format_ordered(directory, file_name, list_dom_H, type_fe_H, H_mesh, iformatted)  
     CALL load_mesh_free_format_ordered(directory, file_name, list_dom_phi, type_fe_phi, phi_mesh, iformatted)
     CALL load_mesh_free_format_ordered(directory_new, file_name_new, list_dom_H, type_fe_H, H_mesh_new, iformatted)  
     CALL load_mesh_free_format_ordered(directory_new, file_name_new, list_dom_phi, type_fe_phi, phi_mesh_new, iformatted)  
! CARO
     write(*,*) 
     write(*,*) 'H_mesh%np=', H_mesh%np,  'phi_mesh%np=', phi_mesh%np
     write(*,*) 'H_mesh_new%np=', H_mesh_new%np,  'phi_mesh_new%np=', phi_mesh_new%np
     write(*,*) 
! CARO
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
        WRITE(*,*) 'm_max_maxwell, m_max_c_maxwell =', m_max, m_max_c
     IF (MOD(m_max,nb_procs)/=0 .OR. m_max_c==0) THEN
        WRITE(*,*) ' BUG '
        STOP
     END IF
     IF (.NOT.ALLOCATED(list_mode_H)) ALLOCATE(list_mode_H(m_max_c))
     list_mode_H = -1
     ALLOCATE(Hn1  (H_mesh%np,  6,  m_max_c))
     ALLOCATE(Hn   (H_mesh%np,  6,  m_max_c))
     ALLOCATE(phin1(phi_mesh%np,2,  m_max_c))   
     ALLOCATE(phin (phi_mesh%np,2,  m_max_c))
     ALLOCATE(Hn1_new   (H_mesh_new%np,  6,  m_max_c))
     ALLOCATE(Hn_new    (H_mesh_new%np,  6,  m_max_c))
     ALLOCATE(phin1_new (phi_mesh_new%np,2,  m_max_c))   
     ALLOCATE(phin_new  (phi_mesh_new%np,2,  m_max_c))
  END IF
  !-----------INTERPOLATION FOR MAXWELL------------------------------------------
 IF (type_pb=='mxw' .OR. type_pb=='mhd' .OR. type_pb=='mxx') THEN
! CARO
     WRITE(*,*) 'list_mode_H avant= ', list_mode_H
     CALL read_restart_maxwell(H_mesh, phi_mesh, time_h, list_mode_H, Hn, Hn1, &
          phin, phin1, file_name, interpol=.true.)
     WRITE(*,*) 'list_mode_H apres= ', list_mode_H

    Hn_new     = -1.d0*Hn
    Hn1_new    = -1.d0*Hn1
    phin_new   = -1.d0*phin
    phin1_new  = -1.d0*phin1
!
  CALL write_restart_maxwell(H_mesh_new, phi_mesh_new, time_h, list_mode_H, &
       Hn_new, Hn1_new, phin_new, phin1_new, file_name_new, 999, 1)
 END IF 
  !------------------------------------------------------------------------------

  !------------------------FIN------------------------------------------------
  CALL MPI_FINALIZE(code)
  !---------------------------------------------------------------------------

END PROGRAM interpol_neg

