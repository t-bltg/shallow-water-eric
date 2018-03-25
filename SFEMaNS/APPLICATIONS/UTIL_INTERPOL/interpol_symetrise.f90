PROGRAM interpol_sym
!
! CN essai de symetriser un champ de restart pour preparer CI non lineaire
!
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
  ! SYM
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:)     :: un_sym, un_m1_sym
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:)     :: pn_sym, pn_m1_sym
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:)     :: incpn_sym, incpn_m1_sym
  INTEGER, ALLOCATABLE, DIMENSION(:)              :: point_sym_u, point_sym_p
  INTEGER, ALLOCATABLE, DIMENSION(:)              :: diff_point_u, diff_point_p
  REAL(KIND=8),DIMENSION(3)                       :: type_sym_u
  REAL(KIND=8),DIMENSION(1)                       :: type_sym_p
  !---------------------------------------------------------------------------

  !Champs pour Maxwell--------------------------------------------------------
  TYPE(mesh_type), TARGET                         :: H_mesh, phi_mesh
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:)     :: Hn, Hn1, phin, phin1
  TYPE(mesh_type), TARGET                         :: H_mesh_new, phi_mesh_new
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:)     :: Hn_new, Hn1_new, phin_new, phin1_new
  ! SYM
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:)     :: Hn_sym, Hn1_sym, phin_sym, phin1_sym
  INTEGER, ALLOCATABLE, DIMENSION(:)              :: point_sym_H, point_sym_phi
  INTEGER, ALLOCATABLE, DIMENSION(:)              :: diff_point_H, diff_point_phi
  REAL(KIND=8),DIMENSION(3)                       :: type_sym_H
  REAL(KIND=8),DIMENSION(1)                       :: type_sym_phi
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
  INTEGER, DIMENSION(:), ALLOCATABLE      :: list_mode_u, list_mode_H
  INTEGER                                 :: k, m, npr_vv, npr_pp, npr_h, npr_phi, m_max, m_max_c, m_max_cr, nb_procs_r
  INTEGER                                 :: code, rang, nb_procs

  REAL(KIND=8), DIMENSION(2) :: r_0, r_1
  REAL(KIND=8)               :: eps =1.d-4
  INTEGER                    :: i, find_point, i_compteur

  !-------------DEBUT PARALLELISATION-----------------------------------------
  CALL MPI_INIT(code)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,code)

  !-------------CONTROL OF DATA INPUTS-------------------------------------------
  OPEN(UNIT = 21, FILE = 'data_interpol_sym', FORM = 'formatted', STATUS = 'unknown')
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
! CARO
     write(*,*) 
     write(*,*) 'pp_mesh%np=', pp_mesh%np,  'vv_mesh%np=', vv_mesh%np
     write(*,*) 'pp_mesh_new%np=', pp_mesh_new%np,  'vv_mesh_new%np=', vv_mesh_new%np
     write(*,*)
! CARO


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
! CARO
     write(*,*) 
     write(*,*) 'H_mesh%np=', H_mesh%np,  'phi_mesh%np=', phi_mesh%np
     write(*,*) 'H_mesh_new%np=', H_mesh_new%np,  'phi_mesh_new%np=', phi_mesh_new%np
     write(*,*) 
! CARO
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
        WRITE(*,*) 'm_max_vitesse, m_max_c_vitesse =', m_max, m_max_c
     IF (MOD(m_max,nb_procs)/=0 .OR. m_max_c==0) THEN
        WRITE(*,*) ' BUG '
        STOP
     END IF
     IF (.NOT.ALLOCATED(list_mode_u)) ALLOCATE(list_mode_u(m_max_c))
     list_mode_u = -1
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
  ! SYM
     ALLOCATE(un_m1_sym   (vv_mesh_new%np, 6, m_max_c))
     ALLOCATE(un_sym      (vv_mesh_new%np, 6, m_max_c))
     ALLOCATE(pn_m1_sym   (pp_mesh_new%np, 2, m_max_c))
     ALLOCATE(pn_sym      (pp_mesh_new%np, 2, m_max_c))
     ALLOCATE(incpn_m1_sym(pp_mesh_new%np, 2, m_max_c))
     ALLOCATE(incpn_sym   (pp_mesh_new%np, 2, m_max_c))
     ALLOCATE(point_sym_u(vv_mesh_new%np),point_sym_p(pp_mesh_new%np))
     ALLOCATE(diff_point_u(vv_mesh_new%np),diff_point_p(pp_mesh_new%np))
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
  ! SYM
     ALLOCATE(Hn1_sym   (H_mesh_new%np,  6,  m_max_c))
     ALLOCATE(Hn_sym    (H_mesh_new%np,  6,  m_max_c))
     ALLOCATE(phin1_sym (phi_mesh_new%np,2,  m_max_c))   
     ALLOCATE(phin_sym  (phi_mesh_new%np,2,  m_max_c))
     ALLOCATE(point_sym_H(H_mesh_new%np),point_sym_phi(phi_mesh_new%np))
     ALLOCATE(diff_point_H(H_mesh_new%np),diff_point_phi(phi_mesh_new%np))
  END IF
  !------------------------------------------------------------------------------

  !-----------INTERPOLATION FOR NAVIER-STOKES------------------------------------
 IF (type_pb=='nst' .OR. type_pb=='mhd' .OR. type_pb=='mxx') THEN
! CARO
     WRITE(*,*) 'list_mode_u avant= ', list_mode_u
     CALL read_restart_ns(vv_mesh, pp_mesh, time_u, list_mode_u, un, un_m1, pn, pn_m1, &
          incpn, incpn_m1, file_name, interpol=.true.)
     WRITE(*,*) 'list_mode_u apres= ', list_mode_u

     type_sym_u(1) = 1.d0
     type_sym_u(2) = 1.d0
     type_sym_u(3) =-1.d0
     type_sym_p(1) = 1.d0
     point_sym_u = 0
     point_sym_p = 0

  DO i=1, vv_mesh%np ! on parcourt la 1/2 boite et on trouve ds la gde boite le bon noeud
     point_sym_u(i) = find_point(vv_mesh_new,vv_mesh%rr(1,i),vv_mesh%rr(2,i))
     un_sym(point_sym_u(i), :, :)     = un(i, :, :)
     un_m1_sym(point_sym_u(i), :, :)  = un_m1(i, :, :)
  END DO

  DO i=1, pp_mesh%np
      point_sym_p(i) = find_point(pp_mesh_new,pp_mesh%rr(1,i),pp_mesh%rr(2,i))
      pn_sym(point_sym_p(i), :, :)       = pn(i, :, :)
      pn_m1_sym(point_sym_p(i), :, :)    = pn_m1(i, :, :)
      incpn_sym(point_sym_p(i), :, :)    = incpn(i, :, :)
      incpn_m1_sym(point_sym_p(i), :, :) = incpn_m1(i, :, :)
  END DO

!
!! ICI suppose ici que les points vv_mesh et vv_mesh_new coincident!
!  DO i= 1, vv_mesh%np
!    un_sym(i, :, :)    = un(i, :, :)
!    un_m1_sym(i, :, :) = un_m1(i, :, :)
!  END DO
!
!  DO i= 1, pp_mesh%np
!      pn_sym(i, :, :)       = pn(i, :, :)
!      pn_m1_sym(i, :, :)    = pn_m1(i, :, :)
!      incpn_sym(i, :, :)    = incpn(i, :, :)
!      incpn_m1_sym(i, :, :) = incpn_m1(i, :, :)
!  END DO
!
! ICI
     i_compteur = vv_mesh%np +1
  DO i=1, vv_mesh%np
     IF (vv_mesh_new%rr(2,point_sym_u(i)).EQ.0.d0) CYCLE
     point_sym_u(i_compteur) = find_point(vv_mesh_new, &
                                  vv_mesh_new%rr(1,point_sym_u(i)), &
                                 -1.d0*vv_mesh_new%rr(2,point_sym_u(i)))
       IF (i_compteur.GT.vv_mesh_new%np) THEN
           write(*,*) 'pb de compteur i_compteur_u= ', i_compteur
           STOP
       ENDIF
     un_sym(point_sym_u(i_compteur), 1:2, :)    = type_sym_u(1)* un_sym(point_sym_u(i), 1:2 ,:)
     un_sym(point_sym_u(i_compteur), 3:4, :)    = type_sym_u(2)* un_sym(point_sym_u(i), 3:4 ,:)
     un_sym(point_sym_u(i_compteur), 5:6, :)    = type_sym_u(3)* un_sym(point_sym_u(i), 5:6 ,:)
     un_m1_sym(point_sym_u(i_compteur), 1:2, :) = type_sym_u(1)* un_m1_sym(point_sym_u(i), 1:2 ,:)
     un_m1_sym(point_sym_u(i_compteur), 3:4, :) = type_sym_u(2)* un_m1_sym(point_sym_u(i), 3:4 ,:)
     un_m1_sym(point_sym_u(i_compteur), 5:6, :) = type_sym_u(3)* un_m1_sym(point_sym_u(i), 5:6 ,:)
! UPDATE
     i_compteur = i_compteur +1
  ENDDO

     i_compteur = pp_mesh%np +1
  DO i=1, pp_mesh%np
     IF (pp_mesh_new%rr(2,point_sym_p(i)).EQ.0.d0) CYCLE
     point_sym_p(i_compteur) = find_point(pp_mesh_new, &
                                  pp_mesh_new%rr(1,point_sym_p(i)), &
                                 -1.d0*pp_mesh_new%rr(2,point_sym_p(i)))
       IF (i_compteur.GT.pp_mesh_new%np) THEN
           write(*,*) 'pb de compteur i_compteur_p= ', i_compteur
           STOP
       ENDIF
     pn_sym(point_sym_p(i_compteur), 1:2, :)       = type_sym_p(1)* pn_sym(point_sym_p(i), 1:2 ,:)
     pn_m1_sym(point_sym_p(i_compteur), 1:2, :)    = type_sym_p(1)* pn_m1_sym(point_sym_p(i), 1:2 ,:)
     incpn_sym(point_sym_p(i_compteur), 1:2, :)    = type_sym_p(1)* incpn_sym(point_sym_p(i), 1:2 ,:)
     incpn_m1_sym(point_sym_p(i_compteur), 1:2, :) = type_sym_p(1)* incpn_m1_sym(point_sym_p(i), 1:2 ,:)
! UPDATE
     i_compteur = i_compteur +1
  ENDDO
! ICI
  DO i=1 , vv_mesh_new%np
     diff_point_u(i) = point_sym_u(i)-i
     IF (point_sym_u(i).eq.0) THEN
        write(*,*) 'pb de numerotation en i_u= ',i, vv_mesh_new%rr(1,i),vv_mesh_new%rr(2,i)
     ENDIF
  ENDDO
  IF (MAXVAL(IABS(diff_point_u)).ne.0) THEN
     write(*,*) 'pas la meme numerotation pour i_u MAXVAL= ', MAXVAL(IABS(diff_point_u))
  ENDIF
  DO i=1 , pp_mesh_new%np
     diff_point_p(i) = point_sym_p(i)-i
     IF (point_sym_p(i).eq.0) THEN
        write(*,*) 'pb de numerotation en i_p= ',i, pp_mesh_new%rr(1,i),pp_mesh_new%rr(2,i)
     ENDIF
  ENDDO
 IF (MAXVAL(IABS(diff_point_p)).ne.0) THEN
    write(*,*) 'pas la meme numerotation pour i_p MAXVAL= ', MAXVAL(IABS(diff_point_p))
 ENDIF
! ICI
!     CALL mesh_interp(pp_mesh_new,pp_mesh_new,pn_sym(:,1,1),pn_new(:,1,1),.true.)
!     DO m = 1, m_max_c
!        DO k = 1, 2
!           CALL mesh_interp(pp_mesh_new,pp_mesh_new,      pn_sym(:,k,m),      pn_new(:,k,m),.false.)
!           CALL mesh_interp(pp_mesh_new,pp_mesh_new,   pn_m1_sym(:,k,m),   pn_m1_new(:,k,m),.false.)
!           CALL mesh_interp(pp_mesh_new,pp_mesh_new,   incpn_sym(:,k,m),   incpn_new(:,k,m),.false.)
!           CALL mesh_interp(pp_mesh_new,pp_mesh_new,incpn_m1_sym(:,k,m),incpn_m1_new(:,k,m),.false.)
!        END DO
!     END DO
!
!     CALL mesh_interp(vv_mesh_new,vv_mesh_new,un_sym(:,1,1),un_new(:,1,1),.true.)
!     DO m = 1, m_max_c
!        DO k =1, 6
!           CALL mesh_interp(vv_mesh_new,vv_mesh_new,   un_sym(:,k,m),   un_new(:,k,m),.false.)
!           CALL mesh_interp(vv_mesh_new,vv_mesh_new,un_m1_sym(:,k,m),un_m1_new(:,k,m),.false.)
!        END DO
!     END DO
! ICI
    un_new       = un_sym
    un_m1_new    = un_m1_sym
    pn_new       = pn_sym
    pn_m1_new    = pn_m1_sym
    incpn_new    = incpn_sym
    incpn_m1_new = incpn_m1_sym
!
     CALL write_restart_ns(vv_mesh_new, pp_mesh_new, time_u, list_mode_u, &
          un_new, un_m1_new, pn_new, pn_m1_new, incpn_new, incpn_m1_new, file_name_new, 999, 1)
 END IF
  !------------------------------------------------------------------------------

  !-----------INTERPOLATION FOR MAXWELL------------------------------------------
 IF (type_pb=='mxw' .OR. type_pb=='mhd' .OR. type_pb=='mxx') THEN
! CARO
     WRITE(*,*) 'list_mode_H avant= ', list_mode_H
     CALL read_restart_maxwell(H_mesh, phi_mesh, time_h, list_mode_H, Hn, Hn1, &
          phin, phin1, file_name, interpol=.true.)
     WRITE(*,*) 'list_mode_H apres= ', list_mode_H

     type_sym_H(1) =-1.d0
     type_sym_H(2) =-1.d0
     type_sym_H(3) = 1.d0
     type_sym_phi(1) = -1.d0
     point_sym_H = 0
     point_sym_phi = 0

! eps=1.d-12, point_sym_H=0, type_sym_H(1:3), type_sym_phi(1)

  DO i=1, H_mesh%np ! on parcourt la 1/2 boite et on trouve ds la gde boite le bon noeud
     point_sym_H(i) = find_point(H_mesh_new,H_mesh%rr(1,i),H_mesh%rr(2,i))
     Hn_sym(point_sym_H(i), :, :)     = Hn(i, :, :)
     Hn1_sym(point_sym_H(i), :, :)    = Hn1(i, :, :)
  END DO

  DO i=1, phi_mesh%np
      point_sym_phi(i) = find_point(phi_mesh_new,phi_mesh%rr(1,i),phi_mesh%rr(2,i))
      phin_sym(point_sym_phi(i), :, :)     = phin(i, :, :)
      phin1_sym(point_sym_phi(i), :, :)    = phin1(i, :, :)
  END DO

! ICI
     i_compteur = H_mesh%np +1
  DO i=1, H_mesh%np
     IF (H_mesh_new%rr(2,point_sym_H(i)).EQ.0.d0) CYCLE
     point_sym_H(i_compteur) = find_point(H_mesh_new, &
                                  H_mesh_new%rr(1,point_sym_H(i)), &
                                 -1.d0*H_mesh_new%rr(2,point_sym_H(i)))
       IF (i_compteur.GT.H_mesh_new%np) THEN
           write(*,*) 'pb de compteur i_compteur_H= ', i_compteur
           STOP
       ENDIF
     Hn_sym(point_sym_H(i_compteur), 1:2, :)  = type_sym_H(1)* Hn_sym(point_sym_H(i), 1:2 ,:)
     Hn_sym(point_sym_H(i_compteur), 3:4, :)  = type_sym_H(2)* Hn_sym(point_sym_H(i), 3:4 ,:)
     Hn_sym(point_sym_H(i_compteur), 5:6, :)  = type_sym_H(3)* Hn_sym(point_sym_H(i), 5:6 ,:)
     Hn1_sym(point_sym_H(i_compteur), 1:2, :) = type_sym_H(1)* Hn1_sym(point_sym_H(i), 1:2 ,:)
     Hn1_sym(point_sym_H(i_compteur), 3:4, :) = type_sym_H(2)* Hn1_sym(point_sym_H(i), 3:4 ,:)
     Hn1_sym(point_sym_H(i_compteur), 5:6, :) = type_sym_H(3)* Hn1_sym(point_sym_H(i), 5:6 ,:)
! UPDATE
     i_compteur = i_compteur +1
  ENDDO

     i_compteur = phi_mesh%np +1
  DO i=1, phi_mesh%np
     IF (phi_mesh_new%rr(2,point_sym_phi(i)).EQ.0.d0) CYCLE
     point_sym_phi(i_compteur) = find_point(phi_mesh_new, &
                                  phi_mesh_new%rr(1,point_sym_phi(i)), &
                                 -1.d0*phi_mesh_new%rr(2,point_sym_phi(i)))
       IF (i_compteur.GT.phi_mesh_new%np) THEN
           write(*,*) 'pb de compteur i_compteur_phi= ', i_compteur
           STOP
       ENDIF
     phin_sym(point_sym_phi(i_compteur), 1:2, :)  = type_sym_phi(1)* phin_sym(point_sym_phi(i), 1:2 ,:)
     phin1_sym(point_sym_phi(i_compteur), 1:2, :) = type_sym_phi(1)* phin1_sym(point_sym_phi(i), 1:2 ,:)
! UPDATE
     i_compteur = i_compteur +1
  ENDDO
! ICI
  DO i=1, H_mesh_new%np
     diff_point_H(i) = point_sym_H(i)-i
     IF (point_sym_H(i).eq.0) THEN
        write(*,*) 'pb de numerotation en i_H= ',i, H_mesh_new%rr(1,i),H_mesh_new%rr(2,i)
     ENDIF
  ENDDO
      IF (MAXVAL(IABS(diff_point_H)).ne.0) THEN
         write(*,*) 'pas la meme numerotation pour i_H MAXVAL= ', MAXVAL(IABS(diff_point_H))
      ENDIF
  DO i= 1 , phi_mesh_new%np
     diff_point_phi(i) = point_sym_phi(i)-i
     IF (point_sym_phi(i).eq.0) THEN
        write(*,*) 'pb de numerotation en i_phi= ',i, phi_mesh_new%rr(1,i),phi_mesh_new%rr(2,i)
     ENDIF
  ENDDO
      IF (MAXVAL(IABS(diff_point_phi)).ne.0) THEN
         write(*,*) 'pas la meme numerotation pour i_phi MAXVAL= ', MAXVAL(IABS(diff_point_phi))
      ENDIF
! ICI
!
!     CALL mesh_interp(phi_mesh_new,phi_mesh_new,phin_sym(:,1,1),phin_new(:,1,1),.true.)
!     DO m = 1, m_max_c
!        DO k = 1, 2
!           CALL mesh_interp(phi_mesh_new,phi_mesh_new,    phin_sym(:,k,m),    phin_new(:,k,m),.false.)
!           CALL mesh_interp(phi_mesh_new,phi_mesh_new,   phin1_sym(:,k,m),   phin1_new(:,k,m),.false.)
!        END DO
!     END DO
!
!     CALL mesh_interp(H_mesh_new,H_mesh_new,Hn_sym(:,1,1),Hn_new(:,1,1),.true.)
!     DO m = 1, m_max_c
!        DO k =1, 6
!           CALL mesh_interp(H_mesh_new,H_mesh_new, Hn_sym(:,k,m), Hn_new(:,k,m),.false.)
!           CALL mesh_interp(H_mesh_new,H_mesh_new,Hn1_sym(:,k,m),Hn1_new(:,k,m),.false.)
!        END DO
!     END DO
! ICI
    Hn_new     = Hn_sym
    Hn1_new    = Hn1_sym
    phin_new   = phin_sym
    phin1_new  = phin1_sym
!
  CALL write_restart_maxwell(H_mesh_new, phi_mesh_new, time_h, list_mode_H, &
       Hn_new, Hn1_new, phin_new, phin1_new, file_name_new, 999, 1)
 END IF 
  !------------------------------------------------------------------------------

  !------------------------FIN------------------------------------------------
  CALL MPI_FINALIZE(code)
  !---------------------------------------------------------------------------

END PROGRAM interpol_sym

  FUNCTION find_point(mesh,r,z)  RESULT(n)
  
  USE def_type_mesh
  USE fem_tn_axi
  USE tn_parallele

  IMPLICIT NONE
  TYPE(mesh_type),                 INTENT(IN) :: mesh !type de maillage
  REAL(KIND=8),                    INTENT(IN) :: r, z
  INTEGER                                     :: n, i, j
  REAL(KIND=8)                                :: dist, dist_min
  
  dist_min = 1.d0
  DO i =1, mesh%np
     dist = sqrt((mesh%rr(1,i)-r)**2 + (mesh%rr(2,i)-z)**2)
     IF (dist .LE. dist_min) THEN
        dist_min = dist
        n = i
     ENDIF
  ENDDO 

  END FUNCTION find_point 

