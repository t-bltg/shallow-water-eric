PROGRAM mhd_prog

  USE def_type_mesh
  USE tn_parallele
  USE fem_tn_NS_MHD
  USE sub_plot
  USE initialisation
  USE boundary
  USE post_processing
  USE chaine_caractere
! mars 6
  USE fem_tn_axi
! mars 6
!TEST PRECESSION
  USE restart

  !---------------------------------------------------------------------------

  IMPLICIT NONE
  include 'mpif.h'

  !Champs pour Navier-Stokes-------------------------------------------------
  TYPE(mesh_type), POINTER                        :: pp_mesh, vv_mesh    
  REAL(KIND=8), POINTER, DIMENSION(:,:,:)         :: un, pn !Vitesse, pression
  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE     :: udif !un - r e_theta=vitesse ds ref wall
  ! (noeuds, type, mode)
  !---------------------------------------------------------------------------

  !Champs pour Maxwell--------------------------------------------------------
  TYPE(mesh_type), POINTER                        :: H_mesh, phi_mesh
  TYPE(interface_type), POINTER                   :: interface_H_mu, interface_H_phi
  REAL(KIND=8), POINTER,      DIMENSION(:,:,:)    :: Hn, phin, vel !Chmp mag, pot scalaire, vitesse
  REAL(KIND=8), POINTER,      DIMENSION(:)        :: sigma_field, mu_H_field
  REAL(KIND=8)                                    :: mu_phi
  !---------------------------------------------------------------------------

  !Variables de couplage------------------------------------------------------
  CHARACTER(len=3)                                :: type_pb
  !---------------------------------------------------------------------------

  !liste des modes------------------------------------------------------------
  INTEGER,      POINTER,      DIMENSION(:)        :: list_mode
  INTEGER,      POINTER,      DIMENSION(:)        :: select_mode_ns !modes that you want to be zero
  INTEGER,      POINTER,      DIMENSION(:)        :: select_mode_mxw !modes that you want to be zero
  !---------------------------------------------------------------------------

  !Noms reserves--------------------------------------------------------------
  INTEGER                                         :: nb_iteration, m_max_c 
  REAL(KIND=8)                                    :: dt, time, Re, Rem
  !---------------------------------------------------------------------------

  !Gestion locale-------------------------------------------------------------
  INTEGER                                         :: code, rang, nb_procs
  INTEGER                                         :: it, k, i, mode, unit, j
  REAL(KIND=8)                                    :: user_time, tps, dummy
  REAL(KIND=8)                                    :: err, norm, normr, normt, normz, pol, tore
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:)       :: vect
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:)         :: vel_loc
  !-------------END OF DECLARATIONS-------------------------------------------

  !-------------POST_PROCESSING-----------------------------------------------
  INTEGER                                         :: freq_restart, freq_en, freq_plot
  INTEGER                                         :: select_ns, select_mxw
  !Declare whatever you want here for your own purpose
! mars 6
  !Taylor-Couette------------------------------------------------------------
  LOGICAL                                         :: if_select 
  LOGICAL                                         :: START_NL=.true.
  !Traitement des fichiers de sortie-----------------------------------------
  CHARACTER(LEN=2)                                :: truc
  INTEGER                                         :: dd_end, ff_end
  CHARACTER(LEN=100)                              :: en_file
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:)         :: e_c_u, e_cm_h
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:)         :: e_c_u_sym, e_c_u_anti, e_cm_h_sym, e_cm_h_anti
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:)         :: e_c_u_dif, e_c_u_sym_dif, e_c_u_anti_dif
  REAL(KIND=8),              DIMENSION(3)         :: x_anemo_v, x_anemo_h
  REAL(KIND=8),              DIMENSION(5)         :: y_anemo_v, y_anemo_h
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:)       :: norme_u, norme_h
  INTEGER , ALLOCATABLE, DIMENSION(:)             :: anemo_v, anemo_h
  CHARACTER(LEN=100)                              :: e_c_u_file, e_cm_h_file, e_c_u_z_file
  CHARACTER(LEN=100)                              :: norme_u_file, norme_h_file
  CHARACTER(LEN=100)                              :: anemometre_v_file, anemometre_h_file
  REAL(KIND=8),DIMENSION(3)                       :: type_sym_u, type_sym_h
  INTEGER , ALLOCATABLE, DIMENSION(:)             :: point_sym_u, point_sym_h
  REAL(KIND=8), DIMENSION(3)                      :: dipole, vort_moyenne
  REAL(KIND=8), DIMENSION(3,3)                    :: quadripole
! mars 6
! cas PRECESSION
  CHARACTER(LEN=64)                                :: file_vel
  REAL(KIND=8), DIMENSION(:),   ALLOCATABLE  :: x

  !---------------------------------------------------------------------------
  
  !-------------DEBUT PARALLELISATION-----------------------------------------
  CALL MPI_INIT(code)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)
  !---------------------------------------------------------------------------

  !-------------INITIALISATION------------------------------------------------
  CALL initial(vv_mesh, pp_mesh, H_mesh, phi_mesh, interface_H_phi, interface_H_mu, list_mode, &
                     un, pn, Hn, phin, vel, mu_H_field, sigma_field, mu_phi, &
                     Re, Rem, time, dt, nb_iteration, m_max_c, type_pb, test_de_convergence)
  !---------------------------------------------------------------------------

  !-------------YOUR OWN INITIALIZATIONS--------------------------------------
  OPEN(UNIT = 21, FILE = 'data', FORM = 'formatted', STATUS = 'unknown')
  CALL read_until(21, 'data_postproc')
  READ(21,*) freq_restart, freq_en, freq_plot
  READ(21,*) x_anemo_v(:)
  READ(21,*) y_anemo_v(:)
  READ(21,*) x_anemo_h(:)
  READ(21,*) y_anemo_h(:)
  CLOSE(21) 
  !OPEN(UNIT = 21, FILE = 'data', FORM = 'formatted', STATUS = 'unknown')
  !CALL read_until(21, 'my_stuff')
  !READ(21,*) stuff_1, stuff_2, whatever 
  !CLOSE(21)
!  IF (type_pb == 'mxx') THEN
!     ALLOCATE(vel_loc(H_mesh%np))
!     WRITE(*,*) 'on fait MXX'
!     !Vous gerez vel ici
!     DO i = 1, m_max_c
!        !List_mode contient les modes physiques
!        IF(list_mode(i)/=0) THEN
!           vel(:,:,i) = 0.d0
!        ENDIF
!     END DO
!     !IF (rang==0) CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, vel(:,3,1), 'vel3.plt')
!  END IF 
  !---------------------------------------------------------------------------
!----------------RESTART AVEC MULTIPLICATION ENERGIE MODES H IMPAIRS POUR Taylor-Couette----
! avril 11
        IF (START_NL) THEN
!       ALLOCATE(x(vv_mesh%np))
!       WRITE(*,*) 'bruit supplementaire'

         DO i = 1, m_max_c
          IF (list_mode(i).LE.5) THEN
           Hn(:,:,i)    = 300.d0*Hn(:,:,i)
           phin(:,:,i)  = 300.d0*phin(:,:,i)
          ENDIF
!! le 21/05/09 pas de select mode
!!         IF (mod((list_mode(i)+1),2)==0) THEN
!            Hn(:,:,i)    = 10.d0*Hn(:,:,i)
!            phin(:,:,i)  = 10.d0*phin(:,:,i)
!!         ENDIF
!! bruit HYDRO PRECESSION
!       IF (list_mode(i).LE.3) THEN
!          CALL RANDOM_SEED
!          CALL RANDOM_NUMBER(x)
!!         un(:,3,i)= 1.d-3*(x - 0.5d0) !les 28 et 29 juin 2010
!          un(:,3,i)=un(:,3,i)+ 1.d-3*(x - 0.5d0)
!       END IF

          ENDDO

        END IF
! avril 11
!----------------RESTART AVEC MULTIPLICATION ENERGIE MODES H IMPAIRS POUR Taylor-Couette---

  !----------------------QUELQUES VERIFICATIONS------------------------------
  ! nst = Navier-Stokes
  ! mxw = Maxwell
  ! mxx = Maxwell + u_restart =.true. (defined in initialisation)
  ! mhd = MHD  
  IF (type_pb == 'nst' .OR. type_pb == 'mhd') THEN 
     !IF (rang==0) CALL plot_scalar_field(pp_mesh%jj,pp_mesh%rr, pn(:,1,1),'pn_init.plt')
     !IF (rang==0) CALL plot_scalar_field(vv_mesh%jj,vv_mesh%rr, un(:,2,1),'un_init.plt')
     norm = norme_H1_champ_par(vv_mesh,list_mode,un)+ 1.d-15
     err = norme_div_par(vv_mesh,list_mode, un)
     IF (rang==0) WRITE(*,*) 'div(un) relative champ init   = ', err/norm
     IF (rang==0) write(*,*) 'div(un) absolue  champ init   = ', err
     CALL trace_profile(vv_mesh, un, 0, freq_plot, list_mode, 'Vin')
  END IF

  IF (type_pb == 'mxw' .OR. type_pb == 'mhd' .OR. type_pb == 'mxx') THEN

   DO i = 1, m_max_c
     CALL norme_interface(H_mesh,phi_mesh,interface_H_phi,mu_H_field,mu_phi,Hn(:,:,i),phin(:,:,i),list_mode(i),err)
     WRITE(*,*) ' erreur on the interface H_phi for initial data', err, 'for the mode ',list_mode(i)
     CALL norme_interface_H_mu(H_mesh,interface_H_mu,mu_H_field,Hn(:,:,i),list_mode(i),err)
     WRITE(*,*) ' erreur on the interface H_mu initial data', err, 'for the mode ',list_mode(i)
   ENDDO
! CARO a remettre
     CALL trace_profile(H_mesh, Hn, 0, freq_plot, list_mode, 'Hin')
! cas PRECESSION mxw pr tracer le champ de vitesse sur maillage en H
!    CALL trace_profile(H_mesh, vel, 0, freq_plot, list_mode, 'Vci')
!    file_vel = 'vel_instationnaire'
!    CALL write_restart_maxwell(H_mesh, phi_mesh, time, list_mode, vel, vel, phin, phin, file_vel, 0, freq_restart)
! cas PRECESSION mxw pr tracer le champ de vitesse sur maillage en H---- fin
     !IF (rang==0)  CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn(:,1,1),'Hn1_init.plt')
     !IF (rang==0)  CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn(:,2,1),'Hn2_init.plt')
     !IF (rang==0)  CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn(:,3,1),'Hn3_init.plt')
     !IF (rang==0)  CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn(:,4,1),'Hn4_init.plt')
     !IF (rang==0)  CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn(:,5,1),'Hn5_init.plt')
     !IF (rang==0)  CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn(:,6,1),'Hn6_init.plt')
     !IF (rang==0)  CALL plot_scalar_field(phi_mesh%jj, phi_mesh%rr, phin(:,1,1),'phin1_init.plt')
     !IF (rang==0)  CALL plot_scalar_field(phi_mesh%jj, phi_mesh%rr, phin(:,2,1),'phin2_init.plt')
     norm = norme_H1_champ_par(H_mesh,list_mode,Hn)+ 1.d-15
     err = norme_div_par(H_mesh,list_mode, Hn)
     WRITE(*,*) 'div(Hn) relative champ init   = ', err/norm
     write(*,*) 'div(Hn) absolue  champ init   = ', err
  END IF
  !---------------------------------------------------------------------------


! mars 6
!**********************OUVERTURE DES FICHIERS DE SORTIE********************************
    !------------FICHIERS DE SORTIE NAVIER-STOKES-----------------------------------

    IF (type_pb=='nst' .OR. type_pb=='mhd') THEN
       !energie
! cas PRECESSION 28/07/09
       !=======
!    type_sym_u(1) = 1.d0
!    type_sym_u(2) = 1.d0
!    type_sym_u(3) =-1.d0
!    ALLOCATE(point_sym_u(vv_mesh%np))
!    CALL calcul_point_sym(vv_mesh, point_sym_u)
! cas PRECESSION 28/07/09
  !------------------------REF WALL-----------------------------------
  ALLOCATE(udif(vv_mesh%np,6,m_max_c))

       en_file = 'e_c_u'
       dd_end = last_c_leng (50, en_file)
       DO i = 1, m_max_c
          dd_end = last_c_leng (50, en_file)
          WRITE(truc,'(i2)') list_mode(i)
          ff_end = start_of_string (truc)
          en_file = en_file(1:dd_end)//'_'//truc(ff_end:)
         DO k = 1, 6
            udif(:,k,i)= un(:,k,i) - vv_exact_rot_axez(k,vv_mesh%rr,list_mode(i),time)
         END DO
       END DO
       e_c_u_file = en_file
       OPEN(UNIT=20,FILE=e_c_u_file, FORM='formatted', STATUS='unknown')
       WRITE(20,*)'#energie par mode'
       WRITE(20,*)'# Temps ec ordre croissant mode -> ecu ecusym ecuanti ecu_d ecus_d ecua_d'
       ALLOCATE(e_c_u(m_max_c),e_c_u_sym(m_max_c),e_c_u_anti(m_max_c))
       ALLOCATE(e_c_u_dif(m_max_c),e_c_u_sym_dif(m_max_c),e_c_u_anti_dif(m_max_c))
!      CALL val_ener(vv_mesh, list_mode, un, e_c_u)
!      CALL val_ener_sym_glob(vv_mesh, list_mode, un, point_sym_u, e_c_u, e_c_u_sym, e_c_u_anti, type_sym_u)
! cas PRECESSION 28/07/09
       CALL val_ener_sym_centrale(vv_mesh, list_mode, un, e_c_u, e_c_u_sym, e_c_u_anti)
       CALL val_ener_sym_centrale(vv_mesh, list_mode, udif, e_c_u_dif, &
                                  e_c_u_sym_dif, e_c_u_anti_dif)
           WRITE(20,103) time, e_c_u, e_c_u_sym, e_c_u_anti, e_c_u_dif,&
                               e_c_u_sym_dif, e_c_u_anti_dif
       CLOSE(20)
       !energie axiale
       !==========
       en_file = 'e_c_u_z'
       dd_end = last_c_leng (50, en_file)
       DO i = 1, m_max_c
          dd_end = last_c_leng (50, en_file)
          WRITE(truc,'(i2)') list_mode(i)
          ff_end = start_of_string (truc)
          en_file = en_file(1:dd_end)//'_'//truc(ff_end:)
       END DO
       e_c_u_z_file = en_file
       OPEN(UNIT=60,FILE=e_c_u_z_file, FORM='formatted', STATUS='unknown')
       WRITE(60,*)'#energie axiale par mode'
       WRITE(60,*)'# Temps    valeur par ordre croissant de mode -> '
       CALL val_ener(vv_mesh, list_mode, un(:,5:6,:), e_c_u)
           WRITE(60,103) time, e_c_u
       CLOSE(60)
       !composante
       !==========
       en_file = 'norme_comp_u'
       dd_end = last_c_leng (50, en_file)
       DO i = 1, m_max_c
          dd_end = last_c_leng (50, en_file)
          WRITE(truc,'(i2)') list_mode(i)
          ff_end = start_of_string (truc)
          en_file = en_file(1:dd_end)//'_'//truc(ff_end:)
       END DO
       norme_u_file = en_file
       OPEN(UNIT=25,FILE=norme_u_file, FORM='formatted', STATUS='unknown')
       WRITE(25,*)'#normes des composantes par mode'
       WRITE(25,*)'# Temps    C1(m=0..m_nax_c) C2 C3 C4 C5 C6  '
       ALLOCATE(norme_u(m_max_c,6))
          DO i = 1, m_max_c
             DO j = 1, 6
                 CALL ns_0(vv_mesh, un(:,j,i), norme_u(i,j))
             ENDDO
          ENDDO
          WRITE(25,103) time, norme_u(:,1), norme_u(:,2), norme_u(:,3),& 
                              norme_u(:,4), norme_u(:,5), norme_u(:,6)
       CLOSE(25)
    
    ENDIF

    !------------FICHIERS DE SORTIE MAXWELL-----------------------------------

    IF (type_pb=='mxw' .OR. type_pb=='mhd' .OR. type_pb=='mxx') THEN
       !energie
       !=======
! cas PRECESSION 28/07/09
!    type_sym_h(1) = 1.d0 ! pr G=2; different pr G=2pi :-1.d0
!    type_sym_h(2) = 1.d0 ! pr G=2; different pr G=2pi :-1.d0
!    type_sym_h(3) =-1.d0 ! pr G=2; different pr G=2pi :1.d0
!    ALLOCATE(point_sym_h(H_mesh%np))
!    CALL calcul_point_sym(H_mesh, point_sym_h)
! cas PRECESSION 28/07/09

       en_file = 'e_cm_h'
       dd_end = last_c_leng (50, en_file)
       DO i = 1, m_max_c
          dd_end = last_c_leng (50, en_file)
          WRITE(truc,'(i2)') list_mode(i)
          ff_end = start_of_string (truc)
          en_file = en_file(1:dd_end)//'_'//truc(ff_end:)
       END DO
       e_cm_h_file = en_file
       OPEN(UNIT=50,FILE=e_cm_h_file, FORM='formatted', STATUS='unknown')
       WRITE(50,*)'#energie par mode dans le conducteur '
       WRITE(50,*)'# Temps    valeur par ordre croissant de mode -> '
       ALLOCATE(e_cm_h(m_max_c),e_cm_h_sym(m_max_c),e_cm_h_anti(m_max_c))
!      CALL val_ener(H_mesh, list_mode, Hn, e_cm_h)
!      CALL val_ener_sym_glob(H_mesh, list_mode, Hn, point_sym_h, e_cm_h, e_cm_h_sym, e_cm_h_anti, type_sym_h)
! cas PRECESSION 28/07/09
       CALL val_ener_sym_centrale(H_mesh, list_mode, Hn, e_cm_h, e_cm_h_sym, e_cm_h_anti)
           WRITE(50,103) time, e_cm_h, e_cm_h_sym, e_cm_h_anti
       CLOSE(50)
       !composante
       !==========
       en_file = 'norme_comp_h'
       dd_end = last_c_leng (50, en_file)
       DO i = 1, m_max_c
          dd_end = last_c_leng (50, en_file)
          WRITE(truc,'(i2)') list_mode(i)
          ff_end = start_of_string (truc)
          en_file = en_file(1:dd_end)//'_'//truc(ff_end:)
       END DO
       norme_h_file = en_file
       OPEN(UNIT=55,FILE=norme_h_file, FORM='formatted', STATUS='unknown')
       WRITE(55,*)'#normes des composantes par mode pour le champ magnetique dans le conducteur'
       WRITE(55,*)'# Temps    C1(m=0..m_nax_c) C2 C3 C4 C5 C6  '
       ALLOCATE(norme_h(m_max_c,6))
          DO i = 1, m_max_c
             DO j = 1, 6
                CALL ns_0(H_mesh, Hn(:,j,i), norme_h(i,j))
             ENDDO
          ENDDO
          WRITE(55,103) time, norme_h(:,1), norme_h(:,2), norme_h(:,3),&
                              norme_h(:,4), norme_h(:,5), norme_h(:,6)
       CLOSE(55)

   
    !Dipole et quadrupole
    !==========
       en_file = 'dip_quad_h'
       OPEN(UNIT=58,FILE=en_file, FORM='formatted', STATUS='unknown')
       WRITE(58,*)'# dipole et quadrupole ds espace phys pr le champ magnetique dans le conducteur'
       WRITE(58,*)'# Temps dip(1:3) quad(3:3)'
         CALL moments(H_mesh, list_mode, Hn, dipole, quadripole)
         WRITE(58,103) time, dipole(:), &
                       quadripole(1,:), quadripole(2,:), quadripole(3,:)
       CLOSE(58)


    ENDIF

    !Anemometre
    !==========
 
    !Navier-Stokes
    ALLOCATE(anemo_v(15))
    IF (type_pb=='nst' .OR. type_pb=='mhd') THEN   

       
       anemo_v(1)  = find_point(vv_mesh,x_anemo_v(1),y_anemo_v(1)) 
       anemo_v(2)  = find_point(vv_mesh,x_anemo_v(2),y_anemo_v(1))  
       anemo_v(3)  = find_point(vv_mesh,x_anemo_v(3),y_anemo_v(1))
       anemo_v(4)  = find_point(vv_mesh,x_anemo_v(1),y_anemo_v(2))
       anemo_v(5)  = find_point(vv_mesh,x_anemo_v(2),y_anemo_v(2))
       anemo_v(6)  = find_point(vv_mesh,x_anemo_v(3),y_anemo_v(2))
       anemo_v(7)  = find_point(vv_mesh,x_anemo_v(1),y_anemo_v(3))    
       anemo_v(8)  = find_point(vv_mesh,x_anemo_v(2),y_anemo_v(3))    
       anemo_v(9)  = find_point(vv_mesh,x_anemo_v(3),y_anemo_v(3))
       anemo_v(10) = find_point(vv_mesh,x_anemo_v(1),y_anemo_v(4))
       anemo_v(11) = find_point(vv_mesh,x_anemo_v(2),y_anemo_v(4))
       anemo_v(12) = find_point(vv_mesh,x_anemo_v(3),y_anemo_v(4))
       anemo_v(13) = find_point(vv_mesh,x_anemo_v(1),y_anemo_v(5))  
       anemo_v(14) = find_point(vv_mesh,x_anemo_v(2),y_anemo_v(5))  
       anemo_v(15) = find_point(vv_mesh,x_anemo_v(3),y_anemo_v(5))

       en_file = 'anemometre_V'
       dd_end = last_c_leng (50, en_file)
       DO i = 1, m_max_c
          dd_end = last_c_leng (50, en_file)
          WRITE(truc,'(i2)') list_mode(i)
          ff_end = start_of_string (truc)
          en_file = en_file(1:dd_end)//'_'//truc(ff_end:)
       END DO
       anemometre_v_file = en_file
       OPEN(UNIT=56,FILE=anemometre_v_file, FORM='formatted', STATUS='unknown')
       WRITE(56,*)'#Anemometre pour le champ de vitesse'
       WRITE(56,*)'#Coordonnees des points :    '
       DO i = 1, 15
           WRITE(56,*) '# point ',i, ' : r=', vv_mesh%rr(1,anemo_v(i)), '; z = ', vv_mesh%rr(2,anemo_v(i)) 
       ENDDO 
       CLOSE(56)
    ENDIF

    !Maxwell
    ALLOCATE(anemo_h(15))
    IF (type_pb=='mxw' .OR. type_pb=='mhd' .OR. type_pb=='mxx') THEN
       anemo_h(1)  = find_point(H_mesh,x_anemo_h(1),y_anemo_h(1)) 
       anemo_h(2)  = find_point(H_mesh,x_anemo_h(2),y_anemo_h(1))
       anemo_h(3)  = find_point(H_mesh,x_anemo_h(3),y_anemo_h(1))
       anemo_h(4)  = find_point(H_mesh,x_anemo_h(1),y_anemo_h(2))
       anemo_h(5)  = find_point(H_mesh,x_anemo_h(2),y_anemo_h(2))
       anemo_h(6)  = find_point(H_mesh,x_anemo_h(3),y_anemo_h(2))
       anemo_h(7)  = find_point(H_mesh,x_anemo_h(1),y_anemo_h(3))
       anemo_h(8)  = find_point(H_mesh,x_anemo_h(2),y_anemo_h(3))
       anemo_h(9)  = find_point(H_mesh,x_anemo_h(3),y_anemo_h(3))
       anemo_h(10) = find_point(H_mesh,x_anemo_h(1),y_anemo_h(4))
       anemo_h(11) = find_point(H_mesh,x_anemo_h(2),y_anemo_h(4))
       anemo_h(12) = find_point(H_mesh,x_anemo_h(3),y_anemo_h(4))
       anemo_h(13) = find_point(H_mesh,x_anemo_h(1),y_anemo_h(5))
       anemo_h(14) = find_point(H_mesh,x_anemo_h(2),y_anemo_h(5))
       anemo_h(15) = find_point(H_mesh,x_anemo_h(3),y_anemo_h(5))

       en_file = 'anemometre_H'
       dd_end = last_c_leng (50, en_file)
       DO i = 1, m_max_c
          dd_end = last_c_leng (50, en_file)
          WRITE(truc,'(i2)') list_mode(i)
          ff_end = start_of_string (truc)
          en_file = en_file(1:dd_end)//'_'//truc(ff_end:)
       END DO
       anemometre_h_file = en_file
       OPEN(UNIT=57,FILE=anemometre_h_file, FORM='formatted', STATUS='unknown')
       WRITE(57,*) '#Anemometre pour le champ magnetique'
       WRITE(57,*) '#Coordonnees des points :    '
       DO i = 1, 15
           WRITE(57,*) '# point ',i, 'r=', H_mesh%rr(1,anemo_h(i)), '; z = ', H_mesh%rr(2,anemo_h(i))
       ENDDO
       CLOSE(57)
    ENDIF

    103 FORMAT(1500(e22.9,2x))
    !------------------------------------------------------------------------------
! mars 6
    !Gestion des modes nuls pour mhd
     OPEN(UNIT = 21, FILE = 'data', FORM = 'formatted', STATUS = 'unknown')
     CALL read_until(21, 'data_select_mode_nuls')
     READ(21,*) if_select
     CLOSE(21)
!----------------MISE A ZERO POUR Taylor-Couette----------------------------
     IF (if_select) THEN
        select_ns = 0
        select_mxw = 0
        DO i = 1, m_max_c
          IF (mod((list_mode(i)+1),2)==0) THEN
             select_ns = select_ns + 1
          ELSE
             select_mxw = select_mxw + 1
          ENDIF
        END DO
        ALLOCATE(select_mode_ns(select_ns),select_mode_mxw(select_mxw))
        select_ns = 0
        select_mxw = 0
        DO i = 1, m_max_c
          IF (mod((list_mode(i)+1),2)==0) THEN
             select_ns = select_ns + 1
             select_mode_ns(select_ns) = i
          ELSE
             select_mxw = select_mxw + 1
             select_mode_mxw(select_mxw) = i
          ENDIF
        END DO
     END IF 
!----------------MISE A ZERO POUR Taylor-Couette----------------------------

     
    !Gestion des modes nuls pour mhd
!********************** FIN OUVERTURE DES FICHIERS DE SORTIE********************************

  !------------------------RESOLUTION ----------------------------------------
  tps = user_time(dummy)
  DO it = 1, nb_iteration  ! DEBUT RESOLUTION

     time = time + dt    
     WRITE(*,*) 'Iteration = ',it,' ;  Time = ', time

     !======Navier Stokes  
     IF (type_pb == 'nst') THEN

        CALL navier_stokes(time, un, pn)

        !TEST pour Navier_stokes
        norm= norme_H1_champ_par(vv_mesh,list_mode,un)+1.d-15
        err = norme_div_par(vv_mesh,list_mode, un)
        IF(rang==0) WRITE(*,*) 'div(un) relative    = ', err/norm
        IF(rang==0) WRITE(*,*) 'div(un) absolue     = ', err
        DO i = 1, m_max_c
           unit = rang*m_max_c+10+i
           WRITE(unit,*) time, &
                    norme_l2_champ(vv_mesh, list_mode(i:i), un(:,:,i:i)),err/norm
        END DO
        DO mode = 1, SIZE(list_mode)
           IF (list_mode(mode)==0) THEN
               CALL nv_0_cn(vv_mesh, un(:,:,mode), pol, tore)
               IF (tore/=0) WRITE(53,*) time, pol, tore
           END IF
        END DO
        normr = norme_L2_champ_par(vv_mesh,list_mode,un(:,1:2,:))
        normt = norme_L2_champ_par(vv_mesh,list_mode,un(:,3:4,:))
        normz = norme_L2_champ_par(vv_mesh,list_mode,un(:,5:6,:))
        IF (rang==0) THEN
           WRITE(52,*) time, SQRT((normr**2+normz**2)/normt**2) !Poloidal/Toroidal, Energy
        END IF
        norm = norme_L2_champ(vv_mesh, list_mode(1:1), un(:,:,1:1))
        normt = norme_L2_champ_par(vv_mesh,list_mode,un)
        IF (normt > 1.d5) THEN
           WRITE(*,*) ' EXPLOSION '
           STOP
        END IF
        IF (list_mode(1)==0) WRITE(51,*) time, norm, normt ! Mode 0, all modes, Energy

     ELSE IF (type_pb == 'mxw' .OR. type_pb == 'mxx') THEN

       CALL maxwell(time, Hn, phin)

        !TEST pour Maxwell
        norm = norme_H1_champ_par(H_mesh,list_mode,Hn)+1.d-15
        err = norme_div_par(H_mesh,list_mode, Hn)
        IF(rang==0) WRITE(*,*) 'div(Hn) relative    = ', err/norm
        IF(rang==0) WRITE(*,*) 'div(Hn) absolue     = ', err
        DO i = 1, m_max_c
           unit = 100 + rang*m_max_c + i
           WRITE(unit,*) time, &
                     norme_l2_champ(H_mesh, list_mode(i:i), Hn(:,:,i:i)),err/norm
        END DO


     ELSE IF (type_pb == 'mhd') THEN
        IF (if_select) THEN
            CALL mhd(time, un, pn, Hn, phin, select_mode_ns, select_mode_mxw)
        ELSE
            CALL mhd(time, un, pn, Hn, phin)
        ENDIF
        !TEST pour Navier_stokes
        normr = norme_L2_champ_par(vv_mesh,list_mode,un(:,1:2,:))
        normt = norme_L2_champ_par(vv_mesh,list_mode,un(:,3:4,:))
        normz = norme_L2_champ_par(vv_mesh,list_mode,un(:,5:6,:))
        DO mode = 1, SIZE(list_mode)
           IF (list_mode(mode)==0) THEN
               CALL nv_0_cn(vv_mesh, un(:,:,mode), pol, tore)
               IF (tore/=0) WRITE(53,*) time, pol , tore
           END IF
        END DO
        IF (rang==0) THEN
           WRITE(52,*) time, SQRT((normr**2+normz**2)/normt**2) !Poloidal/Toroidal, Energy
        END IF
        norm = norme_L2_champ(vv_mesh, list_mode(1:1), un(:,:,1:1))
        normt = norme_L2_champ_par(vv_mesh,list_mode,un)
        IF (normt > 1.d5) THEN
           WRITE(*,*) ' EXPLOSION '
           STOP
        END IF
        IF (list_mode(1)==0) WRITE(51,*) time, norm, normt ! Mode 0, all modes, Energy

        norm = norme_H1_champ_par(vv_mesh,list_mode,un)+1.d-15
        err = norme_div_par(vv_mesh,list_mode, un)
        IF(rang==0) WRITE(*,*) 'div(un) relative    = ', err/norm
        IF(rang==0) WRITE(*,*) 'div(un) absolue     = ', err
        DO i = 1, m_max_c
           unit = rang*m_max_c+10+i
           WRITE(unit,*) time, &
                norme_l2_champ(vv_mesh, list_mode(i:i), un(:,:,i:i)),err/norm
        END DO

        !TEST pour Maxwell
        norm = norme_H1_champ_par(H_mesh,list_mode,Hn)+1.d-15
        err = norme_div_par(H_mesh,list_mode, Hn)
        IF(rang==0) WRITE(*,*) 'div(Hn) relative    = ', err/norm
        IF(rang==0) WRITE(*,*) 'div(Hn) absolue     = ', err
        DO i = 1, m_max_c
           unit = 100 + rang*m_max_c + i
           WRITE(unit,*) time, &
                norme_l2_champ(H_mesh, list_mode(i:i), Hn(:,:,i:i)),err/norm
        END DO

     END IF

!    !===================Put your postprocessing stuff here
!     ! Whatever 
!     IF (mod(it,freq_en) == 0) THEN
!        !CALL your_stuff_here
!     END IF
!     !Traces au format PlotMTV
!     IF (mod(it,freq_plot) == 0) THEN
!        !CALL trace_whatever_you_want_here
!     ENDIF
!     !Ecriture fichier restart de securite
!     IF (mod(it, freq_restart) == 0) THEN
!        CALL  sauvegarde(it,freq_restart)
!     ENDIF

!**************************************************************************************
! mars 6
    IF (type_pb=='nst' .OR. type_pb=='mhd') THEN

        !Ecriture des fichiers de sortie
        !===============================
        IF (mod(it,freq_en) == 0) THEN
           !Energie
          DO i = 1, m_max_c
             DO k = 1, 6
                udif(:,k,i)= un(:,k,i) - vv_exact_rot_axez(k,vv_mesh%rr,list_mode(i),time)
             END DO
          END DO

           OPEN(UNIT=20,FILE=e_c_u_file, FORM='formatted', POSITION = 'append', &
           STATUS='unknown')
!          CALL val_ener(vv_mesh, list_mode, un, e_c_u)
!          CALL val_ener_sym_glob(vv_mesh, list_mode, un, point_sym_u, e_c_u, e_c_u_sym, e_c_u_anti, type_sym_u)
! cas PRECESSION 28/07/09
           CALL val_ener_sym_centrale(vv_mesh, list_mode, un, e_c_u, e_c_u_sym, e_c_u_anti)
           CALL val_ener_sym_centrale(vv_mesh, list_mode, udif, e_c_u_dif, &
                                  e_c_u_sym_dif, e_c_u_anti_dif)
           WRITE(20,103) time, e_c_u, e_c_u_sym, e_c_u_anti, e_c_u_dif,&
                               e_c_u_sym_dif, e_c_u_anti_dif
           CLOSE(20)
           !energie axiale
           OPEN(UNIT=60,FILE=e_c_u_z_file, FORM='formatted', POSITION = 'append', &
           STATUS='unknown')
           CALL val_ener(vv_mesh, list_mode, un(:,5:6,:), e_c_u)
               WRITE(60,103) time, e_c_u
           CLOSE(60)
           !Composantes
           OPEN(UNIT=25,FILE=norme_u_file, FORM='formatted',POSITION = 'append', &
            STATUS='unknown')
           DO i = 1, m_max_c
              DO j = 1, 6
                 CALL ns_0(vv_mesh, un(:,j,i), norme_u(i,j))
              ENDDO
           ENDDO
           WRITE(25,103) time, norme_u(:,1), norme_u(:,2), norme_u(:,3),&
                              norme_u(:,4), norme_u(:,5), norme_u(:,6)
           CLOSE(25)
           !Anemometre 
           OPEN(UNIT=56,FILE=anemometre_v_file, FORM='formatted', POSITION = 'append', &
           STATUS='unknown')
           WRITE(56,103) time, un(anemo_v(:),1,:), un(anemo_v(:),2,:), un(anemo_v(:),3,:), un(anemo_v(:),4,:), &
                               un(anemo_v(:),5,:), un(anemo_v(:),6,:)
         
           CLOSE(56)
        ENDIF
        !Traces au format PlotMTV
        IF (mod(it,freq_plot) == 0) THEN
           CALL trace_profile(vv_mesh, un, it, freq_plot, list_mode, 'Uns')   
           CALL trace_profile(vv_mesh, udif, it, freq_plot, list_mode, 'Udi')   
        ENDIF
        !Ecriture fichier restart de securite
        IF (mod(it,freq_restart) == 0) THEN
           CALL sauvegarde(it,freq_restart)
! sauvegarde ds referentiel du cylindre
           DO i = 1, m_max_c
              DO k = 1, 6
                 udif(:,k,i)= un(:,k,i) - vv_exact_rot_axez(k,vv_mesh%rr,list_mode(i),time)
              END DO
           END DO
           file_vel = 'u_dif_cyl'
           CALL write_restart_ns(vv_mesh, pp_mesh, time, list_mode, udif, udif, pn, pn, &
                                 pn, pn, file_vel, it, freq_restart)
        ENDIF
    ENDIF
!**************************************************************************************
    IF (type_pb=='mxw' .OR. type_pb=='mhd' .OR. type_pb=='mxx') THEN

        !Ecriture des fichiers de sortie
        !===============================
        IF (mod(it,freq_en) == 0) THEN
           !Energie
           OPEN(UNIT=50,FILE=e_cm_h_file, FORM='formatted',POSITION = 'append',&
           STATUS='unknown')
!          CALL val_ener(H_mesh, list_mode, Hn, e_cm_h)
!          CALL val_ener_sym_glob(H_mesh, list_mode, Hn, point_sym_h, e_cm_h, e_cm_h_sym, e_cm_h_anti, type_sym_h)
           CALL val_ener_sym_centrale(H_mesh, list_mode, Hn, e_cm_h, e_cm_h_sym, e_cm_h_anti)
           WRITE(50,103) time, e_cm_h, e_cm_h_sym, e_cm_h_anti
           CLOSE(50)
           !Composantes
           OPEN(UNIT=55,FILE=norme_h_file, FORM='formatted',POSITION = 'append', &
           STATUS='unknown')
           DO i = 1, m_max_c
              DO j = 1, 6
                 CALL ns_0(H_mesh, Hn(:,j,i), norme_h(i,j))
              ENDDO
           ENDDO
           WRITE(55,103) time, norme_h(:,1), norme_h(:,2), norme_h(:,3),&
                              norme_h(:,4), norme_h(:,5), norme_h(:,6)
           CLOSE(55)
           !Anemometre
           OPEN(UNIT=57,FILE=anemometre_h_file, FORM='formatted', POSITION = 'append', &
           STATUS='unknown')
           WRITE(57,103) time, Hn(anemo_h(:),1,:), Hn(anemo_h(:),2,:), Hn(anemo_h(:),3,:), Hn(anemo_h(:),4,:), &
                               Hn(anemo_h(:),5,:), Hn(anemo_h(:),6,:)
           CLOSE(57)

        ENDIF
        !Traces au format PlotMTV
        IF (mod(it,freq_plot) == 0) THEN
           CALL trace_profile(H_mesh, Hn, it, freq_plot, list_mode, 'Hma') 
           CALL trace_profile(phi_mesh, phin, it, freq_plot, list_mode, 'phi') 
    !Dipole et quadrupole
           en_file = 'dip_quad_h'
           OPEN(UNIT=58,FILE=en_file, FORM='formatted', POSITION = 'append', &
           STATUS='unknown')
           CALL moments(H_mesh, list_mode, Hn, dipole, quadripole)
           WRITE(58,103) time, dipole(:), &
                         quadripole(1,:), quadripole(2,:), quadripole(3,:)
           CLOSE(58)

        ENDIF
        IF (mod(it,freq_restart) == 0) THEN
           CALL sauvegarde(it,freq_restart)
! cas PRECESSION mxw pr tracer le champ de vitesse sur maillage en H
!          CALL write_restart_maxwell(H_mesh, phi_mesh, time, list_mode, vel, vel, phin, phin, file_vel, it, freq_restart)
!          CALL trace_profile(H_mesh, vel, it, freq_restart, list_mode, 'VEL')
! cas PRECESSION mxw pr tracer le champ de vitesse sur maillage en H--- fin
        ENDIF
    ENDIF
! mars 6

!**************************************************************************************
     !===================Put your postprocessing stuff here

  ENDDO ! FIN RESOLUTION
  tps = user_time(dummy) - tps
  WRITE(*,*) ' Temps total ', tps
  !---------------------------------------------------------------------------


  !------------------------POST PROCESSING CONVERGENCE------------------------
  IF (test_de_convergence) CALL post_proc_test(time)
  !---------------------------------------------------------------------------

  !------------------------POST PROCESSING NAVIER-STOKES----------------------
  IF (type_pb == 'nst' .OR. type_pb=='mhd') THEN
     !IF(rang==0) CALL plot_scalar_field(vv_mesh%jj, vv_mesh%rr, un(:,1,1),'ur1.plt')
     !IF(rang==0) CALL plot_scalar_field(vv_mesh%jj, vv_mesh%rr, un(:,2,1),'ur2.plt')
     !IF(rang==0) CALL plot_scalar_field(vv_mesh%jj, vv_mesh%rr, un(:,3,1),'ut1.plt')
     !IF(rang==0) CALL plot_scalar_field(vv_mesh%jj, vv_mesh%rr, un(:,4,1),'ut2.plt')
     !IF(rang==0) CALL plot_scalar_field(vv_mesh%jj, vv_mesh%rr, un(:,5,1),'uz1.plt')
     !IF(rang==0) CALL plot_scalar_field(vv_mesh%jj, vv_mesh%rr, un(:,6,1),'uz2.plt')
     !ALLOCATE(vect(2,vv_mesh%np))
     !vect(1,:) = un(:,1,1); vect(2,:) = un(:,5,1)
     !IF(rang==0) CALL plot_arrow_label(vv_mesh%jj, vv_mesh%rr, vect,'vect.plt')
    !Dipole et quadrupole
    !==========
      en_file = 'vort_moyenne'
      OPEN(UNIT=54,FILE=en_file, FORM='formatted', STATUS='unknown')
      WRITE(54,*)'# rotV ds espace phys'
      WRITE(54,*)'# Temps vort_moyenne(1:3)'
        CALL vorticite(vv_mesh, list_mode, un, vort_moyenne)
        WRITE(54,103) time, vort_moyenne(:)
      CLOSE(54)
 END IF
  !---------------------------------------------------------------------------


  !------------------------POST MAXWELL---------------------------------------
  IF (type_pb == 'mxw' .OR. type_pb == 'mhd' .OR. type_pb == 'mxx') THEN
   DO i = 1, m_max_c
     CALL norme_interface(H_mesh,phi_mesh,interface_H_phi,mu_H_field,mu_phi,Hn(:,:,i),phin(:,:,i),list_mode(i),err)
     WRITE(*,*) ' erreur on the interface H_phi', err, 'for the mode ',list_mode(i)
     CALL norme_interface_H_mu(H_mesh,interface_H_mu,mu_H_field,Hn(:,:,i),list_mode(i),err)
     WRITE(*,*) ' erreur on the interface H_mu', err, 'for the mode ',list_mode(i)
   ENDDO
     !IF (rang==0)  CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn(:,1,1),'Hn1.plt')
     !IF (rang==0)  CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn(:,2,1),'Hn2.plt')
     !IF (rang==0)  CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn(:,3,1),'Hn3.plt')
     !IF (rang==0)  CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn(:,4,1),'Hn4.plt')
     !IF (rang==0)  CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn(:,5,1),'Hn5.plt')
     !IF (rang==0)  CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn(:,6,1),'Hn6.plt')
     !IF (rang==0)  CALL plot_scalar_field(phi_mesh%jj, phi_mesh%rr, phin(:,1,1),'phin1.plt')
     !IF (rang==0)  CALL plot_scalar_field(phi_mesh%jj, phi_mesh%rr, phin(:,2,1),'phin2.plt')
  END IF 
  !---------------------------------------------------------------------------


  !----------------SAUVEGARDE-------------------------------------------------
  CALL sauvegarde(it,freq_restart) 
  !---------------------------------------------------------------------------


  !------------------------FIN------------------------------------------------
  CALL MPI_FINALIZE(code)
  !---------------------------------------------------------------------------

END PROGRAM mhd_prog
