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

  !---------------------------------------------------------------------------

  IMPLICIT NONE
  include 'mpif.h'

  !Champs pour Navier-Stokes-------------------------------------------------
  TYPE(mesh_type), POINTER                        :: pp_mesh, vv_mesh    
  REAL(KIND=8), POINTER, DIMENSION(:,:,:)         :: un, pn !Vitesse, pression
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
  REAL(KIND=8)                                    :: err, norm
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:)       :: vect
  !-------------END OF DECLARATIONS-------------------------------------------

  !-------------POST_PROCESSING-----------------------------------------------
  INTEGER                                         :: freq_restart, freq_en, freq_plot
  INTEGER                                         :: select_ns, select_mxw
  !Declare whatever you want here for your own purpose
! mars 6
  !Taylor-Couette------------------------------------------------------------
  LOGICAL                                         :: if_select 
  LOGICAL                                         :: START_NL=.false.
  !Traitement des fichiers de sortie-----------------------------------------
  CHARACTER(LEN=2)                                :: truc
  INTEGER                                         :: dd_end, ff_end
  CHARACTER(LEN=100)                              :: en_file
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:)         :: e_c_u, e_cm_h
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:)         :: e_c_u_sym, e_c_u_anti, e_cm_h_sym, e_cm_h_anti
  REAL(KIND=8),              DIMENSION(3)         :: x_anemo_v, x_anemo_h
  REAL(KIND=8),              DIMENSION(5)         :: y_anemo_v, y_anemo_h
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:)       :: norme_u, norme_h
  INTEGER , ALLOCATABLE, DIMENSION(:)             :: anemo_v, anemo_h
  CHARACTER(LEN=100)                              :: e_c_u_file, e_cm_h_file
  CHARACTER(LEN=100)                              :: norme_u_file, norme_h_file
  CHARACTER(LEN=100)                              :: anemometre_v_file, anemometre_h_file
  REAL(KIND=8),DIMENSION(3)                       :: type_sym_u, type_sym_h
  INTEGER , ALLOCATABLE, DIMENSION(:)             :: point_sym_u, point_sym_h
! mars 6

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
  IF (type_pb == 'mxx') THEN
     !Vous gerez vel ici
     DO i = 1, m_max_c
        !List_mode contient les modes physiques
        IF(list_mode(i)==1) THEN
           vel(:,:,i) = 0.d0
        ENDIF
     END DO
     !IF (rang==0) CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, vel(:,3,1), 'vel3.plt')
  END IF 
  !---------------------------------------------------------------------------
!----------------RESTART AVEC MULTIPLICATION ENERGIE MODES H IMPAIRS POUR Taylor-Couette----
! avril 11
        IF (START_NL) THEN

         DO i = 1, m_max_c
          IF (mod((list_mode(i)+1),2)==0) THEN
            Hn(:,:,i)    = 30.d0*Hn(:,:,i)
            phin(:,:,i)  = 30.d0*phin(:,:,i)
          ENDIF
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
     !CALL trace_profile(vv_mesh, un, 0, freq_plot, list_mode, 'Uns')
     IF (rang==0) CALL plot_scalar_field(pp_mesh%jj,pp_mesh%rr, pn(:,1,1),'pn_init.plt')
     IF (rang==0) CALL plot_scalar_field(vv_mesh%jj,vv_mesh%rr, un(:,2,1),'un_init.plt')
     norm = norme_H1_champ_par(vv_mesh,list_mode,un)+ 1.d-15
     err = norme_div_par(vv_mesh,list_mode, un)
     IF (rang==0) WRITE(*,*) 'div(un) relative champ init   = ', err/norm
     IF (rang==0) write(*,*) 'div(un) absolue  champ init   = ', err
  END IF

  IF (type_pb == 'mxw' .OR. type_pb == 'mhd' .OR. type_pb == 'mxx') THEN

   DO i = 1, m_max_c
     CALL norme_interface(H_mesh,phi_mesh,interface_H_phi,mu_H_field,mu_phi,Hn(:,:,i),phin(:,:,i),list_mode(i),err)
     WRITE(*,*) ' erreur on the interface H_phi for initial data', err, 'for the mode ',list_mode(i)
     CALL norme_interface_H_mu(H_mesh,interface_H_mu,mu_H_field,Hn(:,:,i),list_mode(i),err)
     WRITE(*,*) ' erreur on the interface H_mu initial data', err, 'for the mode ',list_mode(i)
   ENDDO
! CARO a remettre
     !CALL trace_profile(H_mesh, Hn, 0, freq_plot, list_mode, 'Hin')
     !CALL trace_profile(H_mesh, vel, 0, freq_plot, list_mode, 'Vci')
     IF (rang==0)  CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn(:,1,1),'Hn1_init.plt')
     !IF (rang==0)  CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn(:,2,1),'Hn2_init.plt')
     !IF (rang==0)  CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn(:,3,1),'Hn3_init.plt')
     !IF (rang==0)  CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn(:,4,1),'Hn4_init.plt')
     !IF (rang==0)  CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn(:,5,1),'Hn5_init.plt')
     !IF (rang==0)  CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn(:,6,1),'Hn6_init.plt')
     IF (rang==0)  CALL plot_scalar_field(phi_mesh%jj, phi_mesh%rr, phin(:,1,1),'phin1_init.plt')
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
       !=======
     type_sym_u(1) = 1.d0
     type_sym_u(2) = 1.d0
     type_sym_u(3) =-1.d0
     ALLOCATE(point_sym_u(vv_mesh%np))
     CALL calcul_point_sym(vv_mesh, point_sym_u)

       en_file = 'e_c_u'
       dd_end = last_c_leng (20, en_file)
       DO i = 1, m_max_c
          dd_end = last_c_leng (20, en_file)
          WRITE(truc,'(i2)') list_mode(i)
          ff_end = start_of_string (truc)
          en_file = en_file(1:dd_end)//'_'//truc(ff_end:)
       END DO
       e_c_u_file = en_file
       OPEN(UNIT=20,FILE=e_c_u_file, FORM='formatted', STATUS='unknown')
       WRITE(20,*)'#energie par mode'
       WRITE(20,*)'# Temps    valeur par ordre croissant de mode -> '
       ALLOCATE(e_c_u(m_max_c),e_c_u_sym(m_max_c),e_c_u_anti(m_max_c))
       CALL val_ener(vv_mesh, list_mode, un, e_c_u)
       CALL val_ener_sym_glob(vv_mesh, list_mode, un, point_sym_u, e_c_u, e_c_u_sym, e_c_u_anti, type_sym_u)
           WRITE(20,103) time, e_c_u, e_c_u_sym, e_c_u_anti
       CLOSE(20)
       !composante
       !==========
       en_file = 'norme_comp_u'
       dd_end = last_c_leng (20, en_file)
       DO i = 1, m_max_c
          dd_end = last_c_leng (20, en_file)
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
     type_sym_h(1) =-1.d0
     type_sym_h(2) =-1.d0
     type_sym_h(3) = 1.d0
     ALLOCATE(point_sym_h(H_mesh%np))
     CALL calcul_point_sym(H_mesh, point_sym_h)

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
       CALL val_ener_sym_glob(H_mesh, list_mode, Hn, point_sym_h, e_cm_h, e_cm_h_sym, e_cm_h_anti, type_sym_h)
           WRITE(50,103) time, e_cm_h, e_cm_h_sym, e_cm_h_anti
       CLOSE(50)
       !composante
       !==========
       en_file = 'norme_comp_h'
       dd_end = last_c_leng (20, en_file)
       DO i = 1, m_max_c
          dd_end = last_c_leng (20, en_file)
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

    ENDIF
   
    !Anemometre
    !==========
 
    !Navier-Stokes
    ALLOCATE(anemo_v(15))
    IF (type_pb=='nst' .OR. type_pb=='mhd') THEN   
       !anemo_v(1)  = find_point(vv_mesh,x_anemo_v(1),y_anemo_v(1)) 
       !anemo_v(2)  = find_point(vv_mesh,x_anemo_v(2),y_anemo_v(1))  
       !anemo_v(3)  = find_point(vv_mesh,x_anemo_v(3),y_anemo_v(1))
       !anemo_v(4)  = find_point(vv_mesh,x_anemo_v(1),y_anemo_v(2))
       !anemo_v(5)  = find_point(vv_mesh,x_anemo_v(2),y_anemo_v(2))
       !anemo_v(6)  = find_point(vv_mesh,x_anemo_v(3),y_anemo_v(2))
       !anemo_v(7)  = find_point(vv_mesh,x_anemo_v(1),y_anemo_v(3))    
       !anemo_v(8)  = find_point(vv_mesh,x_anemo_v(2),y_anemo_v(3))    
       !anemo_v(9)  = find_point(vv_mesh,x_anemo_v(3),y_anemo_v(3))
       !anemo_v(10) = find_point(vv_mesh,x_anemo_v(1),y_anemo_v(4))
       !anemo_v(11) = find_point(vv_mesh,x_anemo_v(2),y_anemo_v(4))
       !anemo_v(12) = find_point(vv_mesh,x_anemo_v(3),y_anemo_v(4))
       !anemo_v(13) = find_point(vv_mesh,x_anemo_v(1),y_anemo_v(5))  
       !anemo_v(14) = find_point(vv_mesh,x_anemo_v(2),y_anemo_v(5))  
       !anemo_v(15) = find_point(vv_mesh,x_anemo_v(3),y_anemo_v(5))

       !en_file = 'anemometre_V'
       !dd_end = last_c_leng (20, en_file)
       !DO i = 1, m_max_c
       !   dd_end = last_c_leng (20, en_file)
       !   WRITE(truc,'(i2)') list_mode(i)
       !   ff_end = start_of_string (truc)
       !   en_file = en_file(1:dd_end)//'_'//truc(ff_end:)
       !END DO
       !anemometre_v_file = en_file
       !OPEN(UNIT=56,FILE=anemometre_v_file, FORM='formatted', STATUS='unknown')
       !WRITE(56,*)'#Anemometre pour le champ de vitesse'
       !WRITE(56,*)'#Coordonnees des points :    '
       !DO i = 1, 15
       !    WRITE(56,*) '# point ',i, ' : r=', vv_mesh%rr(1,anemo_v(i)), '; z = ', vv_mesh%rr(2,anemo_v(i)) 
       !ENDDO 
       !CLOSE(56)
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
       dd_end = last_c_leng (20, en_file)
       DO i = 1, m_max_c
          dd_end = last_c_leng (20, en_file)
          WRITE(truc,'(i2)') list_mode(i)
          ff_end = start_of_string (truc)
          en_file = en_file(1:dd_end)//'_'//truc(ff_end:)
       END DO
       anemometre_h_file = en_file
       OPEN(UNIT=57,FILE=anemometre_h_file, FORM='formatted', STATUS='unknown')
       WRITE(57,*) '#Anemometre pour le champ magnetique'
       WRITE(57,*) '#Coordonnees des points :    '
       !DO i = 1, 15
       !    WRITE(57,*) '# point ',i, 'r=', H_mesh%rr(1,anemo_h(i)), '; z = ', H_mesh%rr(2,anemo_h(i))
       !ENDDO
       CLOSE(57)
    ENDIF

    103 FORMAT(150(e22.9,2x))
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

!April 19th JLG
!TESTTTTTTT
  iF (it==2) tps = user_time(dummy)
!TESTTTTTTT
     time = time + dt    
     WRITE(*,*) 'Iteration = ',it,' ;  Time = ', time

     !======Navier Stokes  
     IF (type_pb == 'nst') THEN

        CALL navier_stokes(time, un, pn)
write(*,*) ' Navier_stokes  done '
        !TEST pour Navier_stokes
!TESTTTTTTT
!April 19th JLG
!        norm = norme_H1_champ_par(vv_mesh,list_mode,un)+1.d-15
!        err = norme_div_par(vv_mesh,list_mode, un)
!        IF(rang==0) WRITE(*,*) 'div(un) relative    = ', err/norm
!        IF(rang==0) WRITE(*,*) 'div(un) absolue     = ', err
!        DO i = 1, m_max_c
!           unit = rang*m_max_c+10+i
!           WRITE(unit,*) time, norme_l2_champ(vv_mesh, list_mode(i:i), un(:,:,i:i))
!        END DO
!April 19th JLG
!TESTTTTTTT

     ELSE IF (type_pb == 'mxw' .OR. type_pb == 'mxx') THEN

       CALL maxwell(time, Hn, phin)

        !TEST pour Maxwell
        norm = norme_H1_champ_par(H_mesh,list_mode,Hn)+1.d-15
        err = norme_div_par(H_mesh,list_mode, Hn)
        IF(rang==0) WRITE(*,*) 'div(Hn) relative    = ', err/norm
        IF(rang==0) WRITE(*,*) 'div(Hn) absolue     = ', err
        DO i = 1, m_max_c
           unit = 100 + rang*m_max_c + i
           WRITE(unit,*) time, norme_l2_champ(H_mesh, list_mode(i:i), Hn(:,:,i:i))
        END DO


     ELSE IF (type_pb == 'mhd') THEN
        IF (if_select) THEN
            CALL mhd(time, un, pn, Hn, phin, select_mode_ns, select_mode_mxw)
        ELSE
            CALL mhd(time, un, pn, Hn, phin)
        ENDIF
        !TEST pour Navier_stokes
        norm = norme_H1_champ_par(vv_mesh,list_mode,un)+1.d-15
        err = norme_div_par(vv_mesh,list_mode, un)
        IF(rang==0) WRITE(*,*) 'div(un) relative    = ', err/norm
        IF(rang==0) WRITE(*,*) 'div(un) absolue     = ', err
        DO i = 1, m_max_c
           unit = rang*m_max_c+10+i
           WRITE(unit,*) time, &
                norme_l2_champ(vv_mesh, list_mode(i:i), un(:,:,i:i))
        END DO

        !TEST pour Maxwell
        norm = norme_H1_champ_par(H_mesh,list_mode,Hn)+1.d-15
        err = norme_div_par(H_mesh,list_mode, Hn)
        IF(rang==0) WRITE(*,*) 'div(Hn) relative    = ', err/norm
        IF(rang==0) WRITE(*,*) 'div(Hn) absolue     = ', err
        DO i = 1, m_max_c
           unit = 100 + rang*m_max_c + i
           WRITE(unit,*) time, &
                norme_l2_champ(H_mesh, list_mode(i:i), Hn(:,:,i:i))
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
!TESTTTTTTTTTTTT
!April 19th
write(*,*) 'ON passe ici'
!TESTTTTTTTTTTTT
           !Energie
           OPEN(UNIT=20,FILE=e_c_u_file, FORM='formatted', POSITION = 'append', &
           STATUS='unknown')
!          CALL val_ener(vv_mesh, list_mode, un, e_c_u)
           CALL val_ener_sym_glob(vv_mesh, list_mode, un, point_sym_u, e_c_u, e_c_u_sym, e_c_u_anti, type_sym_u)
           WRITE(20,103) time, e_c_u, e_c_u_sym, e_c_u_anti
           CLOSE(20)
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
!TESTTTTTTTTTTTT
!April 19th
write(*,*) 'ON passe ici'
!TESTTTTTTTTTTTT
           CALL trace_profile(vv_mesh, un, it, freq_plot, list_mode, 'Uns')   
        ENDIF
        !Ecriture fichier restart de securite
        IF (mod(it,freq_restart) == 0) THEN
!TESTTTTTTTTTTTT
!April 19th
write(*,*) 'ON passe ici'
!TESTTTTTTTTTTTT
           CALL sauvegarde(it,freq_restart)
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
           CALL val_ener_sym_glob(H_mesh, list_mode, Hn, point_sym_h, e_cm_h, e_cm_h_sym, e_cm_h_anti, type_sym_h)
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
        ENDIF
        IF (mod(it,freq_restart) == 0) THEN
           CALL sauvegarde(it,freq_restart)
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
     IF (rang==0)  CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn(:,1,m_max_c),'Hn1.plt')
     IF (rang==0)  CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn(:,2,m_max_c),'Hn2.plt')
     IF (rang==0)  CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn(:,3,m_max_c),'Hn3.plt')
     IF (rang==0)  CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn(:,4,m_max_c),'Hn4.plt')
     IF (rang==0)  CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn(:,5,m_max_c),'Hn5.plt')
     IF (rang==0)  CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn(:,6,m_max_c),'Hn6.plt')
     IF (rang==0)  CALL plot_scalar_field(phi_mesh%jj, phi_mesh%rr, phin(:,1,m_max_c),'phin1.plt')
     IF (rang==0)  CALL plot_scalar_field(phi_mesh%jj, phi_mesh%rr, phin(:,2,m_max_c),'phin2.plt')
  END IF 
  !---------------------------------------------------------------------------


  !----------------SAUVEGARDE-------------------------------------------------
  CALL sauvegarde(it,freq_restart) 
  !---------------------------------------------------------------------------


  !------------------------FIN------------------------------------------------
  CALL MPI_FINALIZE(code)
  !---------------------------------------------------------------------------

END PROGRAM mhd_prog
