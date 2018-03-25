PROGRAM mhd_prog

  USE def_type_mesh
  USE tn_parallele
  USE fem_tn_NS_MHD
  USE sub_plot
  USE initialisation
  USE boundary
  USE chaine_caractere
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
  INTEGER                                         :: it, k, i, mode, unit
  REAL(KIND=8)                                    :: user_time, tps, dummy
  REAL(KIND=8)                                    :: err, norm
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:)       :: vect
  !-------------END OF DECLARATIONS-------------------------------------------

  !-------------POST_PROCESSING-----------------------------------------------
  INTEGER                                         :: freq_restart, freq_en, freq_plot
  INTEGER                                         :: select_ns, select_mxw
  LOGICAL                                         :: if_select
  !Declare whatever you want here for your own purpose
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
write(*,*) ' ok'

  !-------------YOUR OWN INITIALIZATIONS--------------------------------------
  OPEN(UNIT = 21, FILE = 'data', FORM = 'formatted', STATUS = 'unknown')
  CALL read_until(21, 'data_postproc')
write(*,*) ' ok'
  READ(21,*) freq_restart, freq_en, freq_plot
  CLOSE(21)
  !OPEN(UNIT = 21, FILE = 'data', FORM = 'formatted', STATUS = 'unknown')
  !CALL read_until(21, 'my_stuff')
  !READ(21,*) stuff_1, stuff_2, whatever 
  !CLOSE(21)
  IF (type_pb == 'mxx') THEN !Vous gerez vel ici
     DO i = 1, m_max_c
        IF(list_mode(i)==1) THEN !List_mode contient les modes physiques
           vel(:,:,i) = 0.d0
        ENDIF
     END DO
     IF (rang==0) CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, vel(:,3,1), 'vel3.plt')
  END IF
  ! Here you can do other initializations
  ! like CALL prep_VKN(vv_mesh)
  ! etc ...
  !---------------------------------------------------------------------------


  !----------------------QUELQUES VERIFICATIONS------------------------------
  ! nst = Navier-Stokes
  ! mxw = Maxwell
  ! mxx = Maxwell + u_restart =.true. (defined in initialisation)
  ! mhd = MHD  
  IF (type_pb == 'nst' .OR. type_pb == 'mhd') THEN 
     IF (rang==0) CALL plot_scalar_field(pp_mesh%jj,pp_mesh%rr, pn(:,1,1),'pn_init.plt')
     IF (rang==0) CALL plot_scalar_field(vv_mesh%jj,vv_mesh%rr, un(:,2,1),'un_init.plt')
     norm = norme_H1_champ_par(vv_mesh,list_mode,un)+ 1.d-15
     err = norme_div_par(vv_mesh,list_mode, un)
     IF (rang==0) WRITE(*,*) 'div(un) relative champ init   = ', err/norm
     IF (rang==0) write(*,*) 'div(un) absolue  champ init   = ', err
  END IF

  IF (type_pb == 'mxw' .OR. type_pb == 'mhd' .OR. type_pb == 'mxx') THEN
     CALL norme_interface(H_mesh,phi_mesh,interface_H_phi,mu_H_field,mu_phi,Hn(:,:,1),phin(:,:,1),1,err)
     WRITE(*,*) ' erreur on the interface H_phi for initial data', err
     CALL norme_interface_H_mu(H_mesh,interface_H_mu,mu_H_field,Hn(:,:,1),list_mode(1),err)
     WRITE(*,*) ' erreur on the interface H_mu initial data', err
     IF (rang==0)  CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn(:,1,1),'Hn1_init.plt')
     IF (rang==0)  CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn(:,2,1),'Hn2_init.plt')
     IF (rang==0)  CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn(:,3,1),'Hn3_init.plt')
     IF (rang==0)  CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn(:,4,1),'Hn4_init.plt')
     IF (rang==0)  CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn(:,5,1),'Hn5_init.plt')
     IF (rang==0)  CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn(:,6,1),'Hn6_init.plt')
     IF (rang==0)  CALL plot_scalar_field(phi_mesh%jj, phi_mesh%rr, phin(:,1,1),'phin1_init.plt')
     IF (rang==0)  CALL plot_scalar_field(phi_mesh%jj, phi_mesh%rr, phin(:,2,1),'phin2_init.plt')
     norm = norme_H1_champ_par(H_mesh,list_mode,Hn)+ 1.d-15
     err = norme_div_par(H_mesh,list_mode, Hn)
     WRITE(*,*) 'div(Hn) relative champ init   = ', err/norm
     write(*,*) 'div(Hn) absolue  champ init   = ', err
  END IF
  !---------------------------------------------------------------------------


  !------------------------OPTIONAL HANDLING OF ZERO MODES--------------------
  IF (.NOT.test_de_convergence) THEN
     OPEN(UNIT = 21, FILE = 'data', FORM = 'formatted', STATUS = 'unknown')
     CALL read_until(21, 'data_select_mode_nuls')
     READ(21,*) if_select
     ! Here decide what you want to do
     ! Allocate select_mode_ns, select_mode_mxw
     ! and read these arrays
     CLOSE(21)
  END IF
  !---------------------------------------------------------------------------


  !------------------------RESOLUTION ----------------------------------------
  tps = user_time(dummy)
  DO it = 1, nb_iteration

     time = time + dt    
     WRITE(*,*) 'Iteration = ',it,' ;  Time = ', time

     !======Navier Stokes  
     IF (type_pb == 'nst') THEN
        CALL navier_stokes(time, un, pn)
        !TEST pour Navier_stoke
        norm = norme_H1_champ_par(vv_mesh,list_mode,un)+1.d-15
        err = norme_div_par(vv_mesh,list_mode, un)
        IF(rang==0) WRITE(*,*) 'div(un) relative    = ', err/norm
        IF(rang==0) WRITE(*,*) 'div(un) absolue     = ', err
        DO i = 1, m_max_c
           unit = rang*m_max_c+10+i
           WRITE(unit,*) time, norme_l2_champ(vv_mesh, list_mode(i:i), un(:,:,i:i))
        END DO

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

        CALL mhd(time, un, pn, Hn, phin)
        !CALL mhd(time, un, pn, Hn, phin, select_mode_ns, select_mode_mxw)

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

     !===================Put your postprocessing stuff here
     ! Whatever 
     IF (mod(it,freq_en) == 0) THEN
        !CALL your_stuff_here
     END IF
     !Traces au format PlotMTV
     IF (mod(it,freq_plot) == 0) THEN
        !CALL trace_whatever_you_want_here
     ENDIF
     !Ecriture fichier restart de securite
     IF (mod(it, freq_restart) == 0) THEN
        CALL  sauvegarde(it,freq_restart)
     ENDIF
     !===================Put your postprocessing stuff here

  ENDDO
  tps = user_time(dummy) - tps
  WRITE(*,*) ' Temps total ', tps
  !---------------------------------------------------------------------------


  !------------------------POST PROCESSING CONVERGENCE------------------------
  IF (test_de_convergence) CALL post_proc_test(time)
  !---------------------------------------------------------------------------


  !------------------------POST PROCESSING NAVIER-STOKES----------------------
  IF (type_pb == 'nst' .OR. type_pb=='mhd') THEN
     IF(rang==0) CALL plot_scalar_field(vv_mesh%jj, vv_mesh%rr, un(:,1,1),'ur1.plt')
     IF(rang==0) CALL plot_scalar_field(vv_mesh%jj, vv_mesh%rr, un(:,2,1),'ur2.plt')
     IF(rang==0) CALL plot_scalar_field(vv_mesh%jj, vv_mesh%rr, un(:,3,1),'ut1.plt')
     IF(rang==0) CALL plot_scalar_field(vv_mesh%jj, vv_mesh%rr, un(:,4,1),'ut2.plt')
     IF(rang==0) CALL plot_scalar_field(vv_mesh%jj, vv_mesh%rr, un(:,5,1),'uz1.plt')
     IF(rang==0) CALL plot_scalar_field(vv_mesh%jj, vv_mesh%rr, un(:,6,1),'uz2.plt')
     ALLOCATE(vect(2,vv_mesh%np))
     vect(1,:) = un(:,1,1); vect(2,:) = un(:,5,1)
     IF(rang==0) CALL plot_arrow_label(vv_mesh%jj, vv_mesh%rr, vect,'vect.plt')
  END IF
  !---------------------------------------------------------------------------


  !------------------------POST MAXWELL---------------------------------------
  IF (type_pb == 'mxw' .OR. type_pb == 'mhd' .OR. type_pb == 'mxx') THEN
     CALL norme_interface(H_mesh,phi_mesh,interface_H_phi,mu_H_field,mu_phi,Hn(:,:,1),phin(:,:,1),1,err)
     WRITE(*,*) ' erreur on the interface H_phi', err
     CALL norme_interface_H_mu(H_mesh,interface_H_mu,mu_H_field,Hn(:,:,1),list_mode(1),err)
     WRITE(*,*) ' erreur on the interface H_mu', err
     IF (rang==0)  CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn(:,1,1),'Hn1.plt')
     IF (rang==0)  CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn(:,2,1),'Hn2.plt')
     IF (rang==0)  CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn(:,3,1),'Hn3.plt')
     IF (rang==0)  CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn(:,4,1),'Hn4.plt')
     IF (rang==0)  CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn(:,5,1),'Hn5.plt')
     IF (rang==0)  CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn(:,6,1),'Hn6.plt')
     IF (rang==0)  CALL plot_scalar_field(phi_mesh%jj, phi_mesh%rr, phin(:,1,1),'phin1.plt')
     IF (rang==0)  CALL plot_scalar_field(phi_mesh%jj, phi_mesh%rr, phin(:,2,1),'phin2.plt')
  END IF 
  !---------------------------------------------------------------------------


  !----------------SAUVEGARDE-------------------------------------------------
  CALL sauvegarde(it,freq_restart) 
  !---------------------------------------------------------------------------


  !------------------------FIN------------------------------------------------
  CALL MPI_FINALIZE(code)
  !---------------------------------------------------------------------------

END PROGRAM mhd_prog
