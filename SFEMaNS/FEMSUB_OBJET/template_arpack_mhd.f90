MODULE arpack_mhd

  PUBLIC :: solver_arpack_mhd
  PRIVATE

CONTAINS

  SUBROUTINE solver_arpack_mhd(test_de_convergence,H_mesh,phi_mesh,&
       mode_max_c,dt,list_mode,mu_H_field,interface_H_mu, mu_phi, interface_H_phi)
    USE chaine_caractere
    USE def_type_mesh
    USE initialisation
    IMPLICIT NONE
    LOGICAL,                       INTENT(IN) :: test_de_convergence
    TYPE(mesh_type),               INTENT(IN) :: H_mesh, phi_mesh
    INTEGER,                       INTENT(IN) :: mode_max_c
    REAL(KIND=8),                  INTENT(IN) :: dt, mu_phi
    INTEGER, POINTER, DIMENSION(:),INTENT(IN) :: list_mode
    REAL(KIND=8),DIMENSION(:),     INTENT(IN) :: mu_H_field
    TYPE(interface_type),          INTENT(IN) :: interface_H_mu, interface_H_phi
    !Arpack---------------------------------------------------------------------
    LOGICAL                                   :: if_arpack
    CHARACTER(len=3)                          :: arpack_type
    INTEGER                                   :: code
    !---------------------------------------------------------------------------

    !------------------------ARPACK OR NOT? THAT IS THE QUESTION----------------
    if_arpack = .FALSE.
    IF (.NOT.test_de_convergence) THEN
       OPEN(UNIT = 21, FILE = 'data', FORM = 'formatted', STATUS = 'unknown')
       CALL read_until(21, 'data_arpack')
       READ(21,*) if_arpack
       READ(21,*) arpack_type
       CLOSE(21)
    END IF
    !---------------------------------------------------------------------------


    !------------------------ARPACK RESOLUTION----------------------------------
    IF (if_arpack) THEN
       IF (arpack_type=='nst') THEN
          CALL arpack_navier_stokes
       ELSE IF (arpack_type=='max') THEN
          CALL arpack_maxwell(H_mesh,phi_mesh,mode_max_c,dt,list_mode,&
              mu_H_field,interface_H_mu, mu_phi, interface_H_phi)
          !CALL arpack_maxwell_int_by_parts(H_mesh,phi_mesh,mode_max_c,dt,list_mode,mu_H_field,interface_H_mu)
       END IF
       !----------------SAUVEGARDE-------------------------------------------------
       CALL sauvegarde(0,1)
       !---------------------------------------------------------------------------
       CALL MPI_FINALIZE(code)
       STOP
    END IF
    !---------------------------------------------------------------------------


    RETURN
  END SUBROUTINE solver_arpack_mhd

  !---------------------------------------------------------------------------
  SUBROUTINE arpack_maxwell_int_by_parts(H_mesh, phi_mesh, mode_max_c, dt, list_mode, &
       mu_H_field, interface_H_mu)
    USE initialisation
    USE def_type_mesh
    USE chaine_caractere
    USE sub_plot
    USE fem_tn_NS_MHD
    USE post_processing
    IMPLICIT NONE
    TYPE(mesh_type)                           :: H_mesh, phi_mesh
    INTEGER,                       INTENT(IN) :: mode_max_c
    REAL(KIND=8),                  INTENT(IN) :: dt
    INTEGER, POINTER, DIMENSION(:),INTENT(IN) :: list_mode
    REAL(KIND=8),DIMENSION(:),      INTENT(IN):: mu_H_field
    TYPE(interface_type),      INTENT(IN)     :: interface_H_mu

    CHARACTER(len=2)                          :: WHICH
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: eigen
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: eigen_vect
    REAL(KIND=8)                              :: tol_arpack
    REAL(KIND=8)                              :: rho, theta, rmax, norm, err
    LOGICAL                                   :: redo
    INTEGER :: nb_vp, iter_max, i, n, ndim, nconv, nmax, k
    CHARACTER(LEN=3) :: ifile, jfile
    CHARACTER(LEN=7) :: kfile    
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:):: vect, ev_Hn, rot_Hn
    REAL(KIND=8)                              :: time=0.
    REAL(KIND=8), DIMENSION(4)                :: fpar_sp
    INTEGER,      DIMENSION(4)                :: ipar_sp

    OPEN(UNIT = 21, FILE = 'data', FORM = 'formatted', STATUS = 'unknown')
    CALL read_until(21, 'solver_arpack')
    READ(21,*) nb_vp
    READ(21,*) iter_max
    READ(21,*) tol_arpack
    READ(21,*) which
    READ(21,*) redo
    CLOSE(21)

    WRITE(*,*) 'WHICH=', WHICH
    ndim = 12*H_mesh%np
    ALLOCATE(eigen(nb_vp,2),eigen_vect(ndim,nb_vp), ev_Hn(H_mesh%np,6,1), rot_Hn(H_mesh%np,6,1))
    ALLOCATE(vect(6,H_mesh%np,1))

    OPEN (UNIT=21, FILE = 'data', FORM='formatted', STATUS='unknown')
    CALL read_until(21, 'data_solver_maxwell')
    DO k=1,4
      READ(21,*) ipar_sp(k)
    END DO
    DO k=1,2
      READ(21,*) fpar_sp(k)
    END DO

    ! ATTENTION
    DO i = 1, mode_max_c ! ATTENTION
       ! ATTENTION
       CALL arpack_not_sym(nb_vp, nconv, iter_max, tol_arpack, &
            prodmat_maxwell_int_by_parts, eigen, eigen_vect, which, i, list_mode, redo)
       !Postprocessing

       rmax = -1.d10
       DO n = 1, MIN(nconv,nb_vp)
          rho = SQRT(eigen(n,1)**2+eigen(n,2)**2)
          theta = ATAN2(eigen(n,2),eigen(n,1))
          WRITE(*,*) '****************************************************'
          WRITE(*,'(A,e12.5,A,e12.5,A,i3,A,i3)') 'Partie Reelle ', LOG(rho)/dt, &
               ' et Partie Imaginaire ', theta/dt, ' de la vp ', n, ' pour le mode ',list_mode(i)
          WRITE(*,*) '****************************************************'
          WRITE(ifile,'(I3)') n
          WRITE(jfile,'(I3)') list_mode(i)
          kfile=ifile(start_of_string(ifile):)//'_'//jfile(start_of_string(jfile):)
          
          DO k = 1, 6
             ev_Hn(:,k,1) = eigen_vect((k-1)*H_mesh%np+1:k*H_mesh%np,n)
          END DO

          CALL  write_restart_maxwell_mode_int_by_parts(H_mesh, phi_mesh, time, list_mode(i), &
                ev_Hn, ev_Hn, kfile, 0, 1)

          CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, &
               eigen_vect(1:H_mesh%np,n),&
               'Hrc'//kfile(start_of_string(kfile):last_of_string(kfile))//'.plt')
          CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, &
               eigen_vect(H_mesh%np+1:2*H_mesh%np,n),&
               'Hrs'//kfile(start_of_string(kfile):last_of_string(kfile))//'.plt')
          CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, &
               eigen_vect(2*H_mesh%np+1:3*H_mesh%np,n),&
               'Htc'//kfile(start_of_string(kfile):last_of_string(kfile))//'.plt')
          CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, &
               eigen_vect(3*H_mesh%np+1:4*H_mesh%np,n),&
               'Hts'//kfile(start_of_string(kfile):last_of_string(kfile))//'.plt')
          CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, &
               eigen_vect(4*H_mesh%np+1:5*H_mesh%np,n),&
               'Hzc'//kfile(start_of_string(kfile):last_of_string(kfile))//'.plt')
          CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, &
               eigen_vect(5*H_mesh%np+1:6*H_mesh%np,n),&
               'Hzs'//kfile(start_of_string(kfile):last_of_string(kfile))//'.plt')
          !March 3, 2010, JLG + FL
          CALL calcul_rot_champ_vect(H_mesh, list_mode(i:i), ev_Hn, rot_Hn, ipar_sp, fpar_sp)
          CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, rot_Hn(1:H_mesh%np,1,1),&
               'RotH_rc'//kfile(start_of_string(kfile):last_of_string(kfile))//'.plt')
          CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, rot_Hn(1:H_mesh%np,2,1),&
               'RotH_rs'//kfile(start_of_string(kfile):last_of_string(kfile))//'.plt')
          CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, rot_Hn(1:H_mesh%np,3,1),&
               'RotH_tc'//kfile(start_of_string(kfile):last_of_string(kfile))//'.plt')
          CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, rot_Hn(1:H_mesh%np,4,1),&
               'RotH_ts'//kfile(start_of_string(kfile):last_of_string(kfile))//'.plt')
          CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, rot_Hn(1:H_mesh%np,5,1),&
               'RotH_zc'//kfile(start_of_string(kfile):last_of_string(kfile))//'.plt')
          CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, rot_Hn(1:H_mesh%np,6,1),&
               'RotH_zs'//kfile(start_of_string(kfile):last_of_string(kfile))//'.plt')
          !March 3, 2010, JLG + FL

          CALL norme_interface_H_mu(H_mesh,interface_H_mu,mu_H_field,ev_Hn(:,:,i),list_mode(i),err)
          WRITE(*,'(A,e12.5,A,i3,A,i3)') ' erreur on the interface H_mu', err, &
               ' de la vp ', n, ' pour le mode ',list_mode(i)

          DO k = 1, 6
             ev_Hn(:,k,1) = mu_H_field(:)*eigen_vect((k-1)*H_mesh%np+1:k*H_mesh%np,n)
          END DO
          CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, &
               ev_Hn(1:H_mesh%np,1,1),&
               'Brc'//kfile(start_of_string(kfile):last_of_string(kfile))//'.plt')
          CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, &
               ev_Hn(1:H_mesh%np,2,1),&
               'Brs'//kfile(start_of_string(kfile):last_of_string(kfile))//'.plt')
          CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, &
               ev_Hn(1:H_mesh%np,3,1),&
               'Btc'//kfile(start_of_string(kfile):last_of_string(kfile))//'.plt')
          CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, &
               ev_Hn(1:H_mesh%np,4,1),&
               'Bts'//kfile(start_of_string(kfile):last_of_string(kfile))//'.plt')
          CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, &
               ev_Hn(1:H_mesh%np,5,1),&
               'Bzc'//kfile(start_of_string(kfile):last_of_string(kfile))//'.plt')
          CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, &
               ev_Hn(1:H_mesh%np,6,1),&
               'Bzs'//kfile(start_of_string(kfile):last_of_string(kfile))//'.plt')
!
          IF (LOG(rho).GE.rmax) THEN
             nmax = n
             rmax = LOG(rho)
          END IF

          DO k = 1, 6
             vect(k,:,1) = mu_H_field(:)*eigen_vect((k-1)*H_mesh%np+1:k*H_mesh%np,n)
          END DO
          norm = norme_H1_champ(H_mesh,list_mode,vect)+norme_L2_champ(H_mesh,list_mode,vect)
          err  = norme_div(H_mesh, list_mode, vect)
          WRITE(*,*) '----------------------------------------------------'
          WRITE(*,'(A,e12.5,A,i3,A,i3)') '|div(mu Hn)|_L2/|mu Hn|_H1 = ', err/norm, &
               ' de la vp ', n, ' pour le mode ',list_mode(i)
          WRITE(*,*) '----------------------------------------------------'

       END DO

       IF (nconv>0) THEN
          rho = SQRT(eigen(nmax,1)**2+eigen(nmax,2)**2)
          theta = ATAN2(eigen(nmax,2),eigen(nmax,1))
          WRITE(*,*) ' '
          WRITE(*,'(A,i3,A,i3)') 'nmax(vp de plus gde partie reelle) ', nmax,' pr le mode = ', list_mode(i)
          WRITE(*,'(A,e12.5,A,e12.5)') ' Partie Reelle Max = ', LOG(rho)/dt, &
               ' Partie Imaginaire de la vp Max = ', theta/dt
       END IF
       !End postprocessing
    END DO ! Fin boucle sur les modes

    DEALLOCATE(eigen,eigen_vect,vect)

  END SUBROUTINE arpack_maxwell_int_by_parts
  !-----------------------------------------------------------------

  !-----------------------------------------------------------------
  SUBROUTINE arpack_maxwell(H_mesh, phi_mesh, mode_max_c, dt, list_mode,mu_H_field,interface_H_mu,mu_phi,interface_H_phi)
    USE initialisation
    USE def_type_mesh
    USE chaine_caractere
    USE sub_plot
    USE fem_tn_NS_MHD
    USE post_processing
    IMPLICIT NONE
    TYPE(mesh_type)                           :: H_mesh, phi_mesh
    INTEGER,                       INTENT(IN) :: mode_max_c
    REAL(KIND=8),                  INTENT(IN) :: dt, mu_phi
    INTEGER, POINTER, DIMENSION(:),INTENT(IN) :: list_mode
    REAL(KIND=8),DIMENSION(:),      INTENT(IN):: mu_H_field
    TYPE(interface_type),      INTENT(IN)     :: interface_H_mu,interface_H_phi

    CHARACTER(len=2)                          :: WHICH
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: eigen
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: eigen_vect
    REAL(KIND=8)                              :: tol_arpack
    REAL(KIND=8)                              :: rho, theta, rmax, norm, err, err1
    LOGICAL                                   :: redo
    INTEGER :: nb_vp, iter_max, i, n, ndim, nconv, nmax, k
    CHARACTER(LEN=3) :: ifile, jfile
    CHARACTER(LEN=7) :: kfile    
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: vect, ev_Hn, ev_phin, rot_Hn, grad_phi
    REAL(KIND=8)                               :: time=0.
    REAL(KIND=8), DIMENSION(4)                :: fpar_sp
    INTEGER,      DIMENSION(4)                :: ipar_sp
    OPEN(UNIT = 21, FILE = 'data', FORM = 'formatted', STATUS = 'unknown')
    CALL read_until(21, 'solver_arpack')
    READ(21,*) nb_vp
    READ(21,*) iter_max
    READ(21,*) tol_arpack
    READ(21,*) which
    READ(21,*) redo
    CLOSE(21)

    WRITE(*,*) 'WHICH=', WHICH
    ndim = 12*H_mesh%np + 4*phi_mesh%np
    ALLOCATE(eigen(nb_vp,2),eigen_vect(ndim,nb_vp),ev_Hn(H_mesh%np,6,1),ev_phin(phi_mesh%np,2,1))
    ALLOCATE(rot_Hn(H_mesh%np,6,1))
    ALLOCATE(vect(6,H_mesh%np,1))
    ALLOCATE(grad_phi(phi_mesh%np,6,1))
    ! ATTENTION
    IF (mode_max_c/=1) THEN
       WRITE(*,*) ' BUG: mode_max_c/=1 '
       STOP
    END IF
    DO i = 1, mode_max_c ! ATTENTION
       ! ATTENTION
       CALL arpack_not_sym(nb_vp, nconv, iter_max, tol_arpack, &
            prodmat_maxwell, eigen, eigen_vect, which, i, list_mode, redo)
       !Postprocessing

       rmax = -1.d10
       DO n = 1, MIN(nconv,nb_vp)
          rho = SQRT(eigen(n,1)**2+eigen(n,2)**2)
          theta = ATAN2(eigen(n,2),eigen(n,1))
          WRITE(*,*) '****************************************************'
          WRITE(*,'(A,e12.5,A,e12.5,A,i3,A,i3)') 'Partie Reelle ', LOG(rho)/dt, &
               ' et Partie Imaginaire ', theta/dt, ' de la vp ', n, ' pour le mode ',list_mode(i)
          WRITE(*,*) '****************************************************'
          WRITE(ifile,'(I3)') n
          WRITE(jfile,'(I3)') list_mode(i)
          kfile=ifile(start_of_string(ifile):)//'_'//jfile(start_of_string(jfile):)
          
          DO k = 1, 6
             ev_Hn(:,k,1) = eigen_vect((k-1)*H_mesh%np+1:k*H_mesh%np,n)
          END DO
          ev_phin(:,1,1) = eigen_vect(6*H_mesh%np+1:6*H_mesh%np+phi_mesh%np,n)
          ev_phin(:,2,1) = eigen_vect(6*H_mesh%np+phi_mesh%np+1:6*H_mesh%np+2*phi_mesh%np,n)


          CALL  write_restart_maxwell_mode(H_mesh, phi_mesh, time, list_mode(i), &
                ev_Hn, ev_Hn, ev_phin, ev_phin, kfile, 0, 1)

          CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, &
               eigen_vect(1:H_mesh%np,n),&
               'Hrc'//kfile(start_of_string(kfile):last_of_string(kfile))//'.plt')
          CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, &
               eigen_vect(H_mesh%np+1:2*H_mesh%np,n),&
               'Hrs'//kfile(start_of_string(kfile):last_of_string(kfile))//'.plt')
          CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, &
               eigen_vect(2*H_mesh%np+1:3*H_mesh%np,n),&
               'Htc'//kfile(start_of_string(kfile):last_of_string(kfile))//'.plt')
          CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, &
               eigen_vect(3*H_mesh%np+1:4*H_mesh%np,n),&
               'Hts'//kfile(start_of_string(kfile):last_of_string(kfile))//'.plt')
          CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, &
               eigen_vect(4*H_mesh%np+1:5*H_mesh%np,n),&
               'Hzc'//kfile(start_of_string(kfile):last_of_string(kfile))//'.plt')
          CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, &
               eigen_vect(5*H_mesh%np+1:6*H_mesh%np,n),&
               'Hzs'//kfile(start_of_string(kfile):last_of_string(kfile))//'.plt')
          !March 3, 2010, JLG + FL
          CALL calcul_rot_champ_vect(H_mesh, list_mode(i:i), ev_Hn, rot_Hn, ipar_sp, fpar_sp)
          CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, rot_Hn(1:H_mesh%np,1,1),&
               'RotH_rc'//kfile(start_of_string(kfile):last_of_string(kfile))//'.plt')
          CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, rot_Hn(1:H_mesh%np,2,1),&
               'RotH_rs'//kfile(start_of_string(kfile):last_of_string(kfile))//'.plt')
          CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, rot_Hn(1:H_mesh%np,3,1),&
               'RotH_tc'//kfile(start_of_string(kfile):last_of_string(kfile))//'.plt')
          CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, rot_Hn(1:H_mesh%np,4,1),&
               'RotH_ts'//kfile(start_of_string(kfile):last_of_string(kfile))//'.plt')
          CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, rot_Hn(1:H_mesh%np,5,1),&
               'RotH_zc'//kfile(start_of_string(kfile):last_of_string(kfile))//'.plt')
          CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, rot_Hn(1:H_mesh%np,6,1),&
               'RotH_zs'//kfile(start_of_string(kfile):last_of_string(kfile))//'.plt')
          !March 3, 2010, JLG + FL

          !Recontruction of the magnetic fields
          CALL gradphi(phi_mesh, ev_phin(:,1:2,1:1), list_mode(i:i), 'data', grad_phi)
          CALL plot_two_scalar_field(H_mesh%jj, H_mesh%rr, mu_H_field*eigen_vect(1:H_mesh%np,n),&
          phi_mesh%jj, phi_mesh%rr, grad_phi(:,1,1),'BBrc'//kfile(start_of_string(kfile):last_of_string(kfile))//'.plt')
          CALL plot_two_scalar_field(H_mesh%jj, H_mesh%rr, mu_H_field*eigen_vect(H_mesh%np+1:2*H_mesh%np,n),&
          phi_mesh%jj, phi_mesh%rr, grad_phi(:,2,1),'BBrs'//kfile(start_of_string(kfile):last_of_string(kfile))//'.plt')
          CALL plot_two_scalar_field(H_mesh%jj, H_mesh%rr, mu_H_field*eigen_vect(2*H_mesh%np+1:3*H_mesh%np,n),&
          phi_mesh%jj, phi_mesh%rr, grad_phi(:,3,1),'BBtc'//kfile(start_of_string(kfile):last_of_string(kfile))//'.plt')
          CALL plot_two_scalar_field(H_mesh%jj, H_mesh%rr, mu_H_field*eigen_vect(3*H_mesh%np+1:4*H_mesh%np,n),&
          phi_mesh%jj, phi_mesh%rr, grad_phi(:,4,1),'BBts'//kfile(start_of_string(kfile):last_of_string(kfile))//'.plt')
          CALL plot_two_scalar_field(H_mesh%jj, H_mesh%rr, mu_H_field*eigen_vect(4*H_mesh%np+1:5*H_mesh%np,n),&
          phi_mesh%jj, phi_mesh%rr, grad_phi(:,5,1),'BBzc'//kfile(start_of_string(kfile):last_of_string(kfile))//'.plt')
          CALL plot_two_scalar_field(H_mesh%jj, H_mesh%rr, mu_H_field*eigen_vect(5*H_mesh%np+1:6*H_mesh%np,n),&
          phi_mesh%jj, phi_mesh%rr, grad_phi(:,6,1),'BBzs'//kfile(start_of_string(kfile):last_of_string(kfile))//'.plt')
          !End reconstruction magnetic fields

          CALL norme_interface_H_mu(H_mesh,interface_H_mu,mu_H_field,ev_Hn(:,:,i),list_mode(i),err)
          CALL norme_interface(H_mesh,phi_mesh,interface_H_phi,mu_H_field,&
                   mu_phi,ev_Hn(:,:,i),ev_phin(:,:,i),list_mode(i),err1)
          !WRITE(*,'(A,e12.5,A,i3,A,i3)') ' erreur on the interface H_mu', err, &
          !     ' de la vp ', n, ' pour le mode ',list_mode(i)
          WRITE(*,'(A,e12.5,A,e12.5,A,i3,A,i3)') 'Erreur sur interface H_mu', err, &
               ' et erreur sur interface H_phi', err1, ' de la vp', n, ' pour le mode',list_mode(i)
          DO k = 1, 6
             ev_Hn(:,k,1) = mu_H_field(:)*eigen_vect((k-1)*H_mesh%np+1:k*H_mesh%np,n)
          END DO
          CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, &
               ev_Hn(1:H_mesh%np,1,1),&
               'Brc'//kfile(start_of_string(kfile):last_of_string(kfile))//'.plt')
          CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, &
               ev_Hn(1:H_mesh%np,2,1),&
               'Brs'//kfile(start_of_string(kfile):last_of_string(kfile))//'.plt')
          CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, &
               ev_Hn(1:H_mesh%np,3,1),&
               'Btc'//kfile(start_of_string(kfile):last_of_string(kfile))//'.plt')
          CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, &
               ev_Hn(1:H_mesh%np,4,1),&
               'Bts'//kfile(start_of_string(kfile):last_of_string(kfile))//'.plt')
          CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, &
               ev_Hn(1:H_mesh%np,5,1),&
               'Bzc'//kfile(start_of_string(kfile):last_of_string(kfile))//'.plt')
          CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, &
               ev_Hn(1:H_mesh%np,6,1),&
               'Bzs'//kfile(start_of_string(kfile):last_of_string(kfile))//'.plt')
!
          IF (LOG(rho).GE.rmax) THEN
             nmax = n
             rmax = LOG(rho)
          END IF
! TEST
!         WRITE(*,*) '----------------------------------------------------'
!         WRITE(*,*) 'max mu_H_field= ', MAXVAL(mu_H_field(:))
!         WRITE(*,*) '----------------------------------------------------'
          DO k = 1, 6
             vect(k,:,1) = mu_H_field(:)*eigen_vect((k-1)*H_mesh%np+1:k*H_mesh%np,n)
          END DO
          norm = norme_H1_champ(H_mesh,list_mode,vect)+norme_L2_champ(H_mesh,list_mode,vect)
          err  = norme_div(H_mesh, list_mode, vect)
          WRITE(*,*) '----------------------------------------------------'
          WRITE(*,'(A,e12.5,A,i3,A,i3)') '|div(mu Hn)|_L2/|mu Hn|_H1 = ', err/norm, ' de la vp ', n, ' pour le mode ',list_mode(i)
          WRITE(*,*) '----------------------------------------------------'

       END DO

       IF (nconv>0) THEN
          rho = SQRT(eigen(nmax,1)**2+eigen(nmax,2)**2)
          theta = ATAN2(eigen(nmax,2),eigen(nmax,1))
          WRITE(*,*) ' '
          WRITE(*,'(A,i3,A,i3)') 'nmax(vp de plus gde partie reelle) ', nmax,' pr le mode = ', list_mode(i)
          WRITE(*,'(A,e12.5,A,e12.5)') ' Partie Reelle Max = ', LOG(rho)/dt, &
               ' Partie Imaginaire de la vp Max = ', theta/dt
          WRITE(jfile,'(I3)') list_mode(i)
          kfile=jfile(start_of_string(jfile):)
          CALL plot_scalar_field(phi_mesh%jj, phi_mesh%rr, &
               eigen_vect(6*H_mesh%np+1:6*H_mesh%np+phi_mesh%np,nmax),&
               'phic'//kfile(start_of_string(kfile):last_of_string(kfile))//'.plt')
          CALL plot_scalar_field(phi_mesh%jj, phi_mesh%rr, &
               eigen_vect(6*H_mesh%np+phi_mesh%np+1:6*H_mesh%np+2*phi_mesh%np,nmax),&
               'phis'//kfile(start_of_string(kfile):last_of_string(kfile))//'.plt')
       END IF
       !End postprocessing
    END DO ! Fin boucle sur les modes

    DEALLOCATE(eigen,eigen_vect,vect)

  END SUBROUTINE arpack_maxwell
  !-----------------------------------------------------------------

  !-----------------------------------------------------------------
  SUBROUTINE arpack_navier_stokes
    RETURN
  END SUBROUTINE arpack_navier_stokes
  !-----------------------------------------------------------------

  !-----------------------------------------------------------------
  SUBROUTINE arpack_not_sym(nb_vp, nconv, iter_max, tol_arpack, prodmat, &
       eigen, eigen_vect, which, imode, list_mode, redo)
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    INTERFACE
       SUBROUTINE prodmat(x,b,m,i) 
         INTEGER,                     INTENT(IN)    :: m, i
         REAL(KIND=8), DIMENSION(m),  INTENT(IN)    :: x
         REAL(KIND=8), DIMENSION(m),  INTENT(INOUT) :: b
       END SUBROUTINE prodmat
    END INTERFACE
    INTEGER,                          INTENT(IN)    :: nb_vp, iter_max, imode
    REAL(KIND=8),                     INTENT(INOUT) :: tol_arpack
    REAL(KIND=8), DIMENSION(nb_vp,2), INTENT(OUT)   :: eigen
    REAL(KIND=8), DIMENSION(:,:),     INTENT(OUT)   :: eigen_vect
    INTEGER, POINTER, DIMENSION(:),   INTENT(IN)    :: list_mode
    LOGICAL,                          INTENT(IN)    :: redo

    !ARPACK --------------------------------------------------------------------
    INTEGER                                         :: IDO, INFO                                 
    CHARACTER(len=1)                                :: BMAT
    CHARACTER(len=2)                                :: WHICH
    INTEGER, DIMENSION(11)                          :: IPARAM, IPNTR
    INTEGER                                         :: NDIM, NEV, NCV, LWORKL, NCONV
    REAL(KIND=8)                                    :: SIGMAR, sigmai
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:)         :: RESID, WORKD, WORKL
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:)         :: WORKEV
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:)       :: V, d
    LOGICAL,      ALLOCATABLE, DIMENSION(:)         :: SEL
    INTEGER                                         :: NP, NEV2, i, j, unit_save, &
         ndim_r, nb_vp_r, code, nb_procs, ido_loop
    REAL(KIND=8)                                    :: rnorm
    INTEGER,      ALLOCATABLE, DIMENSION(:)         :: tableau_ido
    LOGICAL                                         :: bidon 
    !-------------END OF DECLARATIONS-------------------------------------------

    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,code)
    ALLOCATE(tableau_ido(nb_procs))

    NEV   = nb_vp
    NDIM = SIZE(eigen_vect,1)
    NCV = 2*NEV+1
    !NCV = NEV + 2 ! Does not work with P2P2
    lworkl  = 3*ncv**2+6*ncv 
    ALLOCATE(V(NDIM,NCV), WORKL(lworkl), WORKD(3*NDIM), RESID(NDIM))
    ALLOCATE(SEL(NCV), d(ncv,3), workev(3*ncv))

    IF (redo) THEN
       unit_save=200+list_mode(imode)
       REWIND(unit_save)
       READ(unit_save) ido, tol_arpack, iparam, ipntr, info, &
            np, rnorm, nconv, nev2, ndim_r, bmat, nb_vp_r
       IF (ndim_r/=ndim .OR. nb_vp_r/=nb_vp) THEN
          WRITE(*,*) ' BUG dans restart de arpack_mhd '
          STOP
       END IF
       READ(unit_save) (resid(i), i = 1, ndim)
       READ(unit_save) (workd(i), i = 1, 3*ndim)
       READ(unit_save) (workl(i), i = 1, lworkl)
       DO j=1,ncv
          READ(unit_save) (v(i,j), i = 1, ndim)
       ENDDO

    ELSE
       IDO=0
       BMAT = 'I'
       INFO=0
       IPARAM(1)=1
       IPARAM(2)=1
       IPARAM(4)=1
       IPARAM(5)=1
       IPARAM(6)=1
       IPARAM(7)=1 
       IPARAM(8)=1 
       IPARAM(9)=1
       IPARAM(10)=1
       IPARAM(11)=1

       !----INITIALIZE SO THAT BC ARE ENFORCED CORRECTLY---------------------------
       WORKD(1:ndim) = 1.d0
       CALL prodmat(WORKD(1:ndim), RESID, ndim, imode)
       WORKD(ndim/2+1:ndim)=WORKD(1:ndim/2)
       INFO = 1
       !---------------------------------------------------------------------------
    END IF

    IPARAM(3)=iter_max
    bidon = .FALSE.
    ido_loop = -99
    !---------------------------------------------------------------------------
    DO WHILE(ido_loop.NE.99)

       IF (.NOT.bidon) THEN
          !CALL DNAUPD(IDO,BMAT,NDIM,WHICH,NEV, TOL_arpack, &
          !     RESID, NCV, V, NDIM, IPARAM, IPNTR,&
          !     WORKD, WORKL, LWORKL, INFO)
          CALL DNAUPD_CHCKPT(IDO,BMAT,NDIM,WHICH,NEV, TOL_arpack, &
               RESID, NCV, V, NDIM, IPARAM, IPNTR,&
               WORKD, WORKL, LWORKL, INFO, &
               np, rnorm, nconv, nev2)
       END IF

       CALL MPI_ALLGATHER(ido, 1,  MPI_INTEGER,  tableau_ido, 1, MPI_INTEGER, &
            MPI_COMM_WORLD, code)

       !GESTION DU ID0=-2
       IF (MINVAL(ABS(tableau_ido+2))==0) THEN
          ido = -2
       END IF

       !GESTION DU ID0=99
       IF (ido==99) THEN
          bidon=.TRUE.
          CALL prodmat(WORKD(IPNTR(1):), WORKD(IPNTR(2):), ndim, imode)
       END IF
       IF (MAXVAL(ABS(tableau_ido-99))/=0) THEN
          ido_loop = -99
       ELSE
          ido_loop = 99
       END IF

       IF(IDO==-1 .OR. IDO==1) THEN
          CALL prodmat(WORKD(IPNTR(1):), WORKD(IPNTR(2):), ndim, imode)
       ENDIF

       IF (ido == -2) THEN
          !
          !            %---------------------------------------------------%
          !            | After maxitr iterations without convergence,      |
          !            | output the computed quantities to the file state. |
          !            %---------------------------------------------------%
          !
          PRINT *,'***** DEBUT SAUVEGARDE # ', list_mode(imode), imode
          unit_save=200+list_mode(imode)
          REWIND(unit_save)
          WRITE(unit_save) ido, tol_arpack, iparam, ipntr, info, &
               np, rnorm, nconv, nev2, ndim, bmat, nb_vp
          WRITE(unit_save) (resid(i), i = 1, ndim)
          WRITE(unit_save) (workd(i), i = 1, 3*ndim)
          WRITE(unit_save) (workl(i), i = 1, lworkl)
          DO j=1,ncv
             WRITE(unit_save) (v(i,j), i = 1, ndim)
          ENDDO
          nconv = -1
          PRINT *,'***** FIN SAUVEGARDE # ', list_mode(imode), imode
          DEALLOCATE(tableau_ido)
          RETURN
       ENDIF

    END DO

    NCONV = IPARAM(5)

    IF (NCONV/=NEV) THEN
       WRITE(*,*) ' Only ', nconv, 'eigenvalues have been computed instead  of', nev
       WRITE(*,*) ' it_max not large enough '
    END IF

    IF(NCONV.NE.0) THEN
       CALL dneupd (.TRUE., 'All', SEL, d, d(1,2), v, NDIM, &
            sigmar, sigmai,  workev, bmat, nDIM, which, nev, tol_arpack, &
            resid, ncv, v, ndim, iparam, ipntr, workd, workl, lworkl, INFO)

       eigen_vect = v(:,1:nconv)
       eigen(:,1) = d(1:nconv,1)
       eigen(:,2) = d(1:nconv,2)
    END IF
    DEALLOCATE(tableau_ido)
    DEALLOCATE(V, WORKL, WORKD, RESID)
    DEALLOCATE(SEL, d, workev)

  END SUBROUTINE arpack_not_sym
  !-----------------------------------------------------------------

  SUBROUTINE write_restart_maxwell_mode_int_by_parts(H_mesh, phi_mesh, time, list_mode, Hn, Hn1, &
       filename, it, freq_restart)
    USE def_type_mesh
    USE chaine_caractere
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    TYPE(mesh_type), TARGET                                    :: H_mesh,phi_mesh     
    REAL(KIND=8),                                   INTENT(IN) :: time 
    INTEGER,                                        INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:),                 INTENT(IN) :: Hn, Hn1
    CHARACTER(LEN=7),                               INTENT(IN) :: filename 
    INTEGER,                                        INTENT(IN) :: it, freq_restart

    OPEN(UNIT = 10, FILE = 'suite_maxwell'//'_'//filename, POSITION='append', &
         FORM = 'unformatted', STATUS = 'replace') 
    WRITE(10) time, H_mesh%np , phi_mesh%np ,1 , 1

    WRITE(10) list_mode
    WRITE(10) Hn(:,:,1)
    WRITE(10) Hn1(:,:,1)
    CLOSE(10)        
  END SUBROUTINE write_restart_maxwell_mode_int_by_parts

  SUBROUTINE write_restart_maxwell_mode(H_mesh, phi_mesh, time, list_mode, Hn, Hn1, &
       phin, phin1, filename, it, freq_restart)

    USE def_type_mesh
    USE chaine_caractere
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    TYPE(mesh_type), TARGET                                    :: H_mesh,phi_mesh     
    REAL(KIND=8),                                   INTENT(IN) :: time 
    INTEGER,                                        INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:),                 INTENT(IN) :: Hn, Hn1
    REAL(KIND=8), DIMENSION(:,:,:),                 INTENT(IN) :: phin, phin1
    CHARACTER(LEN=7),                               INTENT(IN) :: filename 
    INTEGER,                                        INTENT(IN) :: it, freq_restart

             OPEN(UNIT = 10, FILE = 'suite_maxwell'//'_'//filename, POSITION='append', &
                  FORM = 'unformatted', STATUS = 'replace') 
             WRITE(10) time, H_mesh%np , phi_mesh%np ,1 , 1

             WRITE(10) list_mode
             WRITE(10) Hn(:,:,1)
             WRITE(10) Hn1(:,:,1)
             WRITE(10) phin(:,:,1)
             WRITE(10) phin1(:,:,1)
             CLOSE(10)        

  END SUBROUTINE write_restart_maxwell_mode

SUBROUTINE calcul_rot_champ_vect(mesh,list_mode, v, rotv,ipar_sp_in,fpar_sp_in)

!Calul du rotationnel du champ vectoriel v
!mesh : maillage
! v(mesh%np,6,size(list_mode))
! rotv(mesh%np,6,size(list_mode))

      USE def_type_mesh
      USE prep_maill
      USE st_matrix
      USE solve_sp
      USE matrix_type
      USE fem_s_axi
      USE fem_s_axi_M
      USE chaine_caractere

      IMPLICIT NONE
      TYPE(mesh_type), TARGET                                         :: mesh
      INTEGER,DIMENSION(:),         INTENT(IN)                        :: list_mode       
      REAL(KIND=8), DIMENSION(mesh%np,6,size(list_mode)) ,INTENT(IN)  :: v 
      REAL(KIND=8), DIMENSION(mesh%np,6,size(list_mode)) ,INTENT(OUT) :: rotv  !
      INTEGER, OPTIONAL     ,DIMENSION(4),               INTENT(IN) :: ipar_sp_in
      REAL(KIND=8), OPTIONAL,DIMENSION(4),               INTENT(IN)   :: fpar_sp_in

      TYPE(matrice_bloc),                            SAVE     :: mat_mass_sg  !matrice single
      TYPE(matrice_bloc),                            SAVE     :: mat_mass_tot !matrice totale
      LOGICAL,                                       SAVE     :: once = .true.
      INTEGER,      DIMENSION(:), POINTER,           SAVE     :: ia, ja, ia_b, ja_b
      REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:),   SAVE     :: u0
      REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:)  ,   SAVE     :: ff, sol
      INTEGER                                                 :: i, j, k, n_inf, n_sup
!------------DECLARATIONS FOR SOLVE_SPARSEKIT--------------------------------------
      INTEGER,      DIMENSION(4),                SAVE :: ipar_sp
      REAL(KIND=8), DIMENSION(4),                SAVE :: fpar_sp


     
      IF (once) THEN
         once = .false.
       
         CALL st_csr(mesh%jj, mat_mass_sg%ia, mat_mass_sg%ja)
         CALL st_csr_bloc(mat_mass_sg%ia, mat_mass_sg%ja, &
         mat_mass_tot%ia, mat_mass_tot%ja, 6)
        
         ALLOCATE(mat_mass_tot%aa(SIZE(mat_mass_tot%ja)))
         !calcul de la matrice de masse (6*mesh%np) * (6*mesh%np)
         mat_mass_tot%aa = 0.d0
         CALL qs_mass_6comp(mesh, 1.d0, mat_mass_tot%ia, mat_mass_tot%ja, mat_mass_tot%aa)
!-------------READ PARAMETERS FOR PRECONDITIONING VELOCITY MATRIX--------------
         IF (.NOT.PRESENT(ipar_sp_in)) THEN
         OPEN (UNIT=21, FILE ='data', FORM='formatted', STATUS='unknown')
         CALL read_until(21, 'data_solver')
         DO k=1,4;      READ(21,*) ipar_sp(k);    END DO
            DO k=1,2;      READ(21,*) fpar_sp(k);    END DO
         CLOSE(21)
         ELSE
         ipar_sp = ipar_sp_in
         fpar_sp = fpar_sp_in
         ENDIF
         !no reordering for the moment, pb avec precond avoir !
         ipar_sp(4) = 0
         ALLOCATE(u0(mesh%np,6,size(list_mode)))
         ALLOCATE(ff(6*mesh%np,size(list_mode)))
         ALLOCATE(sol(6*mesh%np,size(list_mode)))
         u0 = 0.d0 ; ff = 0.d0 ; sol = 0.d0
      ENDIF

      DO j =1, size(list_mode)
         CALL qs_01_rot(mesh,list_mode(j),v(:,:,j),u0(:,:,j))
      ENDDO

     

      DO j = 1, size(list_mode)
         !changement de format
         DO i=1, 6
            n_inf = 1+(i-1)*mesh%np
            n_sup = i*mesh%np
            ff(n_inf:n_sup,j) = u0(:,i,j)
         ENDDO
         !inversion de la matrice de masse
         CALL solve_splib(ipar_sp, fpar_sp, mat_mass_tot, ff(:,j), sol(:,j))
         !changement de format
         DO i=1, 6
            n_inf = 1+(i-1)*mesh%np
            n_sup = i*mesh%np		
            rotv(:,i,j) = sol(n_inf:n_sup,j)
         ENDDO
      ENDDO

END SUBROUTINE calcul_rot_champ_vect
END MODULE arpack_mhd
