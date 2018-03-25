!
!Authors Jean-Luc Guermond, Raphael Laguerre, Copyrights 2005
!Revised June 2008, Jean-Luc Guermond
!
MODULE update_maxwell

CONTAINS

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------

  SUBROUTINE maxwell_decouple(H_mesh, phi_mesh, INTERFACE, interface_H_mu, Hn, phin, Hn1, phin1, vel, &
       stab, sigma_in, R_fourier, index_fourier, mu_H_field, mu_phi, time, dt, Rem, list_mode, &
       nb_syst_lin, H_phi_per, data_fichier)

    USE def_type_mesh
    USE dyn_line
    USE chaine_caractere
    USE st_matrix
    USE solve_sp
    USE Dir_nodes
    !Jan 29 2007
    USE periodic
    !Jan 29 2007
    USE boundary
    USE fem_tn_NS_MHD
    USE fem_s_axi_M
    USE fem_s_axi
    USE tn_parallele
    USE sft_parallele

    IMPLICIT NONE

    TYPE(mesh_type),                INTENT(IN)     :: H_mesh, phi_mesh
    TYPE(interface_type),           INTENT(IN)     :: INTERFACE, interface_H_mu
    INTEGER,      DIMENSION(:),     INTENT(IN)     :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(INOUT)  :: vel  
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(INOUT)  :: Hn, Hn1
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(INOUT)  :: phin, phin1
    REAL(KIND=8), DIMENSION(3),     INTENT(IN)     :: stab 
    REAL(KIND=8),                   INTENT(IN)     :: R_fourier
    INTEGER,                        INTENT(IN)     :: index_fourier
    REAL(KIND=8),                   INTENT(IN)     :: mu_phi, time, dt, Rem     
    REAL(KIND=8), DIMENSION(:),     INTENT(IN)     :: sigma_in, mu_H_field
    INTEGER,      DIMENSION(:),     INTENT(IN)     :: nb_syst_lin
    CHARACTER(len=64),              INTENT(IN)     :: data_fichier
    !jan 29
    TYPE(periodic_type),            INTENT(IN)     :: H_phi_per
    !jan 29

    REAL(KIND=8), DIMENSION(:), ALLOCATABLE,       SAVE  :: sigma
    TYPE(matrice_bloc), ALLOCATABLE, DIMENSION(:), SAVE  :: H_phi_mat1, H_phi_mat2
    !Feb 2 2007
    TYPE(dyn_real_line),ALLOCATABLE, DIMENSION(:), SAVE  :: precond1, precond2
    !Feb 2 2007
    TYPE(dyn_int_line),                            SAVE  :: js_D_mat1, js_D_mat2
    TYPE(dyn_int_line),   DIMENSION(3),            SAVE  :: js_D_H
    TYPE(dyn_int_line),                            SAVE  :: js_D_phi 
    LOGICAL,                                       SAVE  :: once=.TRUE.
    INTEGER, DIMENSION(:), ALLOCATABLE,            SAVE  :: isave
    INTEGER,                                       SAVE  :: m_max_c, nb_tot_syst_lin
    LOGICAL,                                       SAVE  :: per = .FALSE.
    !------------DECLARATIONS FOR SOLVE_SPLIB--------------------------------------
    REAL(KIND=8), DIMENSION(4), SAVE:: fpar_sp
    INTEGER,      DIMENSION(4), SAVE:: ipar_sp
    !-------------END OF DECLARATIONS----------------------------------------------

    !jan 29 2007      
    TYPE(matrice_bloc),        DIMENSION(2)              :: tampon
    !jan 29 2007

    LOGICAL, ALLOCATABLE, DIMENSION(:,:)                 :: Dir_H      ! Type of BC
    LOGICAL, ALLOCATABLE, DIMENSION(:)                   :: Dir_phi 
    INTEGER, ALLOCATABLE, DIMENSION(:)                   :: list_dirichlet_sides,bord_periodic
    REAL(KIND=8), DIMENSION(H_mesh%np,6)                 :: f_H, rhs_H
    REAL(KIND=8), DIMENSION(phi_mesh%np,2)               :: f_phi, rhs_phi
    REAL(KIND=8), DIMENSION(3*H_mesh%np+phi_mesh%np)     :: fG1, solG1, fG2, solG2
    INTEGER,      DIMENSION(4)                           :: n_Dir
    REAL(KIND=8), DIMENSION(H_mesh%np,6,SIZE(list_mode)) :: NL, H_ext
    REAL(KIND=8), DIMENSION(3)                           :: temps_par
    INTEGER, ALLOCATABLE, DIMENSION(:)                   :: ndir
    !June 7 2008
    INTEGER, POINTER,     DIMENSION(:)                   :: cav_index
    !June 7 2008
    INTEGER        :: mode, j_size, k, p, i, l, n_cour, n, nb_dirichlet_sides, n_bord, side1, side2, n_cav
    REAL(KIND=8)   :: user_time, tps, dummy, nr_vel, tps_tot, tps_cumul
    !April 17th, 2008, JLG
    REAL(KIND=8) :: one_and_half
    DATA one_and_half/1.5d0/
    !April 17th, 2008, JLG
    LOGICAL,                                       SAVE  :: i_alpha
    INTEGER,                                       SAVE  :: n_alpha
    REAL(KIND=8),                                  SAVE  :: ampl_alpha
    INTEGER                                              :: unit
    !CARO mars 2008


    IF (once) THEN

       once = .FALSE.

       !-------------RESCALING DE SIGMA-----------------------------------------------
       ALLOCATE(sigma(SIZE(sigma_in)))
       sigma = sigma_in * Rem
       !------------------------------------------------------------------------------

       !-------------DIMENSIONS-------------------------------------------------------
       m_max_c = SIZE(list_mode)
       !------------------------------------------------------------------------------

       !---------------BOUNDARY CONDITIONS---------------------------------------------        
       OPEN (UNIT=21, FILE = data_fichier, FORM='formatted', STATUS='unknown')


       !---------------FORCING TERM---------------------------------------------
! CARO Forcage alpha

       CALL read_until(21, 'data_alpha')
       READ  (21, *)  i_alpha, ampl_alpha, n_alpha
       WRITE(*,*) 'i_alpha, ampl_alpha, n_alpha= ',i_alpha, ampl_alpha, n_alpha

       CALL read_until(21, 'data_periodic')
       READ  (21, *)  n_bord
       IF (n_bord>0) THEN
          ALLOCATE(bord_periodic(2*n_bord)) 
          DO n= 1, n_bord
             READ  (21, *)  bord_periodic(2*n-1), bord_periodic(2*n) 
          END DO
       END IF

       CALL read_until(21, 'data_condlim_maxwell')
       ALLOCATE (Dir_H(3,MAXVAL(H_mesh%sides)))
       ALLOCATE (Dir_phi(MAXVAL(phi_mesh%sides)))
       Dir_H   = .FALSE.
       Dir_phi = .FALSE.
       DO k = 1, 3
          READ(21,*) nb_dirichlet_sides
          IF (nb_dirichlet_sides .LE. 0) THEN
             READ(21,*)
          ELSE
             ALLOCATE(list_dirichlet_sides(nb_dirichlet_sides))
             READ(21,*) list_dirichlet_sides
             DO i = 1, nb_dirichlet_sides ! Check coherence of data
                IF (MINVAL(ABS(H_mesh%sides - list_dirichlet_sides(i))) /= 0) THEN
                   WRITE(*,*) ' BUG list_dirichlet_sides not coherent with H mesh'
                   STOP
                END IF
             END DO
             Dir_H(k,list_dirichlet_sides) = .TRUE.
             ! Verify coeherence with Periodic BCs
             DO n = 1, n_bord
                IF (MINVAL(ABS(list_dirichlet_sides - bord_periodic(2*n-1))) == 0 .OR.  &
                     MINVAL(ABS(list_dirichlet_sides - bord_periodic(2*n))) == 0) THEN
                   WRITE(*,*) ' Donnees Dirichlet et Periodiques incoherentes pour H'
                   STOP
                END IF
             END DO
             ! End Verify coeherence with Periodic BCs
             DEALLOCATE(list_dirichlet_sides)
          END IF
       ENDDO

       READ(21,*) nb_dirichlet_sides
       IF (nb_dirichlet_sides .LE. 0) THEN
          READ(21,*)
       ELSE
          ALLOCATE(list_dirichlet_sides(nb_dirichlet_sides))
          READ(21,*) list_dirichlet_sides
          DO i = 1, nb_dirichlet_sides ! Check coherence of data
             IF (MINVAL(ABS(phi_mesh%sides - list_dirichlet_sides(i))) /= 0) THEN
                WRITE(*,*) ' BUG list_dirichlet_sides not coherent with phi mesh'
                STOP
             END IF
          END DO
          Dir_phi(list_dirichlet_sides) = .TRUE.
          ! Verify coeherence with Periodic BCs
          DO n = 1, n_bord
             IF (MINVAL(ABS(list_dirichlet_sides - bord_periodic(2*n-1))) == 0 .OR.  &
                  MINVAL(ABS(list_dirichlet_sides - bord_periodic(2*n))) == 0) THEN
                WRITE(*,*) ' Donnees Dirichlet et Periodiques incoherentes pour phi'
                STOP
             END IF
          END DO
          ! End Verify coeherence with Periodic BCs
          DEALLOCATE(list_dirichlet_sides)
       END IF
       IF (ALLOCATED(bord_periodic)) DEALLOCATE(bord_periodic)

       !remplissage des tableaux des noeuds contraints, cond Dirichlet
       CALL dirichlet_nodes(H_mesh%jjs, H_mesh%sides, Dir_H(1,:), js_D_H(1)%DIL)
       CALL dirichlet_nodes(H_mesh%jjs, H_mesh%sides, Dir_H(2,:), js_D_H(2)%DIL)
       CALL dirichlet_nodes(H_mesh%jjs, H_mesh%sides, Dir_H(3,:), js_D_H(3)%DIL)
       CALL dirichlet_nodes(phi_mesh%jjs, phi_mesh%sides, Dir_phi(:), js_D_phi%DIL)

       DEALLOCATE(Dir_H,Dir_phi)
       !Jan 29 2007
       IF (H_phi_per%n_bord/=0) THEN
          per = .TRUE.
       ELSE
          per = .FALSE.
       END IF
       !Jan 29 2007

       !------------------------------------------------------------------------------

       !-------------BOUNDARY CONDITIONS----------------------------------------------
       ! nombre de points contraints pour chaque composante
       n_Dir(1) = SIZE(js_D_H(1)%DIL)
       n_Dir(2) = SIZE(js_D_H(2)%DIL)
       n_Dir(3) = SIZE(js_D_H(3)%DIL)
       n_Dir(4) = SIZE(js_D_phi%DIL)
       IF (SIZE(js_D_phi%DIL)==0 .AND. R_fourier.LE.0.d0) THEN
          DEALLOCATE(js_D_phi%DIL)
          ALLOCATE(js_D_phi%DIL(1))
          js_D_phi%DIL(1)=phi_mesh%np
          n_DIR(4) = 1
       END IF

       !Detect cavities FOR DURAND'S SOLUTION, June 7, 2008
       CALL detect_cavities(INTERFACE, phi_mesh, cav_index)
       n_cav = MAXVAL(cav_index)
       IF (n_cav == 1) THEN
          WRITE(*,*) ' I HAVE DETECTED ', 1, ' cavity '
       ELSE IF (n_cav > 1) THEN
          WRITE(*,*) ' I HAVE DETECTED ', MAXVAL(cav_index), ' cavities'
       END IF
       IF (n_cav .GE. 1) THEN
          ALLOCATE(ndir(SIZE(js_D_phi%DIL)))
          ndir = js_D_phi%DIL
          DEALLOCATE(js_D_phi%DIL)
          n_Dir(4) = n_Dir(4) + n_cav
          ALLOCATE(js_D_phi%DIL(n_Dir(4)))
          js_D_phi%DIL(1:n_Dir(4)-n_cav)=ndir
          DEALLOCATE(ndir)
          DO i = 1, n_cav !Loop on the number of cavities
             DO k = 1, phi_mesh%me
                IF (phi_mesh%i_d(k)==cav_index(i)) THEN
                   n = phi_mesh%jj(1,k)
                   EXIT
                END IF
             END DO
             js_D_phi%DIL(n_Dir(4)-n_cav+i) = n
             WRITE(*,'(A,2(e15.5,3x))') ' Phi problem regularized at ', phi_mesh%rr(:,n)
          END DO
       END IF
       !Detect cavities FOR DURAND'S SOLUTION, June 7, 2008

       !remplissage du tableau global des noeuds pour la matrice
       ALLOCATE(js_D_mat1%DIL(n_Dir(1)+n_Dir(2)+n_Dir(3)+n_Dir(4)))
       ALLOCATE(js_D_mat2%DIL(n_Dir(1)+n_Dir(2)+n_Dir(3)+n_Dir(4)))
       !comp Hr
       js_D_mat1%DIL(1:n_Dir(1))                   = js_D_H(1)%DIL(:)
       js_D_mat2%DIL(1:n_Dir(1))                   = js_D_H(1)%DIL(:)
       !comp Hth
       n_cour = n_Dir(1)
       js_D_mat1%DIL(n_cour+1:n_cour+n_Dir(2)) = js_D_H(2)%DIL(:) + H_mesh%np
       js_D_mat2%DIL(n_cour+1:n_cour+n_Dir(2)) = js_D_H(2)%DIL(:) + H_mesh%np
       !comp Hz
       n_cour = n_Dir(1)+n_Dir(2)
       js_D_mat1%DIL(n_cour+1:n_cour+n_Dir(3)) = js_D_H(3)%DIL(:) + 2*H_mesh%np
       js_D_mat2%DIL(n_cour+1:n_cour+n_Dir(3)) = js_D_H(3)%DIL(:) + 2*H_mesh%np
       !comp phi
       n_cour = n_Dir(1)+n_Dir(2)+n_Dir(3)
       js_D_mat1%DIL(n_cour+1:n_cour+n_Dir(4)) = js_D_phi%DIL(:) + 3*H_mesh%np 
       js_D_mat2%DIL(n_cour+1:n_cour+n_Dir(4)) = js_D_phi%DIL(:) + 3*H_mesh%np                 
       !------------------------------------------------------------------------------

       !-------------MATRIX STRUCTURING FOR MAGNETIC UNKNOWNS-------------------------
       !CALL st_scr_maxwell_tmag3d_decouple(H_mesh, phi_mesh, interface, tampon(1)%ia, tampon(1)%ja)
       CALL st_scr_maxwell_tmag3d_decouple_mu(H_mesh, phi_mesh, INTERFACE, interface_H_mu, tampon(1)%ia, tampon(1)%ja)
       !------------------------------------------------------------------------------

       !-------------MATRIX ALLOCATION------------------------------------------------
       ALLOCATE(tampon(2)%ia(SIZE(tampon(1)%ia)), tampon(2)%ja(SIZE(tampon(1)%ja)))
       tampon(2)%ia = tampon(1)%ia
       tampon(2)%ja = tampon(1)%ja
       ALLOCATE(tampon(1)%aa(SIZE(tampon(1)%ja)))
       ALLOCATE(tampon(2)%aa(SIZE(tampon(2)%ja)))
       !H_phi_mat1 pour les comps H1, H4, H5 et phi1
       ALLOCATE(H_phi_mat1(m_max_c))
       !H_phi_mat2 pour les comps H2, H3, H6 et phi2
       ALLOCATE(H_phi_mat2(m_max_c))
       !------------------------------------------------------------------------------

       !-------------MATRIX COMPUTATION-----------------------------------------------
       !Feb 2 2007
       ALLOCATE(precond1(m_max_c), precond2(m_max_c))
       !Feb 2 2007
       ALLOCATE(isave(2*m_max_c))
       DO i = 1, m_max_c !Boucle sur les modes
          IF (i==1) THEN
             IF (per) THEN
                CALL st_csr_per(H_phi_per, tampon(1)%ia, tampon(1)%ja, H_phi_mat1(1)%ia, H_phi_mat1(1)%ja)
                ALLOCATE(H_phi_mat2(1)%ia(SIZE(H_phi_mat1(1)%ia)), H_phi_mat2(1)%ja(SIZE(H_phi_mat1(1)%ja)))
                H_phi_mat2(1)%ia = H_phi_mat1(1)%ia
                H_phi_mat2(1)%ja = H_phi_mat1(1)%ja
             ELSE
                ALLOCATE(H_phi_mat1(1)%ia(SIZE(tampon(1)%ia)), H_phi_mat1(1)%ja(SIZE(tampon(1)%ja)))
                H_phi_mat1(1)%ia = tampon(1)%ia
                H_phi_mat1(1)%ja = tampon(1)%ja
                ALLOCATE(H_phi_mat2(1)%ia(SIZE(tampon(1)%ia)), H_phi_mat2(1)%ja(SIZE(tampon(1)%ja)))
                H_phi_mat2(1)%ia = tampon(1)%ia
                H_phi_mat2(1)%ja = tampon(1)%ja
             END IF
          ELSE
             ALLOCATE(H_phi_mat1(i)%ia(SIZE(H_phi_mat1(i-1)%ia)), H_phi_mat1(i)%ja(SIZE(H_phi_mat1(i-1)%ja)))
             H_phi_mat1(i)%ia = H_phi_mat1(i-1)%ia
             H_phi_mat1(i)%ja = H_phi_mat1(i-1)%ja
             ALLOCATE(H_phi_mat2(i)%ia(SIZE(H_phi_mat2(i-1)%ia)), H_phi_mat2(i)%ja(SIZE(H_phi_mat2(i-1)%ja)))
             H_phi_mat2(i)%ia = H_phi_mat2(i-1)%ia
             H_phi_mat2(i)%ja = H_phi_mat2(i-1)%ja
          END IF

          mode = list_mode(i)
          isave(2*i-1) = -(2*i-1)
          isave(2*i) = -2*i
          ALLOCATE(H_phi_mat1(i)%aa(SIZE(H_phi_mat1(i)%ja)))
          ALLOCATE(H_phi_mat2(i)%aa(SIZE(H_phi_mat2(i)%ja)))

          tampon(1)%aa  =0.d0
          tampon(2)%aa  =0.d0

          CALL mat_maxwell(H_mesh,phi_mesh,INTERFACE, interface_H_mu, &
               mode, tampon(1)%ia, tampon(1)%ja, tampon(1)%aa, tampon(2)%aa, &
               mu_H_field, mu_phi, one_and_half/dt, stab, sigma, R_fourier, index_fourier)
          !Take care of discontinuous mu
          IF (interface_H_mu%mes >0) THEN
             CALL mat_maxwell_mu(H_mesh, interface_H_mu, mode, stab, &
                  tampon(1)%ia, tampon(1)%ja, tampon(1)%aa, tampon(2)%aa, mu_H_field, sigma)
          END IF
          IF (per) THEN
             CALL bc_per_M(H_phi_per, tampon(1), H_phi_mat1(i))
             CALL bc_per_M(H_phi_per, tampon(2), H_phi_mat2(i))
          ELSE
             H_phi_mat1(i)%aa = tampon(1)%aa
             H_phi_mat2(i)%aa = tampon(2)%aa
          END IF

          !Feb 2 2007
          ALLOCATE(precond1(i)%DRL(SIZE(H_phi_mat1(i)%ia)-1))
          ALLOCATE(precond2(i)%DRL(SIZE(H_phi_mat2(i)%ia)-1))
          CALL precond_mat(H_phi_mat1(i)%ia, H_phi_mat1(i)%aa, precond1(i)%DRL)
          CALL precond_mat(H_phi_mat2(i)%ia, H_phi_mat2(i)%aa, precond2(i)%DRL)
          !Feb 2 2007

          CALL Dirichlet_M (js_D_mat1%DIL, H_phi_mat1(i)%ia, H_phi_mat1(i)%ja, H_phi_mat1(i)%aa)
          CALL Dirichlet_M (js_D_mat2%DIL, H_phi_mat2(i)%ia, H_phi_mat2(i)%ja, H_phi_mat2(i)%aa)
       ENDDO
       isave = isave - nb_syst_lin(2)   ! nb_syst_lin(2) = syst lin deja alloues ailleurs
       nb_tot_syst_lin = nb_syst_lin(1) ! nb_syst_lin(1) = nb tot syst lin deja alloues ailleurs
       DEALLOCATE(tampon(1)%ia,tampon(1)%ja,tampon(1)%aa,tampon(2)%ia,tampon(2)%ja,tampon(2)%aa)
       !------------------------------------------------------------------------------

       !-------------READ PARAMETERS FOR PRECONDITIONING MATRIX-----------------------
       CALL read_until(21, 'data_solver_maxwell')
       DO k=1,4;      READ(21,*) ipar_sp(k)
       END DO
       DO k=1,2;      READ(21,*) fpar_sp(k)
       END DO
       !------------------------------------------------------------------------------

       CLOSE(21)
    ENDIF ! fin du .once.

    tps_tot = user_time(dummy)
    tps_cumul = 0
    !-------------TRANSPORT TERM---------------------------------------------------
    tps = user_time(dummy)
    nr_vel = norme_L2_champ_par(H_mesh,list_mode,vel)
    H_ext = 2*Hn-Hn1
    IF (nr_vel .LE. 1.d-10) THEN
       NL = 0.d0
    ELSE
       CALL FFT_PAR_CROSS_PROD(vel, H_ext, NL, temps_par)
    ENDIF
! CARO
          DO i = 1, m_max_c
            IF(i_alpha) THEN
              IF((n_alpha==1)) THEN ! alpha VKS
!             IF((n_alpha==1).AND.(list_mode(i)==0)) THEN ! alpha VKS
                DO k = 3,4
                  DO n = 1, H_mesh%np
!                   NL(n,k,i) = NL(n,k,i) + ampl_alpha * alpha(H_mesh%rr(:,n)) * H_ext(n,k,i)
                    NL(n,k,i) = NL(n,k,i) + &
                                (ampl_alpha/(1.d0+SUM(H_ext(100,:,i)**2))) * &
                                alpha(H_mesh%rr(:,n)) * H_ext(n,k,i)
                  ENDDO
                ENDDO
              unit = 30+list_mode(i)
              WRITE(unit,*) time, ampl_alpha/(1.d0+SUM(H_ext(100,:,i:i)**2))
              ELSE IF((n_alpha==2).AND.(list_mode(i)==1)) THEN ! alpha Karlsruhe pr m=1!
                DO k = 1,4
                  DO n = 1, H_mesh%np
                    NL(n,k,i) = NL(n,k,i) + ampl_alpha * alpha_K(H_mesh%rr(:,n)) * H_ext(n,k,i)
                  ENDDO
                ENDDO
              ELSE IF((n_alpha==3)) THEN ! alpha EsR
                DO k = 1,4
                  DO n = 1, H_mesh%np
                    NL(n,k,i) = NL(n,k,i) + ampl_alpha * alpha_EsR(H_mesh%rr(:,n)) * H_ext(n,k,i)
                  ENDDO
                ENDDO

              ENDIF
            ENDIF
          ENDDO
    tps = user_time(dummy) - tps; tps_cumul=tps_cumul+tps
    !WRITE(*,*) ' Tps NLS_SFT Maxwell', tps
    !------------------------------------------------------------------------------


    !-------------SOLUTION OF LINEAR SYSTEMS---------------------------------------
    DO i = 1, m_max_c

       mode = list_mode(i)

       !-------------SOURCES TERMS----------------------------------------------------
       tps = user_time(dummy)

       DO k = 1, 6
          rhs_H (:,k) = mu_H_field*(4*Hn(:,k,i)-Hn1(:,k,i))/(2*dt)
       END DO
       DO k = 1, 2
          rhs_phi (:,k) = mu_phi*(4*phin(:,k,i)-phin1(:,k,i))/(2*dt)
       END DO
       CALL courant(H_mesh,phi_mesh,INTERFACE,sigma,mu_phi,mu_H_field,time,mode, &
            rhs_H, rhs_phi, NL(:,:,i), f_H, f_phi)
       !Takes care of discontinuous mu
       IF (interface_H_mu%mes >0) THEN
          !write(*,*) 'TTTTTTTTTTTTTTTTTTTT courant_mu enleve'
          CALL  courant_mu(H_mesh,interface_H_mu,sigma,mu_H_field,time,mode,f_H,&
                           NL(:,:,i))
       END IF
       tps = user_time(dummy) - tps; tps_cumul=tps_cumul+tps
       !WRITE(*,*) ' Tps courant', tps
       !------------------------------------------------------------------------------

       !-------------INTERFACE INTEGRAL-----------------------------------------------
       tps = user_time(dummy)
       !write(*,*) 'TTTTTTTTTTTTTTTTTTTT surf_int enleve'
       CALL surf_int(H_mesh,phi_mesh,INTERFACE,interface_H_mu,sigma,mu_phi,mu_H_field,time,mode,f_H,f_phi, &
            R_fourier, index_fourier, rhs_phi)
       tps = user_time(dummy) - tps; tps_cumul=tps_cumul+tps
       !WRITE(*,*) ' Tps surf_int', tps
       !------------------------------------------------------------------------------

       !-------------FORMAT CHANGE----------------------------------------------------
       tps = user_time(dummy)
       fG1(:H_mesh%np)                = f_H(:,1)
       fG1(H_mesh%np+1:2*H_mesh%np)   = f_H(:,4)
       fG1(2*H_mesh%np+1:3*H_mesh%np) = f_H(:,5)
       fG1(3*H_mesh%np+1:)            = f_phi(:,1)

       fG2(:H_mesh%np)                = f_H(:,2)
       fG2(H_mesh%np+1:2*H_mesh%np)   = f_H(:,3)
       fG2(2*H_mesh%np+1:3*H_mesh%np) = f_H(:,6)
       fG2(3*H_mesh%np+1:)            = f_phi(:,2)

       solG1(:H_mesh%np)                = H_ext(:,1,i)
       solG1(H_mesh%np+1:2*H_mesh%np)   = H_ext(:,4,i)
       solG1(2*H_mesh%np+1:3*H_mesh%np) = H_ext(:,5,i)
       solG1(3*H_mesh%np+1:)            = 2*phin(:,1,i) - phin1(:,1,i)

       solG2(:H_mesh%np)                = H_ext(:,2,i)
       solG2(H_mesh%np+1:2*H_mesh%np)   = H_ext(:,3,i)
       solG2(2*H_mesh%np+1:3*H_mesh%np) = H_ext(:,6,i)
       solG2(3*H_mesh%np+1:)            = 2*phin(:,2,i) - phin1(:,2,i)
       tps = user_time(dummy) - tps; tps_cumul=tps_cumul+tps
       !WRITE(*,*) ' Tps change format', tps
       !-------------------------------------------------------------------------------

       !Jan 29 2007
       !Traitement des conditions aux limites periodiques
       IF (per) THEN
          CALL bc_per(H_phi_per, fG1)
          CALL bc_per(H_phi_per, fG2)
       END IF
       !Jan 29 2007

       !Feb 2 2007
       !Preconditionnement
       fG1 = fG1 * precond1(i)%DRL 
       fG2 = fG2 * precond2(i)%DRL 
       !Feb 2 2007

       !-------------DIRICHLET BOUNDARY CONDITIONS-------------------------------------
       tps = user_time(dummy)
       fG1(            js_D_H(1)%DIL) = Hexact(1, H_mesh%rr(:,js_D_H(1)%DIL), mode, mu_H_field, time)
       fG1(H_mesh%np  +js_D_H(2)%DIL) = Hexact(4, H_mesh%rr(:,js_D_H(2)%DIL), mode, mu_H_field, time)
       fG1(2*H_mesh%np+js_D_H(3)%DIL) = Hexact(5, H_mesh%rr(:,js_D_H(3)%DIL), mode, mu_H_field, time)
       fG1(3*H_mesh%np+js_D_phi%DIL)  = Phiexact(1, phi_mesh%rr(:,js_D_phi%DIL), mode, mu_phi, time)

       fG2(            js_D_H(1)%DIL) = Hexact(2, H_mesh%rr(:,js_D_H(1)%DIL), mode, mu_H_field, time)
       fG2(H_mesh%np  +js_D_H(2)%DIL) = Hexact(3, H_mesh%rr(:,js_D_H(2)%DIL), mode, mu_H_field, time)
       fG2(2*H_mesh%np+js_D_H(3)%DIL) = Hexact(6, H_mesh%rr(:,js_D_H(3)%DIL), mode, mu_H_field, time)
       fG2(3*H_mesh%np+js_D_phi%DIL)  = Phiexact(2, phi_mesh%rr(:,js_D_phi%DIL), mode, mu_phi, time)
       tps = user_time(dummy) - tps; tps_cumul=tps_cumul+tps
       !WRITE(*,*) ' Tps bcs', tps
       !------------------------------------------------------------------------------

       !-------------SOLVING THE LINEAR SYSTEMS----------------------------------------
       tps = user_time(dummy)
       CALL solve_splib(ipar_sp, fpar_sp, H_phi_mat1(i), fG1, solG1, isave(2*i-1), nb_tot_syst_lin)
       !            WRITE(*,*) 'Pb1: T precond ', fpar_sp(3),' T solving ', fpar_sp(4), 'for mode ', mode
       CALL solve_splib(ipar_sp, fpar_sp, H_phi_mat2(i), fG2, solG2, isave(2*i))
       !            WRITE(*,*) 'Pb2: T precond ', fpar_sp(3),' T solving ', fpar_sp(4), 'for mode ', mode
       tps = user_time(dummy) - tps; tps_cumul=tps_cumul+tps
       !WRITE(*,*) ' Tps solve Maxwell', tps
       !-------------------------------------------------------------------------------


       !-------------UPDATE------------------------------------------------------------
       tps = user_time(dummy)
       Hn1(:,:,i) = Hn(:,:,i)
       phin1(:,:,i) = phin(:,:,i)

       Hn(:,1,i)   = solG1(1:H_mesh%np) 
       Hn(:,4,i)   = solG1(H_mesh%np+1:2*H_mesh%np)
       Hn(:,5,i)   = solG1(2*H_mesh%np+1:3*H_mesh%np)
       phin(:,1,i) = solG1(3*H_mesh%np+1:3*H_mesh%np+phi_mesh%np)

       Hn(:,2,i)   = solG2(1:H_mesh%np) 
       Hn(:,3,i)   = solG2(H_mesh%np+1:2*H_mesh%np)
       Hn(:,6,i)   = solG2(2*H_mesh%np+1:3*H_mesh%np)
       phin(:,2,i) = solG2(3*H_mesh%np+1:3*H_mesh%np+phi_mesh%np)
       tps = user_time(dummy) - tps; tps_cumul=tps_cumul+tps
       !WRITE(*,*) ' Tps update', tps
       !------------------------------------------------------------------------------ 

    ENDDO
    isave= ABS(isave)
    tps_tot = user_time(dummy) - tps_tot
    WRITE(*,'(A,2(f13.3,2x))') ' Tps boucle en temps Maxwell', tps_tot, tps_cumul
    WRITE(*,*) ' TIME = ', time, '========================================'

  END SUBROUTINE maxwell_decouple


  SUBROUTINE  mat_maxwell(H_mesh, phi_mesh, INTERFACE, interface_H_mu, mode, ia,  &
       ja, aa1, aa2, mu_H_field, mu_phi, c_mass, stab, sigma, R_fourier, index_fourier)
    USE def_type_mesh
    USE Dir_nodes
    USE gauss_points
    USE boundary

    IMPLICIT NONE    
    TYPE(mesh_type),            INTENT(IN)    :: H_mesh
    TYPE(mesh_type),            INTENT(IN)    :: phi_mesh
    TYPE(interface_type),       INTENT(IN)    :: INTERFACE, interface_H_mu 
    INTEGER,                    INTENT(IN)    :: mode   
    INTEGER,      DIMENSION(:), INTENT(IN)    :: ia, ja
    REAL(KIND=8), DIMENSION(:), INTENT(INOUT) :: aa1, aa2
    REAL(KIND=8),               INTENT(IN)    :: mu_phi, c_mass
    REAL(KIND=8), DIMENSION(3), INTENT(IN)    :: stab
    REAL(KIND=8), DIMENSION(:), INTENT(IN)    :: sigma, mu_H_field
    REAL(KIND=8), OPTIONAL :: R_fourier
    INTEGER,      OPTIONAL :: index_fourier
    REAL(KIND=8), DIMENSION(phi_mesh%gauss%n_ws,phi_mesh%gauss%l_Gs)   :: w_cs
    REAL(KIND=8), DIMENSION(2,  H_mesh%gauss%n_w, phi_mesh%gauss%l_Gs, H_mesh%mes) :: dw_cs

    INTEGER :: m, l, ms, ls, ni, nj, k, h, k1, h1, i, j, i_b, j_b, H_bloc_size, p, &
         n_ws1, n_ws2, nj1, nj2, n_w2, n_w1, m1, m2, ki, kj,ib,jb, ms1, ms2
    REAL(KIND=8) :: x, y, z, hm1, stab_div, stab_colle_H_phi
    REAL(KIND=8) :: b, b1, b2, b3, b4, b5, b6, ray, error
    LOGICAL :: mark=.FALSE.

    !REAL(KIND=8), DIMENSION(H_mesh%gauss%n_w,H_mesh%gauss%n_w) :: aij, bij   
    REAL(KIND=8), DIMENSION(9,H_mesh%gauss%n_w,H_mesh%gauss%n_w)  :: TH
    REAL(KIND=8), DIMENSION(phi_mesh%gauss%n_w,phi_mesh%gauss%n_w):: TPhi

    !MATRICES POUR LES TERMES DE VOLUMES c_mass*mu_H*H + Rot((1/sigma)Rot(H)) - Grad(Div(H))
    !                                    -c_mass*mu_phi*Lap(Phi)
    !========================================================================
    !Le probleme est decouple en deux sous groupes de variables : 
    !H1, H4, H5 et Phi1 d'une part et H2, H3, H6 et Phi2 d'autre part.
    !Les matrices (symetriques sans terme de bord) s'ecrivent : 

    !MATRICE 1 :: 
    ! (------------------------) 
    ! ( TH1 | TH2 | TH3 |      )   H1
    ! (     | TH4 | TH5 |      )   H4
    ! (           | TH6 |      )   H5
    ! (                 | TPhi )   Phi1
    ! (------------------------)

    !MATRICE 2 (TH2 => TH8 et TH5 => TH9:: 
    ! (------------------------) 
    ! ( TH1 | TH8 | TH3 |      )   H2
    ! (     | TH4 | TH9 |      )   H3
    ! (           | TH6 |      )   H6
    ! (                 | TPhi )   Phi2
    ! (------------------------)
    !=========================================================================


    REAL(KIND=8), DIMENSION(9,H_mesh%gauss%n_w,H_mesh%gauss%n_w)      :: Hsij
    REAL(KIND=8), DIMENSION(phi_mesh%gauss%n_w,phi_mesh%gauss%n_w)    :: Phisij
    REAL(KIND=8), DIMENSION(6,phi_mesh%gauss%n_w,phi_mesh%gauss%n_w)  :: Sij

    ! MATRICES POUR LES TERMES DE BORDS Hsij et Phisij
    !=================================================
    ! (--------------------------------------------------------------------)
    ! ( Hsij(1)        |        Hsij(2) | Hsij(4)        || Sij(1)         )
    ! (        Hsij(1) | Hsij(3)        |        Hsij(4) ||        Sij(2)  )
    ! (--------------------------------------------------------------------)
    ! (                | Hsij(5)        |                ||        Sij(3)  )
    ! (                |        Hsij(5) |                || Sij(4)         )
    ! (--------------------------------------------------------------------)
    ! ( Hsij(7)        |        Hsij(9) | Hsij(6)        || Sij(5)         )             
    ! (        Hsij(7) | Hsij(8)        |        Hsij(6) ||        Sij(6)  ) 
    ! (====================================================================)
    ! ( Sij'(1)        |        Sij'(3) | Sij'(5)        || Phisij         )
    ! (        Sij'(2) | Sij'(4)        |        Sij'(6) ||        Phisij  )
    ! (------------------------------------------------------------------- )
    !
    ! L'autre partie des termes croises est la symetrique de la premiere
    ! juste apres le calcul du terme de bord dissymetrique    

    REAL(KIND=8), DIMENSION(phi_mesh%gauss%n_w,phi_mesh%gauss%n_w) :: b7ij

    !fonctions de forme propres a H_mesh
    REAL(KIND=8), DIMENSION(:,:),     POINTER :: ww_H
    !derivees des fonctions de forme propres a H_mesh
    REAL(KIND=8), DIMENSION(:,:,:,:), POINTER :: dw_H
    !jacobien pour H
    REAL(KIND=8), DIMENSION(:,:),     POINTER :: rj_H
    !fonctions de forme propres a phi_mesh
    REAL(KIND=8), DIMENSION(:,:),     POINTER :: ww_phi  
    !derivees des fonctions de forme propres a phi_mesh
    REAL(KIND=8), DIMENSION(:,:,:,:), POINTER :: dw_phi
    !jacobien pour phi
    REAL(KIND=8), DIMENSION(:,:),     POINTER :: rj_phi   

    REAL(KIND=8), DIMENSION(2,phi_mesh%gauss%l_Gs) :: gauss1, gauss2
    INTEGER  :: ls1, ls2, n, n1, n2, ni1, ni2, nip
    REAL(KIND=8) :: ref, diff, d1, d2, mu_H, c_mu_H, c_mu_phi, h2, muhl, dzmuhl, drmuhl, c_div
    !June 8 2008
    REAL(KIND=8) :: c_sym=0.d0 ! Symmetrization of the bilinear form
    !June 8 2008
    !REAL(KIND=8), DIMENSION(6) :: lap
    !REAL(KIND=8), DIMENSION(2,6) :: a

    c_mu_phi = c_mass*mu_phi 

    ww_H   => H_mesh%gauss%ww
    dw_H   => H_mesh%gauss%dw
    rj_H   => H_mesh%gauss%rj 
    ww_phi => phi_mesh%gauss%ww
    dw_phi => phi_mesh%gauss%dw
    rj_phi => phi_mesh%gauss%rj 


    stab_div = stab(1) 
    stab_colle_H_phi = stab(2) 

    !==Block on H
    DO m = 1, H_mesh%me

       TH = 0.d0

       DO l = 1, H_mesh%gauss%l_G

          !--------On calcule le rayon du point gauss
          !Feb 8 2007, muhl
          muhl = SUM(mu_H_field(H_mesh%jj(:,m))*ww_H(:,l))
          drmuhl = SUM(mu_H_field(H_mesh%jj(:,m))*dw_H(1,:,l,m))
          dzmuhl = SUM(mu_H_field(H_mesh%jj(:,m))*dw_H(2,:,l,m))
          c_mu_H   = c_mass*muhl
          !Feb 8 2007, muhl
          !June 7 2008, normalization
          c_div    = stab_div/muhl**2
          !June 7 2008, normalization

          ray = 0
          DO ni = 1, H_mesh%gauss%n_w;  i = H_mesh%jj(ni,m)
             ray = ray + H_mesh%rr(1,i)*ww_H(ni,l)
          END DO

          DO ni = 1, H_mesh%gauss%n_w     
             DO nj = 1, H_mesh%gauss%n_w

                ! mu_H * <bi,bj> + <Div bi,Div bj> + <(1/sigma) Rot bi,Rot bj>
                TH(1,ni,nj) = TH(1,ni,nj) +  rj_H(l,m) * ray* ( &
                     c_mu_H*ww_H(ni,l)*ww_H(nj,l) & 
                     + (dw_H(2,ni,l,m)*dw_H(2,nj,l,m) + mode**2/ray**2*ww_H(ni,l)*ww_H(nj,l))/sigma(m) &
                                !DIVERGENCE, June 8 2008
                     + c_div*(muhl*(ww_H(ni,l)/ray+dw_H(1,ni,l,m)) + ww_H(ni,l)*drmuhl) &
                     *(muhl*(ww_H(nj,l)/ray+dw_H(1,nj,l,m)) + ww_H(nj,l)*drmuhl))
                !+ stab_div*(ww_H(ni,l)*ww_H(nj,l)/ray**2+dw_H(1,ni,l,m)*dw_H(1,nj,l,m) &
                !+ 1/ray*(ww_H(ni,l)*dw_H(1,nj,l,m)+ww_H(nj,l)*dw_H(1,ni,l,m))))
                !                       

                TH(2,ni,nj) = TH(2,ni,nj)+ rj_H(l,m) * ray* (  &
                     mode/ray**2 * ww_H(ni,l)*(ww_H(nj,l)+ray*dw_H(1,nj,l,m))/sigma(m) &
                                !DIVERGENCE , June 8 2008
                     + c_div*mode/ray*(muhl*(ww_H(ni,l)/ray+dw_H(1,ni,l,m)) + ww_H(ni,l)*drmuhl)*muhl*ww_H(nj,l))
                !+ stab_div*mode/ray*(ww_H(ni,l)/ray+dw_H(1,ni,l,m))*ww_H(nj,l))
                !             

                TH(8,ni,nj) = TH(8,ni,nj)+ rj_H(l,m) * ray* (  &
                     - mode/ray**2 * ww_H(ni,l)*(ww_H(nj,l)+ray*dw_H(1,nj,l,m))/sigma(m) &
                                !DIVERGENCE, June 8 2008
                     - c_div*mode/ray*(muhl*(ww_H(ni,l)/ray+dw_H(1,ni,l,m)) + ww_H(ni,l)*drmuhl)*muhl*ww_H(nj,l))
                !-stab_div*mode/ray*(ww_H(ni,l)/ray+dw_H(1,ni,l,m))*ww_H(nj,l))
                !           

                TH(3,ni,nj) = TH(3,ni,nj)+ rj_H(l,m) * ray* ( &
                     - dw_H(2,ni,l,m)*dw_H(1,nj,l,m)/sigma(m) &
                                !DIVERGENCE, June 8 2008
                     + c_div*(muhl*(ww_H(ni,l)/ray+dw_H(1,ni,l,m)) + ww_H(ni,l)*drmuhl)*&
                     (muhl*dw_H(2,nj,l,m) + ww_H(nj,l)*dzmuhl))
                !+ stab_div*(ww_H(ni,l)/ray+dw_H(1,ni,l,m))*dw_H(2,nj,l,m))
                !        

                TH(4,ni,nj) = TH(4,ni,nj) + rj_H(l,m) * ray* ( &
                     c_mu_H*ww_H(ni,l)*ww_H(nj,l)  & 
                     + (dw_H(2,ni,l,m)*dw_H(2,nj,l,m) &
                     + 1/ray**2 *(ww_H(ni,l)+ray*dw_H(1,ni,l,m))*(ww_H(nj,l)+ray*dw_H(1,nj,l,m)))/sigma(m) &
                                !DIVERGENCE, June 8 2008
                     +c_div*muhl**2*mode**2/ray**2*ww_H(ni,l)*ww_H(nj,l))
                !+stab_div*mode**2/ray**2*ww_H(ni,l)*ww_H(nj,l)) 
                !

                TH(5,ni,nj) = TH(5,ni,nj)  + rj_H(l,m) * ray* (&
                     + mode/ray*dw_H(2,ni,l,m)*ww_H(nj,l)/sigma(m) &
                                !DIVERGENCE, June 8 2008
                     +c_div*mode/ray*muhl*ww_H(ni,l)*(muhl*dw_H(2,nj,l,m) + ww_H(nj,l)*dzmuhl))
                !+stab_div*mode/ray*ww_H(ni,l)*dw_H(2,nj,l,m))
                !    

                TH(9,ni,nj) = TH(9,ni,nj)  + rj_H(l,m) * ray* (&
                     - mode/ray*dw_H(2,ni,l,m)*ww_H(nj,l)/sigma(m) &
                                !DIVERGENCE, June 8 2008
                     - c_div*mode/ray*muhl*ww_H(ni,l)*(muhl*dw_H(2,nj,l,m) + ww_H(nj,l)*dzmuhl))
                !- stab_div*mode/ray*ww_H(ni,l)*dw_H(2,nj,l,m))             
                !

                TH(6,ni,nj) = TH(6,ni,nj) + rj_H(l,m) * ray* ( &
                     c_mu_H*ww_H(ni,l)*ww_H(nj,l)  &
                     + (mode**2/ray**2*ww_H(ni,l)*ww_H(nj,l) + dw_H(1,ni,l,m)*dw_H(1,nj,l,m))/sigma(m) &
                                !DIVERGENCE, June 8 2008
                     + c_div*(muhl*dw_H(2,ni,l,m) + ww_H(ni,l)*dzmuhl) &
                     *(muhl*dw_H(2,nj,l,m) + ww_H(nj,l)*dzmuhl))
                !+ stab_div*dw_H(2,ni,l,m)*dw_H(2,nj,l,m))
                !               
             ENDDO
          END DO

       END DO

       DO ki= 1, 3  
          DO ni = 1, H_mesh%gauss%n_w 
             i = H_mesh%jj(ni, m)
             ib = i + (ki-1)*H_mesh%np 
             DO kj = 1, 3
                DO nj = 1, H_mesh%gauss%n_w
                   j = H_mesh%jj(nj, m)
                   jb = j + (kj-1)*H_mesh%np 
                   DO p = ia(ib),  ia(ib+1) - 1
                      IF (ja(p) == jb) THEN

                         IF   ((ki == 1) .AND. (kj == 1)) THEN
                            aa1(p) = aa1(p) + TH(1,ni,nj)
                            aa2(p) = aa2(p) + TH(1,ni,nj)

                         ELSEIF   ((ki == 1) .AND. (kj == 2)) THEN
                            aa1(p) = aa1(p) + TH(2,ni,nj)
                            aa2(p) = aa2(p) + TH(8,ni,nj)

                         ELSEIF   ((ki == 2) .AND. (kj == 1)) THEN  
                            aa1(p) = aa1(p) + TH(2,nj,ni)
                            aa2(p) = aa2(p) + TH(8,nj,ni)

                         ELSEIF  ((ki == 1) .AND. (kj == 3)) THEN
                            aa1(p) = aa1(p) + TH(3,ni,nj)
                            aa2(p) = aa2(p) + TH(3,ni,nj)                          

                         ELSEIF  ((ki == 3) .AND. (kj == 1)) THEN
                            aa1(p) = aa1(p) + TH(3,nj,ni)  
                            aa2(p) = aa2(p) + TH(3,nj,ni)  

                         ELSEIF  ((ki == 2) .AND. (kj == 2)) THEN
                            aa1(p) = aa1(p) + TH(4,ni,nj)
                            aa2(p) = aa2(p) + TH(4,ni,nj)

                         ELSEIF   ((ki == 2) .AND. (kj == 3)) THEN
                            aa1(p) = aa1(p) + TH(5,ni,nj)
                            aa2(p) = aa2(p) + TH(9,ni,nj)

                         ELSEIF   ((ki == 3) .AND. (kj == 2)) THEN
                            aa1(p) = aa1(p) + TH(5,nj,ni)   
                            aa2(p) = aa2(p) + TH(9,nj,ni)

                         ELSEIF  ((ki == 3) .AND. (kj == 3)) THEN  
                            aa1(p) = aa1(p) + TH(6,ni,nj)   
                            aa2(p) = aa2(p) + TH(6,ni,nj) 

                         ENDIF

                      ENDIF
                   END DO
                END DO
             END DO
          END DO
       END DO

    END DO


    !==Block on phi  

    DO m = 1,phi_mesh%me


       TPhi = 0.d0

       DO l = 1, phi_mesh%gauss%l_G

          !--------On calcule le rayon du point gauss
          ray = 0
          DO ni = 1, phi_mesh%gauss%n_w;  i = phi_mesh%jj(ni,m)
             ray = ray + phi_mesh%rr(1,i)*ww_phi(ni,l)
          END DO

          DO ni = 1, phi_mesh%gauss%n_w
             DO nj = 1, phi_mesh%gauss%n_w

                !mu_phi * <Grad bi, Grad bj>
                TPhi(ni,nj) = TPhi(ni,nj) + rj_phi(l,m) * ray* (c_mu_phi) & 
                     *(dw_phi(1,ni,l,m)*dw_phi(1,nj,l,m)+dw_phi(2,ni,l,m)*dw_phi(2,nj,l,m) &
                     +mode**2/ray**2*ww_phi(ni,l)*ww_phi(nj,l))
             ENDDO
          END DO

       END DO

       !TEST      
       !TPhi = 0.d0
       !TEST

       DO ni = 1, phi_mesh%gauss%n_w 
          i = phi_mesh%jj(ni, m)
          ib = i + 3*H_mesh%np 
          DO nj = 1, phi_mesh%gauss%n_w
             j = phi_mesh%jj(nj, m)
             jb = j + 3*H_mesh%np  
             DO p = ia(ib),  ia(ib+1) - 1
                IF (ja(p) == jb) THEN
                   aa1(p) = aa1(p) + TPhi(ni,nj)
                   aa2(p) = aa2(p) + TPhi(ni,nj)
                   EXIT
                ENDIF

             END DO
          END DO
       END DO

    END DO

    !*********************************************************************************
    !--------------------TERMES SUR L'INTERFACE SIGMA----------------------------------
    !**********************************************************************************

    WRITE(*,*) 'Assemblage interface'
    CALL gauss(phi_mesh)
    n_ws1 = H_mesh%gauss%n_ws
    n_ws2 = phi_mesh%gauss%n_ws
    n_w1  = H_mesh%gauss%n_w
    n_w2  = phi_mesh%gauss%n_w

    IF (H_mesh%gauss%n_ws == n_ws) THEN

       DO ms = 1, interface%mes

          ms2 = interface%mesh2(ms)
          m2 = phi_mesh%neighs(ms2)
          ms1 = interface%mesh1(ms)
          m1 = H_mesh%neighs(ms1)

          ref = 1.d-8+SUM((H_mesh%rr(:,H_mesh%jjs(1,ms1)) - H_mesh%rr(:,H_mesh%jjs(2,ms1)))**2)
          diff = SUM((H_mesh%rr(:,H_mesh%jjs(1,ms1)) - phi_mesh%rr(:,phi_mesh%jjs(1,ms2)))**2)
          IF (diff/ref .LT. 1.d-10) THEN ! 1 = 1
             w_cs = wws
          ELSE                ! 1 = 2
             DO ls = 1, l_Gs
                w_cs(1,ls)= wws(2,ls)
                w_cs(2,ls)= wws(1,ls)
                w_cs(3,ls)= wws(3,ls) 
                WRITE(*,*) ' Ouaps! oder of shape functions changed?'
             END DO
          END IF

          DO ls = 1, l_Gs
             gauss2(1,ls) = SUM(phi_mesh%rr(1,phi_mesh%jjs(:,ms2))*phi_mesh%gauss%wws(:,ls))
             gauss2(2,ls) = SUM(phi_mesh%rr(2,phi_mesh%jjs(:,ms2))*phi_mesh%gauss%wws(:,ls))
             gauss1(1,ls) = SUM(  H_mesh%rr(1,  H_mesh%jjs(:,ms1))*  H_mesh%gauss%wws(:,ls))
             gauss1(2,ls) = SUM(  H_mesh%rr(2,  H_mesh%jjs(:,ms1))*  H_mesh%gauss%wws(:,ls))
          END DO

          DO ls2 = 1, l_Gs
             ref = SQRT(1.d-8+SUM(gauss2(:,ls2)**2))
             mark = .FALSE.
             DO ls1 = 1, l_Gs
                diff = SQRT(SUM((gauss2(:,ls2)-gauss1(:,ls1))**2))
                IF (diff .LT. 1.d-10) THEN
                   dw_cs(:,:,ls2,ms1) =  H_mesh%gauss%dw_s(:,:,ls1,ms1)
                   mark = .TRUE.
                   EXIT
                END IF
             END DO
             IF (.NOT.mark) WRITE(*,*) ' BUG '
          END DO

       END DO

    ELSE
       DO ms = 1, interface%mes

          ms2 = interface%mesh2(ms)
          m2 = phi_mesh%neighs(ms2)
          ms1 = interface%mesh1(ms)
          m1 = H_mesh%neighs(ms1)

          ref = 1.d-8+SUM((H_mesh%rr(:,H_mesh%jjs(1,ms1)) - H_mesh%rr(:,H_mesh%jjs(2,ms1)))**2)
          diff = SUM((H_mesh%rr(:,H_mesh%jjs(1,ms1)) - phi_mesh%rr(:,phi_mesh%jjs(1,ms2)))**2)
          IF (diff/ref .LT. 1.d-10) THEN ! 1 = 1
             DO ls = 1, l_Gs
                w_cs(1,ls)= wws(1,ls)+0.5*wws(3,ls)
                w_cs(2,ls)= wws(2,ls)+0.5*wws(3,ls)
                w_cs(3,ls)= 0  
             END DO
          ELSE                ! 1 = 2
             DO ls = 1, l_Gs
                w_cs(1,ls)= wws(2,ls)+0.5*wws(3,ls)
                w_cs(2,ls)= wws(1,ls)+0.5*wws(3,ls)
                w_cs(3,ls)= 0 
                WRITE(*,*) ' Ouaps! oder of shape functions changed?'
             END DO
          END IF

          DO ls = 1, l_Gs
             dw_cs(1,:,ls,ms1) = H_mesh%gauss%dw(1,:,1,m1)
             dw_cs(2,:,ls,ms1) = H_mesh%gauss%dw(2,:,1,m1)
          END DO

       END DO
    END IF

    error = 0
    DO ms = 1, interface%mes

       ms2 = interface%mesh2(ms)
       ms1 = interface%mesh1(ms)
       m2 = phi_mesh%neighs(ms2)
       m1 =   H_mesh%neighs(ms1)
       mu_H = SUM(mu_H_field(H_mesh%jj(:,m1)))/H_mesh%gauss%n_w
       hm1 = stab_colle_H_phi*(((mu_phi+mu_H)/mu_H)/SUM(rjs(:,ms2)))
 
       !====================================================================================
       !------------------------------------TERMES SUR LE BLOC H----------------------------
       !====================================================================================

       !-------------------------------hm1 (bi x ni) . (bj x nj)----------------------------
       !====================================================================================

       Hsij = 0.d0  
       DO ls = 1, l_Gs
          !--------On calcule le rayon du point gauss
          ray = SUM(phi_mesh%rr(1,phi_mesh%jjs(:,ms2))* phi_mesh%gauss%wws(:,ls))
          x = hm1*rjs(ls,ms2)*ray 

          DO ni = 1, n_ws1
             DO nj = 1, n_ws1
                y = x * w_cs(ni,ls)*w_cs(nj,ls)
                Hsij(1,ni,nj) = Hsij(1,ni,nj) + y*(rnorms(2,ls,ms2)**2) 
                Hsij(4,ni,nj) = Hsij(4,ni,nj) - y*rnorms(1,ls,ms2)*rnorms(2,ls,ms2)                        
                Hsij(5,ni,nj) = Hsij(5,ni,nj) + y                                                
                Hsij(6,ni,nj) = Hsij(6,ni,nj) + y*(rnorms(1,ls,ms2)**2)
             ENDDO
          ENDDO

       ENDDO


       !TEST
       !Hsij = 0.d0
       !Hsij = Hsij / hm1
       !TEST

       DO ki= 1, 3 
          DO ni = 1, n_ws1 
             i = interface%jjs1(ni,ms)
             ib = i + (ki-1)*H_mesh%np 
             DO kj = 1, 3
                DO nj = 1, n_ws1
                   j = interface%jjs1(nj,ms)
                   jb = j + (kj-1)*H_mesh%np 
                   DO p = ia(ib),  ia(ib+1) - 1
                      IF (ja(p) == jb) THEN
                         IF  ((ki == 1) .AND. (kj == 1)) THEN
                            aa1(p) = aa1(p) + Hsij(1,ni,nj)
                            aa2(p) = aa2(p) + Hsij(1,ni,nj)
                            EXIT                                                                                 
                         ELSEIF  ((ki == 1) .AND. (kj == 3)) THEN
                            aa1(p) = aa1(p) + Hsij(4,ni,nj)
                            aa2(p) = aa2(p) + Hsij(4,ni,nj)
                            EXIT                                                                                
                         ELSEIF  ((ki == 3) .AND. (kj == 1)) THEN
                            aa1(p) = aa1(p) + Hsij(4,nj,ni)
                            aa2(p) = aa2(p) + Hsij(4,nj,ni)
                            EXIT                                                                                 
                         ELSEIF  ((ki == 2) .AND. (kj == 2)) THEN
                            aa1(p) = aa1(p) + Hsij(5,ni,nj)
                            aa2(p) = aa2(p) + Hsij(5,ni,nj)
                            EXIT                                                                                
                         ELSEIF  ((ki == 3) .AND. (kj == 3)) THEN
                            aa1(p) = aa1(p) + Hsij(6,ni,nj)
                            aa2(p) = aa2(p) + Hsij(6,ni,nj)
                            EXIT 
                         ENDIF

                      ENDIF
                   END DO
                END DO
             END DO
          END DO
       END DO

       !====================================================================================
       !------------------------(1/sigma) (Rot bj) . (bi x ni)------------------------------
       !====================================================================================


       Hsij = 0.d0

       DO ls = 1, phi_mesh%gauss%l_Gs

          !--------On calcule le rayon du point gauss
          ray = SUM(phi_mesh%rr(1,phi_mesh%jjs(:,ms2))* phi_mesh%gauss%wws(:,ls))
          x = rjs(ls,ms2)*ray/sigma(m1)

          !terme sans derivees
          DO ni = 1,n_ws1
             DO nj = 1, n_ws1
                y = x*w_cs(ni,ls)*w_cs(nj,ls)
                Hsij(2,ni,nj) = Hsij(2,ni,nj) + y * (-mode/ray)*(-rnorms(1,ls,ms2))
                Hsij(3,ni,nj) = Hsij(3,ni,nj) + y *   mode/ray *(-rnorms(1,ls,ms2))
                Hsij(5,ni,nj) = Hsij(5,ni,nj) + y * (-1/ray)   *(-rnorms(1,ls,ms2))
                Hsij(8,ni,nj) = Hsij(8,ni,nj) + y * (-mode/ray)*(-rnorms(2,ls,ms2))
                Hsij(9,ni,nj) = Hsij(9,ni,nj) + y *   mode/ray *(-rnorms(2,ls,ms2))
             ENDDO
          ENDDO

       ENDDO

       !TEST
       !Hsij = 0.d0
       !TEST


       DO ki= 1, 3  
          DO ni = 1, n_ws1
             i = interface%jjs1(ni,ms)
             ib = i + (ki-1)*H_mesh%np 
             DO kj = 1, 3
                DO nj = 1, n_ws1
                   j = interface%jjs1(nj,ms)
                   jb = j + (kj-1)*H_mesh%np 
                   DO p = ia(ib),  ia(ib+1) - 1
                      IF (ja(p) == jb) THEN

                         IF  ( (ki == 2) .AND. (kj == 1)) THEN
                            aa1(p) = aa1(p) + Hsij(2,ni,nj)
                            aa2(p) = aa2(p) + Hsij(3,ni,nj)
                            EXIT
                         ELSEIF  ((ki == 2) .AND. (kj == 2))  THEN  
                            aa1(p) = aa1(p) + Hsij(5,ni,nj)
                            aa2(p) = aa2(p) + Hsij(5,ni,nj)
                            EXIT
                         ELSEIF  ( (ki == 2) .AND. (kj == 3)) THEN
                            aa1(p) = aa1(p) + Hsij(8,ni,nj)
                            aa2(p) = aa2(p) + Hsij(9,ni,nj) 
                            EXIT
                         ENDIF
                      ENDIF
                   END DO

                END DO
             END DO
          END DO
       END DO

       !Feb 2 2007
       Hsij=c_sym*Hsij !SYM
       DO ki= 1, 3
          DO ni = 1, n_ws1
             i = interface%jjs1(ni,ms)
             ib = i + (ki-1)*H_mesh%np
             DO kj = 1, 3
                DO nj = 1, n_ws1
                   j = interface%jjs1(nj,ms)
                   jb = j + (kj-1)*H_mesh%np
                   DO p = ia(ib),  ia(ib+1) - 1
                      IF (ja(p) == jb) THEN

                         IF  ( (kj == 2) .AND. (ki == 1)) THEN
                            aa1(p) = aa1(p) + Hsij(2,nj,ni)
                            aa2(p) = aa2(p) + Hsij(3,nj,ni)
                            EXIT
                         ELSEIF  ((kj == 2) .AND. (ki == 2))  THEN
                            aa1(p) = aa1(p) + Hsij(5,nj,ni)
                            aa2(p) = aa2(p) + Hsij(5,nj,ni)
                            EXIT
                         ELSEIF  ( (kj == 2) .AND. (ki == 3)) THEN
                            aa1(p) = aa1(p) + Hsij(8,nj,ni)
                            aa2(p) = aa2(p) + Hsij(9,nj,ni)
                            EXIT
                         ENDIF
                      ENDIF
                   END DO

                END DO
             END DO
          END DO
       END DO
       !feb 2 2007


       Hsij = 0.d0

       DO ls = 1, phi_mesh%gauss%l_Gs

          !--------On calcule le rayon du point gauss
          ray = SUM(phi_mesh%rr(1,phi_mesh%jjs(:,ms2))* phi_mesh%gauss%wws(:,ls))
          x = rjs(ls,ms2)*ray /sigma(m1)

          !termes avec derivees
          DO ni = 1,n_ws1
             y = x*w_cs(ni,ls)
             DO nj = 1, n_w1
                Hsij(1,ni,nj) = Hsij(1,ni,nj) + y*(-dw_cs(2,nj,ls,ms1))*(-rnorms(2,ls,ms2))
                Hsij(4,ni,nj) = Hsij(4,ni,nj) + y*  dw_cs(1,nj,ls,ms1) *(-rnorms(2,ls,ms2))
                Hsij(5,ni,nj) = Hsij(5,ni,nj) + &
                     y*(-dw_cs(2,nj,ls,ms1)*(-rnorms(2,ls,ms2))-dw_cs(1,nj,ls,ms1)*(-rnorms(1,ls,ms2)))
                Hsij(6,ni,nj) = Hsij(6,ni,nj) + y*(-dw_cs(1,nj,ls,ms1))*(-rnorms(1,ls,ms2))
                Hsij(7,ni,nj) = Hsij(7,ni,nj) + y*  dw_cs(2,nj,ls,ms1) *(-rnorms(1,ls,ms2))
             ENDDO
          ENDDO

       ENDDO

       !TEST
       !Hsij = 0.d0
       !TEST

       DO ki= 1, 3  
          DO ni = 1, n_ws1
             i = interface%jjs1(ni,ms)
             !i = H_mesh%jjs(ni,ms1)
             ib = i + (ki-1)*H_mesh%np 
             DO kj = 1, 3
                DO nj = 1, n_w1
                   j = H_mesh%jj(nj,m1)
                   jb = j + (kj-1)*H_mesh%np

                   DO p = ia(ib),  ia(ib+1) - 1
                      IF (ja(p) == jb) THEN

                         IF  ((ki == 1) .AND. (kj == 1))  THEN
                            aa1(p) = aa1(p) + Hsij(1,ni,nj)
                            aa2(p) = aa2(p) + Hsij(1,ni,nj)
                            EXIT
                         ELSEIF  ((ki == 1) .AND. (kj == 3)) THEN
                            aa1(p) = aa1(p) + Hsij(4,ni,nj)
                            aa2(p) = aa2(p) + Hsij(4,ni,nj)
                            EXIT
                         ELSEIF  ((ki == 2) .AND. (kj == 2))  THEN  
                            aa1(p) = aa1(p) + Hsij(5,ni,nj)
                            aa2(p) = aa2(p) + Hsij(5,ni,nj)
                            EXIT
                         ELSEIF  ((ki == 3) .AND. (kj == 3)) THEN  
                            aa1(p) = aa1(p) + Hsij(6,ni,nj)   
                            aa2(p) = aa2(p) + Hsij(6,ni,nj) 
                            EXIT
                         ELSEIF  ((ki == 3) .AND. (kj == 1)) THEN
                            aa1(p) = aa1(p) + Hsij(7,ni,nj)
                            aa2(p) = aa2(p) + Hsij(7,ni,nj)
                            EXIT
                         ENDIF

                      ENDIF
                   END DO
                END DO
             END DO
          END DO
       END DO

       !Feb 2 2007
       Hsij=c_sym*Hsij !SYM
       DO ki = 1, 3
          DO ni = 1, n_w1
             i = H_mesh%jj(ni,m1)
             ib = i + (ki-1)*H_mesh%np
             DO kj= 1, 3
                DO nj = 1, n_ws1
                   j = interface%jjs1(nj,ms)
                   jb = j + (kj-1)*H_mesh%np
                   DO p = ia(ib),  ia(ib+1) - 1
                      IF (ja(p) == jb) THEN

                         IF  ((kj == 1) .AND. (ki == 1))  THEN
                            aa1(p) = aa1(p) + Hsij(1,nj,ni)
                            aa2(p) = aa2(p) + Hsij(1,nj,ni)
                            EXIT
                         ELSEIF  ((kj == 1) .AND. (ki == 3)) THEN
                            aa1(p) = aa1(p) + Hsij(4,nj,ni)
                            aa2(p) = aa2(p) + Hsij(4,nj,ni)
                            EXIT
                         ELSEIF  ((kj == 2) .AND. (ki == 2))  THEN
                            aa1(p) = aa1(p) + Hsij(5,nj,ni)
                            aa2(p) = aa2(p) + Hsij(5,nj,ni)
                            EXIT
                         ELSEIF  ((kj == 3) .AND. (ki == 3)) THEN
                            aa1(p) = aa1(p) + Hsij(6,nj,ni)
                            aa2(p) = aa2(p) + Hsij(6,nj,ni)
                            EXIT
                         ELSEIF  ((kj == 3) .AND. (ki == 1)) THEN
                            aa1(p) = aa1(p) + Hsij(7,nj,ni)
                            aa2(p) = aa2(p) + Hsij(7,nj,ni)
                            EXIT
                         ENDIF

                      ENDIF
                   END DO
                END DO
             END DO
          END DO
       END DO
       !Feb 2 2007


       !====================================================================================
       !------------------------------------TERMES SUR LE BLOC PHI--------------------------
       !====================================================================================

       !------------------------hm1 (Grad(phi_i) x ni).(Grad(phi_j) x nj)-------------------
       !====================================================================================

       Phisij = 0.d0

       DO ls = 1, phi_mesh%gauss%l_Gs

          !--------On calcule le rayon du point gauss
          ray = SUM(phi_mesh%rr(1,phi_mesh%jjs(:,ms2))* phi_mesh%gauss%wws(:,ls))
          x = hm1*rjs(ls,ms2)*ray 

          !terme sans derivee
          DO ni=1, n_ws2
             DO nj=1, n_ws2
                Phisij(ni,nj) = Phisij(ni,nj) + x*mode**2/ray**2*wws(ni,ls)*wws(nj,ls) 
             ENDDO
          ENDDO

       ENDDO

       !TEST
       !Phisij = 0.d0
       !Phisij = Phisij/hm1
       !TEST

       DO ni = 1, n_ws2
          i =  interface%jjs2(ni,ms)
          ib = i + 3*H_mesh%np 
          DO nj = 1, n_ws2
             j = interface%jjs2(nj,ms)
             jb = j + 3*H_mesh%np 
             DO p = ia(ib),  ia(ib+1) - 1
                IF (ja(p) == jb) THEN                   
                   aa1(p) = aa1(p) + Phisij(ni,nj)  
                   aa2(p) = aa2(p) + Phisij(ni,nj) 
                   EXIT
                ENDIF
             END DO
          END DO
       END DO

       Phisij = 0.d0

       DO ls = 1, l_Gs

          !--------On calcule le rayon du point gauss
          ray = SUM(phi_mesh%rr(1,phi_mesh%jjs(:,ms2))* phi_mesh%gauss%wws(:,ls))
          x = hm1*rjs(ls,ms2)*ray 

          !terme avec derivee
          DO ni = 1, n_w2
             DO nj = 1, n_w2
                Phisij(ni,nj) = Phisij(ni,nj) + x*( &
                     (dw_s(2,ni,ls,ms2)*rnorms(1,ls,ms2) - dw_s(1,ni,ls,ms2)*rnorms(2,ls,ms2))* &
                     (dw_s(2,nj,ls,ms2)*rnorms(1,ls,ms2) - dw_s(1,nj,ls,ms2)*rnorms(2,ls,ms2)))
             ENDDO
          ENDDO
       ENDDO

       !Phisij = 0.d0
       !Phisij = Phisij/hm1
       !TEST

       DO ni = 1, n_w2
          i = phi_mesh%jj(ni, m2)
          ib = i + 3*H_mesh%np 
          DO nj = 1, n_w2
             j = phi_mesh%jj(nj, m2)
             jb = j + 3*H_mesh%np 
             DO p = ia(ib),  ia(ib+1) - 1
                IF (ja(p) == jb) THEN
                   aa1(p) = aa1(p) + Phisij(ni,nj)
                   aa2(p) = aa2(p) + Phisij(ni,nj)
                   EXIT
                ENDIF
             END DO
          END DO
       END DO

       !====================================================================================   
       !------------------------------------TERMES CROISES----------------------------------
       !====================================================================================

       !====================================================================================
       !------------------------hm1 (bi x ni) . (Grad(phi_j) x nj)--------------------------
       !------------------      + hm1(Grad(phi_i) x ni).(bj x nj)---------------------------
       !====================================================================================

       Sij = 0.d0

       DO ls = 1, l_Gs

          !--------On calcule le rayon du point gauss
          ray = SUM(phi_mesh%rr(1,phi_mesh%jjs(:,ms2))* phi_mesh%gauss%wws(:,ls))
          x = hm1*rjs(ls,ms2)*ray 

          !terme sans derivee
          DO ni = 1, n_ws1
             DO nj = 1, n_ws2
                Sij(3,ni,nj) = Sij(3,ni,nj) + x*(mode/ray)*w_cs(ni,ls)*wws(nj,ls)
             ENDDO
          ENDDO
       ENDDO
       Sij(4,:,:) = -Sij(3,:,:)

       !TEST
       !Sij = 0.d0
       !Sij = Sij /hm1
       !TEST

       ki = 2
       DO ni = 1, n_ws1 
          i = interface%jjs1(ni,ms)
          ib = i + (ki-1)*H_mesh%np 
          DO nj = 1, n_ws2
             j = interface%jjs2(nj,ms)
             jb = j + 3*H_mesh%np 
             DO p = ia(ib),  ia(ib+1) - 1
                IF (ja(p) == jb) THEN
                   aa1(p) = aa1(p) + Sij(3,ni,nj) 
                   aa2(p) = aa2(p) + Sij(4,ni,nj)  
                   EXIT
                ENDIF
             END DO
          END DO
       ENDDO

       !TEST SYM
       !Feb 2 2003
       !Sij = 0.d0
       kj = 2
       DO ni = 1, n_ws2
          i = interface%jjs2(ni,ms)
          ib = i + 3*H_mesh%np
          DO nj = 1, n_ws1
             j = interface%jjs1(nj,ms)
             jb = j + (kj-1)*H_mesh%np
             DO p = ia(ib),  ia(ib+1) - 1
                IF (ja(p) == jb) THEN
                   aa1(p) = aa1(p) + Sij(3,nj,ni)
                   aa2(p) = aa2(p) + Sij(4,nj,ni)
                   EXIT
                ENDIF
             END DO
          END DO
       ENDDO
       !Feb 2 2003
       !TEST SYM

       Sij = 0.d0

       DO ls = 1, l_Gs

          !--------On calcule le rayon du point gauss
          ray = SUM(phi_mesh%rr(1,phi_mesh%jjs(:,ms2))* phi_mesh%gauss%wws(:,ls))
          x = hm1*rjs(ls,ms2)*ray 

          !terme avec derivee
          DO ni = 1, n_ws1
             y = x * w_cs(ni,ls)
             DO nj = 1, n_w2
                Sij(1,ni,nj) = Sij(1,ni,nj) + &
                     y*(-dw_s(1,nj,ls,ms2)*rnorms(2,ls,ms2)**2 + dw_s(2,nj,ls,ms2)*rnorms(1,ls,ms2)*rnorms(2,ls,ms2))
                Sij(5,ni,nj) = Sij(5,ni,nj) + & 
                     y*(-dw_s(2,nj,ls,ms2)*rnorms(1,ls,ms2)**2 + dw_s(1,nj,ls,ms2)*rnorms(1,ls,ms2)*rnorms(2,ls,ms2)) 
             ENDDO
          ENDDO
       ENDDO

       !TEST
       !Sij = 0.d0
       !Sij = Sij /hm1
       !TEST

       DO ki= 1, 3 
          DO ni = 1, n_ws1 
             i = interface%jjs1(ni,ms)
             ib = i + (ki-1)*H_mesh%np 
             DO nj = 1, n_w2
                j = phi_mesh%jj(nj,m2)
                jb = j + 3*H_mesh%np
                DO p = ia(ib),  ia(ib+1) - 1
                   IF (ja(p) == jb) THEN
                      IF (ki == 1)  THEN
                         aa1(p) = aa1(p) + Sij(1,ni,nj)
                         aa2(p) = aa2(p) + Sij(1,ni,nj)
                         EXIT
                      ELSEIF  (ki == 3)  THEN
                         aa1(p) = aa1(p) + Sij(5,ni,nj)
                         aa2(p) = aa2(p) + Sij(5,ni,nj)
                         EXIT
                      ENDIF
                   ENDIF
                END DO
             END DO
          END DO
       END DO
       !TEST SYM
       !Feb 2 2003
       !Sij = 0.d0
       DO ni = 1, n_w2
          i = phi_mesh%jj(ni,m2)
          ib = i + 3*H_mesh%np
          DO kj=1,3
             DO nj = 1, n_ws1
                j = interface%jjs1(nj,ms)
                jb = j + (kj-1)*H_mesh%np
                DO p = ia(ib),  ia(ib+1) - 1
                   IF (ja(p) == jb) THEN
                      IF (kj == 1)  THEN
                         aa1(p) = aa1(p) + Sij(1,nj,ni)
                         aa2(p) = aa2(p) + Sij(1,nj,ni)
                         EXIT
                      ELSEIF  (kj == 3)  THEN
                         aa1(p) = aa1(p) + Sij(5,nj,ni)
                         aa2(p) = aa2(p) + Sij(5,nj,ni)
                         EXIT
                      ENDIF
                   ENDIF
                ENDDO
             END DO
          END DO
       ENDDO
       !TEST SYM
       !Feb 2 2003

       !====================================================================================
       !----------------------(1/sigma) (Rot bj).(Grad(phi_i) x ni)-------------------------
       !====================================================================================
       !        GOTO 200


       Sij = 0.d0

       DO ls = 1, l_Gs

          !--------On calcule le rayon du point gauss
          ray = SUM(phi_mesh%rr(1,phi_mesh%jjs(:,ms2))* phi_mesh%gauss%wws(:,ls))
          x = rjs(ls,ms2)*ray/sigma(m1)

          !terme sans derivee
          DO ni = 1, n_ws2
             DO nj = 1, n_ws1
                y = x * wws(ni,ls)*w_cs(nj,ls)
                Sij(1,ni,nj) = Sij(1,ni,nj) + y*( mode/ray)**2*rnorms(1,ls,ms2)
                Sij(3,ni,nj) = Sij(3,ni,nj) + y*( mode/ray**2)*rnorms(1,ls,ms2)
                Sij(4,ni,nj) = Sij(4,ni,nj) + y*(-mode/ray**2)*rnorms(1,ls,ms2)
                Sij(5,ni,nj) = Sij(5,ni,nj) + y*( mode/ray)**2*rnorms(2,ls,ms2)
             ENDDO
          ENDDO

       ENDDO

       !TEST
       !Sij = 0.d0
       !TEST

       DO ni = 1, n_ws2 
          i = interface%jjs2(ni,ms)
          ib = i + 3*H_mesh%np
          DO kj =1,3
             DO nj = 1, n_ws1
                j = interface%jjs1(nj,ms)
                jb = j + (kj-1)*H_mesh%np
                DO p = ia(ib),  ia(ib+1) - 1
                   IF (ja(p) == jb) THEN
                      IF (kj == 1)  THEN
                         aa1(p) = aa1(p) + Sij(1,ni,nj)
                         aa2(p) = aa2(p) + Sij(1,ni,nj)
                         EXIT                          
                      ELSEIF (kj == 2)  THEN
                         aa1(p) = aa1(p) + Sij(3,ni,nj)
                         aa2(p) = aa2(p) + Sij(4,ni,nj)  
                         EXIT
                      ELSEIF  (kj == 3)  THEN
                         aa1(p) = aa1(p) + Sij(5,ni,nj)
                         aa2(p) = aa2(p) + Sij(5,ni,nj)
                         EXIT
                      ENDIF
                   ENDIF
                END DO
             END DO
          END DO
       END DO

       !Feb 2 2007
       Sij = c_sym*Sij !SYM
       DO ki =1,3
          DO ni = 1, n_ws1
             i = interface%jjs1(ni,ms)
             ib = i + (ki-1)*H_mesh%np
             DO nj = 1, n_ws2
                j = interface%jjs2(nj,ms)
                jb = j + 3*H_mesh%np
                DO p = ia(ib),  ia(ib+1) - 1
                   IF (ja(p) == jb) THEN
                      IF (ki == 1)  THEN
                         aa1(p) = aa1(p) + Sij(1,nj,ni)
                         aa2(p) = aa2(p) + Sij(1,nj,ni)
                         EXIT
                      ELSEIF (ki == 2)  THEN
                         aa1(p) = aa1(p) + Sij(3,nj,ni)
                         aa2(p) = aa2(p) + Sij(4,nj,ni)
                         EXIT
                      ELSEIF  (ki == 3)  THEN
                         aa1(p) = aa1(p) + Sij(5,nj,ni)
                         aa2(p) = aa2(p) + Sij(5,nj,ni)
                         EXIT
                      ENDIF
                   ENDIF
                END DO
             END DO
          END DO
       END DO
       !Feb 2 2007

       Sij = 0.d0

       DO ls = 1, l_Gs

          !--------On calcule le rayon du point gauss
          ray = SUM(phi_mesh%rr(1,phi_mesh%jjs(:,ms2))* phi_mesh%gauss%wws(:,ls))
          x = rjs(ls,ms2)*ray/sigma(m1)

          !terme avec derivee de bi seulement
          DO ni = 1, n_ws2
             y =  x*wws(ni,ls)*mode/ray
             DO nj = 1, n_w1 
                Sij(3,ni,nj) = Sij(3,ni,nj) + &
                     y*(dw_cs(2,nj,ls,ms1)*rnorms(2,ls,ms2) + dw_cs(1,nj,ls,ms1)*rnorms(1,ls,ms2))
             ENDDO
          ENDDO
       ENDDO
       Sij(4,:,:) = -Sij(3,:,:)
       !TEST
       !Sij = 0.d0
       !TEST

       kj=2
       DO ni = 1, n_ws2 
          i = interface%jjs2(ni,ms)
          ib = i + 3*H_mesh%np 
          DO nj = 1, n_w1
             j = H_mesh%jj(nj,m1)
             jb = j + (kj-1)*H_mesh%np
             DO p = ia(ib),  ia(ib+1) - 1
                IF (ja(p) == jb) THEN
                   aa1(p) = aa1(p) + Sij(3,ni,nj)
                   aa2(p) = aa2(p) + Sij(4,ni,nj)  
                   EXIT   
                ENDIF
             END DO
          END DO
       END DO
       !Feb 2 2007
       Sij = c_sym*Sij !SYM
       ki=2
       DO ni = 1, n_w1
          i = H_mesh%jj(ni,m1)
          ib = i + (ki-1)*H_mesh%np
          DO nj = 1, n_ws2
             j = interface%jjs2(nj,ms)
             jb = j + 3*H_mesh%np
             DO p = ia(ib),  ia(ib+1) - 1
                IF (ja(p) == jb) THEN
                   aa1(p) = aa1(p) + Sij(3,nj,ni)
                   aa2(p) = aa2(p) + Sij(4,nj,ni)
                   EXIT
                ENDIF
             END DO
          END DO
       END DO
       !Feb 2 2007

       Sij = 0.d0


       DO ls = 1, l_Gs

          !--------On calcule le rayon du point gauss
          ray = SUM(phi_mesh%rr(1,phi_mesh%jjs(:,ms2))* phi_mesh%gauss%wws(:,ls))
          x = rjs(ls,ms2)*ray/sigma(m1)

          !terme avec derivee de phi et derivee de bi
          DO ni = 1, n_w2
             y =  x*(dw_s(2,ni,ls,ms2)*rnorms(1,ls,ms2) - dw_s(1,ni,ls,ms2)*rnorms(2,ls,ms2))
             DO nj = 1, n_w1
                Sij(1,ni,nj) = Sij(1,ni,nj) +   y *dw_cs(2,nj,ls,ms1) 
                Sij(5,ni,nj) = Sij(5,ni,nj) + (-y)*dw_cs(1,nj,ls,ms1)
             ENDDO
          ENDDO

       ENDDO

       !TEST
       !Sij = 0.d0
       !TEST

       DO ni = 1, n_w2 
          i = phi_mesh%jj(ni,m2)
          ib = i + 3*H_mesh%np
          DO kj=1,3
             DO nj = 1, n_w1
                j = H_mesh%jj(nj,m1)
                jb = j + (kj-1)*H_mesh%np
                DO p = ia(ib),  ia(ib+1) - 1
                   IF (ja(p) == jb) THEN             
                      IF (kj == 1)  THEN
                         aa1(p) = aa1(p) + Sij(1,ni,nj)
                         aa2(p) = aa2(p) + Sij(1,ni,nj)
                         EXIT
                      ELSEIF  (kj == 3)  THEN
                         aa1(p) = aa1(p) + Sij(5,ni,nj)
                         aa2(p) = aa2(p) + Sij(5,ni,nj)
                         EXIT
                      ENDIF
                   ENDIF
                END DO
             END DO
          END DO
       END DO

       !Feb 2 2007
       Sij=c_sym*Sij !SYM
       DO ki=1,3
          DO ni = 1, n_w1
             i = H_mesh%jj(ni,m1)
             ib = i + (ki-1)*H_mesh%np
             DO nj = 1, n_w2
                j = phi_mesh%jj(nj,m2)
                jb = j + 3*H_mesh%np
                DO p = ia(ib),  ia(ib+1) - 1
                   IF (ja(p) == jb) THEN
                      IF (ki == 1)  THEN
                         aa1(p) = aa1(p) + Sij(1,nj,ni)
                         aa2(p) = aa2(p) + Sij(1,nj,ni)
                         EXIT
                      ELSEIF  (ki == 3)  THEN
                         aa1(p) = aa1(p) + Sij(5,nj,ni)
                         aa2(p) = aa2(p) + Sij(5,nj,ni)
                         EXIT
                      ENDIF
                   ENDIF
                END DO
             END DO
          END DO
       END DO
       !Feb 2 2007
       !==================

       !(use .true. for convergence tests)
       !June 6 2008, I put back (.true.) always.
       !Works much better when mu is discontinuous.
       !Mars 22 2007
       !IF (test_de_convergence) THEN 
       IF (.true.) THEN
          !Mars 22 2007 
          !Enforcing weak continuity on the normal components
          Hsij   = 0.d0 
          Sij    = 0.d0
          Phisij = 0.d0

          DO ls = 1, l_Gs

             !Feb 8 2007, muhl
             muhl = SUM(mu_H_field(interface%jjs1(1:n_ws1,ms))*w_cs(1:n_ws1,ls))
             !Feb 8 2007, muhl
             ray = 0.d0
             DO ni = 1, n_ws2;  i = phi_mesh%jjs(ni,ms2)
                ray = ray + phi_mesh%rr(1,i)* phi_mesh%gauss%wws(ni,ls)
             END DO

             !ray = ray*hm1*rjs(ls,ms2)
             !June 8, 2008, Normalization
             ray = ray*hm1*rjs(ls,ms2)/muhl**2
             !June 8, 2008, Normalization
             DO ni = 1, n_ws1
                DO nj = 1, n_ws1
                   x = muhl**2*w_cs(ni,ls)*w_cs(nj,ls)*ray
                   Hsij(1,ni,nj) = Hsij(1,ni,nj) + x*rnorms(1,ls,ms2)**2
                   Hsij(4,ni,nj) = Hsij(4,ni,nj) + x*rnorms(1,ls,ms2)*rnorms(2,ls,ms2)
                   Hsij(6,ni,nj) = Hsij(6,ni,nj) + x*rnorms(2,ls,ms2)**2
                END DO

                DO nj = 1, n_w2
                   x = muhl*mu_phi*w_cs(ni,ls)*(dw_s(1,nj,ls,ms2)*rnorms(1,ls,ms2) +&
                        dw_s(2,nj,ls,ms2)*rnorms(2,ls,ms2))*ray
                   Sij(1,ni,nj) = Sij(1,ni,nj) - x*rnorms(1,ls,ms2)
                   Sij(5,ni,nj) = Sij(5,ni,nj) - x*rnorms(2,ls,ms2)
                ENDDO
             ENDDO

             DO ni = 1, n_w2
                DO nj = 1, n_w2
                   x = mu_phi**2*(dw_s(1,ni,ls,ms2)*rnorms(1,ls,ms2) + dw_s(2,ni,ls,ms2)*rnorms(2,ls,ms2))* &
                        (dw_s(1,nj,ls,ms2)*rnorms(1,ls,ms2) + dw_s(2,nj,ls,ms2)*rnorms(2,ls,ms2))*ray
                   Phisij(ni,nj) = Phisij(ni,nj) + x
                ENDDO
             ENDDO

          END DO
          Sij(2,:,:) = Sij(1,:,:)
          Sij(6,:,:) = Sij(5,:,:)

          DO ni = 1, n_ws1; i = H_mesh%jjs(ni,ms1)
             DO ki= 1, 3, 2
                ib = i + (ki-1)*H_mesh%np 
                DO nj = 1, n_ws1; j = H_mesh%jjs(nj,ms1)
                   DO kj = 1, 3, 2
                      jb = j + (kj-1)*H_mesh%np 
                      DO p = ia(ib), ia(ib+1) - 1
                         IF (ja(p) == jb) THEN
                            IF (ki*kj==1) THEN
                               aa1(p) = aa1(p) + Hsij(1,ni,nj)
                               aa2(p) = aa2(p) + Hsij(1,ni,nj)
                            ELSE IF (ki*kj==9) THEN
                               aa1(p) = aa1(p) + Hsij(6,ni,nj)
                               aa2(p) = aa2(p) + Hsij(6,ni,nj)
                            ELSE IF (ki*kj==3) THEN
                               aa1(p) = aa1(p) + Hsij(4,ni,nj)
                               aa2(p) = aa2(p) + Hsij(4,ni,nj)
                            ELSE
                               WRITE(*,*) ' BUG '; STOP
                            END IF
                            EXIT
                         ENDIF
                      END DO
                   END DO
                END DO

                DO nj = 1, n_w2; j = phi_mesh%jj(nj,m2)
                   jb = j + 3*H_mesh%np 
                   DO p = ia(ib),  ia(ib+1) - 1
                      IF (ja(p) == jb) THEN  
                         aa1(p) = aa1(p) + Sij(2*ki-1,ni,nj)
                         aa2(p) = aa2(p) + Sij(2*ki-1,ni,nj)
                         EXIT  
                      ENDIF
                   END DO
                END DO
             ENDDO
          ENDDO

          DO ni = 1, n_w2; i = phi_mesh%jj(ni,m2)
             ib = i + 3*H_mesh%np
             DO nj = 1, n_ws1; j = H_mesh%jjs(nj,ms1)
                DO kj = 1, 3, 2
                   jb = j + (kj-1)*H_mesh%np 
                   DO p = ia(ib),  ia(ib+1) - 1
                      IF (ja(p) == jb) THEN
                         aa1(p) = aa1(p) + Sij(2*kj-1,nj,ni)
                         aa2(p) = aa2(p) + Sij(2*kj-1,nj,ni)
                         EXIT
                      ENDIF
                   END DO
                END DO
             END DO

             DO nj = 1, n_w2; j = phi_mesh%jj(nj,m2)
                jb = j + 3*H_mesh%np 
                DO p = ia(ib),  ia(ib+1) - 1
                   IF (ja(p) == jb) THEN
                      aa1(p) = aa1(p) + Phisij(ni,nj)  
                      aa2(p) = aa2(p) + Phisij(ni,nj)
                      EXIT
                   ENDIF
                END DO
             END DO
          END DO

       END IF
       !FIN TEST

    ENDDO


    !=========================================================
    !--- Artificial boundary condition: d(phi)/dR + (1/R)*phi = 0
    !=========================================================

    IF (.NOT.PRESENT(index_fourier) .OR. .NOT.PRESENT(R_fourier)) RETURN
    IF (R_fourier.LE.0.d0) RETURN
    WRITE(*,*) ' Assembling the Fourier condition'
    DO ms = 1, phi_mesh%mes
       IF (phi_mesh%sides(ms) /= index_fourier) CYCLE ! Not on the artificial boundary

       Phisij = 0.d0

       DO ls = 1, phi_mesh%gauss%l_Gs

          !--------On calcule le rayon du point gauss
          ray = SUM(phi_mesh%rr(1,phi_mesh%jjs(:,ms))* phi_mesh%gauss%wws(:,ls))

          x = c_mu_phi*rjs(ls,ms)*ray/R_fourier

          DO ni=1,  phi_mesh%gauss%n_ws
             DO nj=1,  phi_mesh%gauss%n_ws
                Phisij(ni,nj) = Phisij(ni,nj) + x*wws(ni,ls)*wws(nj,ls) 
             ENDDO
          ENDDO

       ENDDO

       DO ni = 1, phi_mesh%gauss%n_ws
          i =  phi_mesh%jjs(ni,ms)
          ib = i + 3*H_mesh%np 
          DO nj = 1, phi_mesh%gauss%n_ws
             j = phi_mesh%jjs(nj,ms)
             jb = j + 3*H_mesh%np 
             DO p = ia(ib),  ia(ib+1) - 1
                IF (ja(p) == jb) THEN                   
                   aa1(p) = aa1(p) + Phisij(ni,nj)
                   aa2(p) = aa2(p) + Phisij(ni,nj)
                   EXIT
                ENDIF
             END DO
          END DO
       END DO
    END DO

  END SUBROUTINE mat_maxwell


  SUBROUTINE courant(H_mesh,phi_mesh,INTERFACE,sigma,mu_phi,mu_H_field,time,mode,&
       rhs_H,rhs_phi,nl,src_H,src_phi)
    !forcage faisant intervenir J, volumique et interface pour H et phi
    !pour le probleme en entier

    USE def_type_mesh
    USE gauss_points
    USE boundary

    IMPLICIT NONE
    TYPE(mesh_type),                       INTENT(IN)   :: H_mesh, phi_mesh
    TYPE(interface_type),                  INTENT(IN)   :: INTERFACE
    REAL(KIND=8),                          INTENT(IN)   :: mu_phi, time
    REAL(KIND=8), DIMENSION(H_mesh%me),    INTENT(IN)   :: sigma
    REAL(KIND=8), DIMENSION(H_mesh%np),    INTENT(IN)   :: mu_H_field
    INTEGER,                               INTENT(IN)   :: mode
    REAL(KIND=8), DIMENSION(H_mesh%np,6),  INTENT(IN)   :: nl
    REAL(KIND=8), DIMENSION(H_mesh%np,6),  INTENT(IN)   :: rhs_H
    REAL(KIND=8), DIMENSION(phi_mesh%np,2),INTENT(IN)   :: rhs_phi  
    REAL(KIND=8), DIMENSION(:,:),          INTENT(INOUT):: src_H, src_phi

    REAL(KIND=8), DIMENSION(phi_mesh%gauss%n_ws,phi_mesh%gauss%l_Gs)   :: w_cs   
    REAL(KIND=8), DIMENSION(2) :: gaussp
    REAL(KIND=8) :: j_courant, rot, x, ray
    INTEGER :: m, l, i, i_b, ni, n, k, k1, ms, ls, n_ws1, n_ws2, ms1, ms2, H_bloc_size, n_w2, m1
    REAL(KIND=8), DIMENSION(6)             :: JsolH_anal, rhs_Hl
    REAL(KIND=8), DIMENSION(2)             :: rhs_phil
    REAL(KIND=8), DIMENSION(2,4)           :: rhs_dphil
    REAL(KIND=8), DIMENSION(2,H_mesh%gauss%n_w) :: dwH
    REAL(KIND=8), DIMENSION(2,phi_mesh%gauss%n_w) :: dwphi
    REAL(KIND=8), DIMENSION(2,phi_mesh%gauss%n_w) :: src_phil
    REAL(KIND=8) :: ray_rjl, moderay2, muhl

    !forcage volumique  
    !attention on comprime le calcul sur les points de Gauss et integration !!
    !j/sigma *(Rot(b))

    src_H=0
    src_phi=0 
    DO m = 1, H_mesh%me     
       DO l = 1, H_mesh%gauss%l_G
          !Feb 8 2007, muhl
          muhl=SUM(mu_H_field(H_mesh%jj(:,m))*H_mesh%gauss%ww(:,l))
          !Feb 8 2007, muhl
          dwH = H_mesh%gauss%dw(:,:,l,m)
          !--------On calcule le rayon du point gauss

          JsolH_anal = 0.d0
          rhs_Hl = 0.d0
          gaussp = 0.d0 
          DO ni = 1, H_mesh%gauss%n_w;  i = H_mesh%jj(ni,m)
             gaussp = gaussp + H_mesh%rr(:,i)*H_mesh%gauss%ww(ni,l)
             JsolH_anal(:) = JsolH_anal(:) + NL(i,:)*H_mesh%gauss%ww(ni,l)
             rhs_Hl(:) = rhs_Hl(:) + rhs_H(i,:)*H_mesh%gauss%ww(ni,l)
          ENDDO
          ray = gaussp(1)
          ray_rjl = H_mesh%gauss%rj(l,m)*ray

          DO k=1, 6
             JsolH_anal(k) = muhl*JsolH_anal(k) + &
                  Jexact_gauss(k, gaussp, mode, mu_phi, sigma(m), muhl, time)/sigma(m)
          END DO

          DO ni = 1,H_mesh%gauss%n_w

             i = H_mesh%jj(ni,m)

             !--------Composante r------
             src_H(i,1) = src_H(i,1)+ ray_rjl &
                  *(JsolH_anal(3)*dwH(2,ni) &
                  + mode/ray*JsolH_anal(6)*H_mesh%gauss%ww(ni,l) &
                  + rhs_Hl(1)*H_mesh%gauss%ww(ni,l))

             src_H(i,2) = src_H(i,2)+ ray_rjl  &
                  *(JsolH_anal(4)*dwH(2,ni) &
                  - mode/ray*JsolH_anal(5)*H_mesh%gauss%ww(ni,l) &
                  + rhs_Hl(2)*H_mesh%gauss%ww(ni,l))   

             !--------Composante theta------
             src_H(i,3) = src_H(i,3)+ ray_rjl  &
                  * (-JsolH_anal(1)*dwH(2,ni)  &
                  + 1/ray*JsolH_anal(5)*(ray*dwH(1,ni) + H_mesh%gauss%ww(ni,l)) &
                  + rhs_Hl(3)*H_mesh%gauss%ww(ni,l)) 

             src_H(i,4) = src_H(i,4)+ ray_rjl &
                  * (-JsolH_anal(2)*dwH(2,ni) &
                  + 1/ray*JsolH_anal(6)*(ray*dwH(1,ni) + H_mesh%gauss%ww(ni,l)) &
                  + rhs_Hl(4)*H_mesh%gauss%ww(ni,l))

             !--------Composante z------
             src_H(i,5) = src_H(i,5)+ ray_rjl* &
                  (-mode/ray*JsolH_anal(2)*H_mesh%gauss%ww(ni,l) &
                  - JsolH_anal(3)*dwH(1,ni) &
                  + rhs_Hl(5)*H_mesh%gauss%ww(ni,l))

             src_H(i,6) = src_H(i,6)+ ray_rjl* &
                  (mode/ray*JsolH_anal(1)*H_mesh%gauss%ww(ni,l) &
                  - JsolH_anal(4)*dwH(1,ni) &
                  + rhs_Hl(6)*H_mesh%gauss%ww(ni,l)) 
          ENDDO

       END DO
    END DO

    DO m = 1, phi_mesh%me
       src_phil=0
       DO l = 1, phi_mesh%gauss%l_G
          dwphi = phi_mesh%gauss%dw(:,:,l,m)
          !--------On calcule le rayon du point gauss
          rhs_dphil=0
          rhs_phil=0
          ray = 0
          DO ni = 1, phi_mesh%gauss%n_w;  i = phi_mesh%jj(ni,m)
             ray = ray + phi_mesh%rr(1,i)*phi_mesh%gauss%ww(ni,l)
             rhs_phil(:) = rhs_phil(:) + rhs_phi(i,:)*phi_mesh%gauss%ww(ni,l)
             DO k =1 ,2
                rhs_dphil(:,k) = rhs_dphil(:,k) + rhs_phi(i,:)*dwphi(k,ni)
             END DO
          END DO
          ray_rjl = phi_mesh%gauss%rj(l,m)*ray
          moderay2 = (mode/ray)**2

          DO ni = 1, phi_mesh%gauss%n_w

             src_phil(1,ni) =  src_phil(1,ni) + ray_rjl* &
                  (rhs_dphil(1,1)*dwphi(1,ni) + &
                  moderay2*rhs_phil(1)*phi_mesh%gauss%ww(ni,l) + &
                  rhs_dphil(1,2)*dwphi(2,ni))

             src_phil(2,ni) =  src_phil(2,ni) + ray_rjl* &
                  (rhs_dphil(2,1)*dwphi(1,ni) + &
                  moderay2*rhs_phil(2)*phi_mesh%gauss%ww(ni,l) + &
                  rhs_dphil(2,2)*dwphi(2,ni))
          END DO

       END DO
       DO ni = 1, phi_mesh%gauss%n_w
          i = phi_mesh%jj(ni,m)
          src_phi(i,:) = src_phi(i,:) + src_phil(:,ni) 
       END DO
    END DO

    !==Interface
    !forcage sur l'interface
    !attention on comprime le calcul sur les points de Gauss et integration !!
    !j/sigma*(b x nc + grad(phi) x nv)

    CALL gauss(phi_mesh)

    n_ws1 = H_mesh%gauss%n_ws
    n_ws2 = phi_mesh%gauss%n_ws
    n_w2 = phi_mesh%gauss%n_w

    H_bloc_size = H_mesh%np

    IF (H_mesh%gauss%n_ws == n_ws) THEN
       w_cs = wws
    ELSE    
       DO ls = 1, l_Gs
          w_cs(1,ls)= wws(1,ls)+0.5*wws(3,ls)
          w_cs(2,ls)= wws(2,ls)+0.5*wws(3,ls)
          w_cs(3,ls)=0
       ENDDO
    END IF

    DO ms = 1, interface%mes

       ms2 = interface%mesh2(ms)
       ms1 = interface%mesh1(ms)
       m = phi_mesh%neighs(ms2)           
       m1 = H_mesh%neighs(ms1)

       DO ls = 1,l_Gs
          !Feb 9 2007, muhl
          muhl=SUM(mu_H_field(interface%jjs1(1:n_ws1,ms))*w_cs(1:n_ws1,ls))
          !Feb 9 2007, muhl

          !--------On calcule le rayon du point gauss
          ray = 0
          DO ni = 1, n_ws2;  i = phi_mesh%jjs(ni,interface%mesh2(ms))
             ray = ray + phi_mesh%rr(1,i)* wws(ni,ls)
          END DO

          gaussp = 0.d0    
          DO ni=1, n_ws2
             i=phi_mesh%jjs(ni,ms2)
             gaussp = gaussp + phi_mesh%rr(:,i)*phi_mesh%gauss%wws(ni,ls)
          ENDDO

          DO k=1, 6
             JsolH_anal(k) = Jexact_gauss(k, gaussp, mode, mu_phi ,sigma(m1), muhl, time)/sigma(m1) &
                  + muhl * SUM(NL(H_mesh%jjs(1:n_ws1,ms1),k)*w_cs(1:n_ws1,ls))
          ENDDO

          !---------forcage interface pour H            

          DO ni=1, n_ws1
             i = interface%jjs1(ni,ms)
             src_H(i,1) = src_H(i,1)+rjs(ls,ms2)*ray*( &
                  -JsolH_anal(3)*w_cs(ni,ls)*(-rnorms(2,ls,ms2)))

             src_H(i,2) = src_H(i,2)+rjs(ls,ms2)*ray*( &
                  -JsolH_anal(4)*w_cs(ni,ls)*(-rnorms(2,ls,ms2)))

             src_H(i,3) = src_H(i,3)+rjs(ls,ms2)*ray*( &
                  JsolH_anal(1)*w_cs(ni,ls)*(-rnorms(2,ls,ms2)) &
                  -JsolH_anal(5)*w_cs(ni,ls)*(-rnorms(1,ls,ms2)))

             src_H(i,4) = src_H(i,4)+rjs(ls,ms2)*ray*( &
                  JsolH_anal(2)*w_cs(ni,ls)*(-rnorms(2,ls,ms2))  &
                  -JsolH_anal(6)*w_cs(ni,ls)*(-rnorms(1,ls,ms2)))

             src_H(i,5) = src_H(i,5)+rjs(ls,ms2)*ray*( &
                  JsolH_anal(3)*w_cs(ni,ls)*(-rnorms(1,ls,ms2)))

             src_H(i,6) = src_H(i,6)+rjs(ls,ms2)*ray*( &
                  JsolH_anal(4)*w_cs(ni,ls)*(-rnorms(1,ls,ms2)))
          ENDDO

          !---------forcage interface pour phi            
          !terme sans derivee de phi
          DO ni=1,n_ws2           
             i = interface%jjs2(ni,ms)
             !attention si on force sur l'axe, il faut retirer les 1/ray
             !There was a BUG here. There was w_cs instead of wws
             src_phi(i,1) = src_phi(i,1)+rjs(ls,ms2)*( &
                  - mode*JsolH_anal(2)*wws(ni,ls) * rnorms(2,ls,ms2) &
                  + mode*JsolH_anal(6)*wws(ni,ls) * rnorms(1,ls,ms2))

             src_phi(i,2) = src_phi(i,2)+rjs(ls,ms2)*( &
                  + mode*JsolH_anal(1)*wws(ni,ls) * rnorms(2,ls,ms2) &
                  - mode*JsolH_anal(5)*wws(ni,ls) * rnorms(1,ls,ms2))

          ENDDO

          !terme avec derivee de phi
          DO ni=1,n_w2           
             i = phi_mesh%jj(ni,m)
             src_phi(i,1) = src_phi(i,1)+rjs(ls,ms2)*ray*( &
                  + JsolH_anal(3) *(dw_s(2,ni,ls,ms2) * rnorms(1,ls,ms2)&
                  -dw_s(1,ni,ls,ms2) * rnorms(2,ls,ms2)))  

             src_phi(i,2) = src_phi(i,2)+rjs(ls,ms2)*ray*( &
                  + JsolH_anal(4)*(dw_s(2,ni,ls,ms2) * rnorms(1,ls,ms2)&
                  -dw_s(1,ni,ls,ms2) * rnorms(2,ls,ms2)))

          ENDDO

       END DO
    END DO

  END SUBROUTINE courant



  SUBROUTINE surf_int(H_mesh,phi_mesh,INTERFACE,interface_H_mu,sigma,mu_phi,mu_H_field,time,mode,src_H, src_phi, &
       R_fourier, index_fourier, rhs)
    !calcul du forcage a la frontiere exterieure

    USE def_type_mesh
    USE boundary

    IMPLICIT NONE
    TYPE(mesh_type),              INTENT(IN)    :: H_mesh, phi_mesh
    TYPE(interface_type),         INTENT(IN)    :: INTERFACE, interface_H_mu
    REAL(KIND=8),                 INTENT(IN)    :: mu_phi,  time
    REAL(KIND=8),DIMENSION(H_mesh%me),INTENT(IN):: sigma
    REAL(KIND=8),DIMENSION(H_mesh%np),INTENT(IN):: mu_H_field 
    INTEGER,                      INTENT(IN)    :: mode   
    REAL(KIND=8), DIMENSION(:,:), INTENT(INOUT) :: src_H, src_phi
    REAL(KIND=8), DIMENSION(phi_mesh%np,2), OPTIONAL :: rhs
    REAL(KIND=8),               OPTIONAL :: R_fourier
    INTEGER,                    OPTIONAL :: index_fourier

    REAL(KIND=8), DIMENSION(6) :: gls, ngls
    REAL(KIND=8)               :: gdxn, ray, x, y1, y2, muhl
    INTEGER                    :: ms, ls, ns, i_b, i, k, k1, m, n, ni
    INTEGER                    ::  ms22, ms21
    LOGICAL,      DIMENSION(H_mesh%mes)    :: virgin1 !=.true.
    LOGICAL,      DIMENSION(phi_mesh%mes)  :: virgin2 !=.true.
    REAL(KIND=8), DIMENSION(2)             :: gaussp
    REAL(KIND=8), DIMENSION(6)             :: EsolH_anal, Esolphi_anal
    REAL(KIND=8), DIMENSION(6)             :: ff 

    virgin1=.TRUE.
    virgin2=.TRUE.
    virgin1(interface%mesh1) = .FALSE.
    virgin2(interface%mesh2) = .FALSE.
    virgin1(interface_H_mu%mesh1) = .FALSE.
    virgin1(interface_H_mu%mesh2) = .FALSE.


    DO ms = 1,H_mesh%mes
       IF (.NOT.virgin1(ms)) CYCLE

       m = H_mesh%neighs(ms)
       DO ls = 1, H_mesh%gauss%l_Gs
          !Feb 8 2007, mmuhl
          muhl = SUM(mu_H_field(H_mesh%jjs(:,ms))*H_mesh%gauss%wws(:,ls))
          !Feb 8 2007, mmuhl

          !--------On calcule le rayon du point gauss
          ray = 0
          DO ni = 1, H_mesh%gauss%n_ws;  i = H_mesh%jjs(ni,ms)
             ray = ray + H_mesh%rr(1,i)* H_mesh%gauss%wws(ni,ls)
          END DO

          IF (ray.LT.1.d-15) CYCLE !ATTENTION Axe

          gaussp = 0.d0    
          DO ns=1, H_mesh%gauss%n_ws
             i=H_mesh%jjs(ns,ms)
             gaussp = gaussp + H_mesh%rr(:,i)*H_mesh%gauss%wws(ns,ls)
          ENDDO

          DO k=1, 6
             EsolH_anal(k) = Eexact_gauss(k,gaussp,mode,mu_phi,sigma(m),muhl, time)
          ENDDO

          !forcage pour la frontiere de H 
          ! - E.(b x nc)

          DO ns=1, H_mesh%gauss%n_ws
             i = H_mesh%jjs(ns,ms)
             src_H(i,1) = src_H(i,1)-H_mesh%gauss%rjs(ls,ms)*ray*( &
                  -EsolH_anal(3)*H_mesh%gauss%wws(ns,ls)* &
                  (H_mesh%gauss%rnorms(2,ls,ms)))

             src_H(i,2) = src_H(i,2)-H_mesh%gauss%rjs(ls,ms)*ray*( &
                  -EsolH_anal(4)*H_mesh%gauss%wws(ns,ls)* &
                  (H_mesh%gauss%rnorms(2,ls,ms)))

             src_H(i,3) = src_H(i,3)-H_mesh%gauss%rjs(ls,ms)*ray*( &
                  EsolH_anal(1)*H_mesh%gauss%wws(ns,ls)* &
                  (H_mesh%gauss%rnorms(2,ls,ms)) - &
                  EsolH_anal(5)*H_mesh%gauss%wws(ns,ls) * &
                  (H_mesh%gauss%rnorms(1,ls,ms)))

             src_H(i,4) = src_H(i,4)-H_mesh%gauss%rjs(ls,ms)*ray*( &
                  EsolH_anal(2)*H_mesh%gauss%wws(ns,ls) * &
                  (H_mesh%gauss%rnorms(2,ls,ms)) - &
                  EsolH_anal(6)*H_mesh%gauss%wws(ns,ls) * &
                  (H_mesh%gauss%rnorms(1,ls,ms)))

             src_H(i,5) = src_H(i,5)-H_mesh%gauss%rjs(ls,ms)*ray*( &
                  EsolH_anal(3)*H_mesh%gauss%wws(ns,ls)*  &
                  (H_mesh%gauss%rnorms(1,ls,ms)))

             src_H(i,6) = src_H(i,6)-H_mesh%gauss%rjs(ls,ms)*ray*( &
                  EsolH_anal(4)*H_mesh%gauss%wws(ns,ls) * &
                  (H_mesh%gauss%rnorms(1,ls,ms)))

          ENDDO

       ENDDO
    ENDDO

    !==boundary phi_mesh
    DO ms = 1,phi_mesh%mes

       IF (PRESENT(index_fourier)) THEN
          IF (phi_mesh%sides(ms)==index_fourier) CYCLE
       END IF

       IF (.NOT.virgin2(ms)) CYCLE
       m = phi_mesh%neighs(ms)

       DO ls = 1, phi_mesh%gauss%l_Gs

          !--------On calcule le rayon du point gauss
          ray = 0
          DO ni = 1, phi_mesh%gauss%n_ws;  i = phi_mesh%jjs(ni,ms)
             ray = ray + phi_mesh%rr(1,i)* phi_mesh%gauss%wws(ni,ls)
          END DO

          IF (ray.LT.1.d-10) CYCLE !ATTENTION Axe

          gaussp = 0.d0    
          DO ns=1, phi_mesh%gauss%n_ws
             i=phi_mesh%jjs(ns,ms)
             gaussp = gaussp + phi_mesh%rr(:,i)*phi_mesh%gauss%wws(ns,ls)
          ENDDO

          DO k=1, 6
             ! Here sigma and mu_H_field should not intervene
             ! I put boggus values, Feb 8 2007, Jean-Luc Guermond       
             Esolphi_anal(k) = Eexact_gauss(k,gaussp,mode,mu_phi,sigma(1),mu_H_field(1),time)
          ENDDO

          !forcage pour la frontiere de phi 
          ! - E.(grad(phi) x nv)

          DO ns=1, phi_mesh%gauss%n_ws           
             i = phi_mesh%jjs(ns,ms) 
             DO n = 1, phi_mesh%gauss%n_w
                IF (phi_mesh%jj(n,m) == i) EXIT
             END DO

             !attention si on force sur l'axe on retire les 1/ray 
             src_phi(i,1) = src_phi(i,1)-phi_mesh%gauss%rjs(ls,ms)*ray*( &
                  +Esolphi_anal(3)*(phi_mesh%gauss%dw_s(2,n,ls,ms)*phi_mesh%gauss%rnorms(1,ls,ms) &
                  -phi_mesh%gauss%dw_s(1,n,ls,ms)*phi_mesh%gauss%rnorms(2,ls,ms))) &
                  -phi_mesh%gauss%rjs(ls,ms)*(-mode*Esolphi_anal(2)*phi_mesh%gauss%wws(ns,ls)*phi_mesh%gauss%rnorms(2,ls,ms) &
                  +mode*Esolphi_anal(6)*phi_mesh%gauss%wws(ns,ls)*phi_mesh%gauss%rnorms(1,ls,ms))

             src_phi(i,2) = src_phi(i,2)-phi_mesh%gauss%rjs(ls,ms)*ray*( &
                  +Esolphi_anal(4)*(phi_mesh%gauss%dw_s(2,n,ls,ms)*phi_mesh%gauss%rnorms(1,ls,ms) &
                  -phi_mesh%gauss%dw_s(1,n,ls,ms)*phi_mesh%gauss%rnorms(2,ls,ms))) & 
                  -phi_mesh%gauss%rjs(ls,ms)*( mode*Esolphi_anal(1)*phi_mesh%gauss%wws(ns,ls)*phi_mesh%gauss%rnorms(2,ls,ms) &
                  -mode*Esolphi_anal(5)*phi_mesh%gauss%wws(ns,ls)*phi_mesh%gauss%rnorms(1,ls,ms))


          ENDDO

       ENDDO
    ENDDO

    !=========================================================
    !--- Artificial boundary condition: d(phi_t)/dR + (1/R)*phi_t = 0
    !=========================================================

    !IF (.NOT.present(index_fourier) .OR. .NOT.present(R_fourier)) RETURN
    !IF (R_fourier.le.0.d0) RETURN
    !DO ms = 1, phi_mesh%mes
    !   IF (phi_mesh%sides(ms) /= index_fourier) CYCLE ! Not on the artificial boundary

    !   DO ls = 1, phi_mesh%gauss%l_Gs

    !--------On calcule le rayon du point gauss
    !      ray = SUM(phi_mesh%rr(1,phi_mesh%jjs(:,ms))* phi_mesh%gauss%wws(:,ls))
    !      x = phi_mesh%gauss%rjs(ls,ms)*ray/R_fourier
    !      y1 = x* SUM(rhs(phi_mesh%jjs(:,ms),1)* phi_mesh%gauss%wws(:,ls))
    !      y2 = x* SUM(rhs(phi_mesh%jjs(:,ms),2)* phi_mesh%gauss%wws(:,ls))
    !      DO ns =1, phi_mesh%gauss%n_ws    
    !         src_phi(1,phi_mesh%jjs(ns,ms)) = src_phi(1,phi_mesh%jjs(ns,ms)) + &
    !               y1*phi_mesh%gauss%wws(ns,ls)
    !         src_phi(2,phi_mesh%jjs(ns,ms)) = src_phi(2,phi_mesh%jjs(ns,ms)) + &
    !               y2*phi_mesh%gauss%wws(ns,ls)
    !      ENDDO
    !      
    !   ENDDO
    !END DO

  END SUBROUTINE surf_int

  SUBROUTINE  mat_maxwell_mu(H_mesh, interface_H_mu, mode, stab, ia, ja, aa1, aa2, mu_H_field, sigma)
    USE def_type_mesh
    USE gauss_points

    IMPLICIT NONE
    TYPE(mesh_type),            INTENT(IN)    :: H_mesh
    TYPE(interface_type),       INTENT(IN)    :: interface_H_mu
    INTEGER,                    INTENT(IN)    :: mode
    REAL(KIND=8), DIMENSION(3), INTENT(IN)    :: stab
    INTEGER,      DIMENSION(:), INTENT(IN)    :: ia, ja
    REAL(KIND=8), DIMENSION(:), INTENT(INOUT) :: aa1, aa2
    REAL(KIND=8), DIMENSION(:), INTENT(IN)    :: sigma, mu_H_field

    INTEGER :: m, l, ms, ls, ni, nj, k, h, k1, h1, i, j, i_b, j_b, H_bloc_size, p, &
         n_ws1, n_ws2, nj1, nj2, n_w2, n_w1, m1, m2, ki, kj,ib,jb, ms1, ms2
    REAL(KIND=8) :: x, y, z, norm, hm1
    REAL(KIND=8) :: b, b1, b2, b3, b4, b5, b6, ray, stab_colle_H_mu
    LOGICAL :: mark=.FALSE.

    REAL(KIND=8), DIMENSION(9,H_mesh%gauss%n_w,H_mesh%gauss%n_w,2,2)      :: Hsij, Gsij

    ! MATRICES POUR LES TERMES DE BORDS Hsij et Gsij 
    !=================================================
    ! (--------------------------------------------------)
    ! ( Hsij(1) +G     | GSij(2)        | Hsij(4) +G     )
    ! ( Hsij(1) +G     | GSij(3)        | Hsij(4) +G     )
    ! (--------------------------------------------------)
    ! ( Hsij(2)        | Hsij(5) +G     | Hsij(8)        )
    ! ( Hsij(3)        | Hsij(5) +G     | Hsij(9)        )
    ! (--------------------------------------------------)
    ! ( Hsij(7) +G     |  GSij(8)       | Hsij(6) +G     )
    ! ( Hsij(7) +G     |  GSij(9)       | Hsij(6) +G     )
    ! (==================================================)


    REAL(KIND=8), DIMENSION(H_mesh%gauss%n_ws,H_mesh%gauss%l_Gs)   :: w_cs
    REAL(KIND=8), DIMENSION(2,  H_mesh%gauss%n_w, H_mesh%gauss%l_Gs, H_mesh%mes) :: dw_cs
    REAL(KIND=8), DIMENSION(2,  H_mesh%gauss%n_w) :: dwsi,dwsj
    REAL(KIND=8), DIMENSION(2,H_mesh%gauss%l_Gs)  :: gauss1, gauss2
    REAL(KIND=8), DIMENSION(2)                    :: normi, normj
    INTEGER,      DIMENSION(2)                    :: jjsi, jjsj
    REAL(KIND=8), DIMENSION(H_mesh%gauss%n_ws)    :: wwsi, wwsj 
    INTEGER                                       :: n_wsi, n_wsj, ci, cj, n_wi, n_wj

    INTEGER      :: ls1, ls2, n, n1, n2, ni1, ni2, nip
    REAL(KIND=8) :: ref, diff, d1, d2, mu_H, c_mu_H, h2, muhl1, muhl2, muhi, muhj, sigmai, sigmaj
    ! June 14 2008
    REAL(KIND=8) :: c_sym = 0.d0 ! (c_sym=1.d0 symmetrizes the bilinear form)
    ! June 14 2008


    stab_colle_H_mu = stab(3)

    !**********************************************************************************
    !--------------------TERMES SUR L'INTERFACE SIGMA_MU-------------------------------
    !**********************************************************************************

    WRITE(*,*) 'Assemblage interface'
    CALL gauss(H_mesh)
    n_ws1 = H_mesh%gauss%n_ws
    n_ws2 = H_mesh%gauss%n_ws
    n_w1  = H_mesh%gauss%n_w
    n_w2  = H_mesh%gauss%n_w

    DO ms = 1, interface_H_mu%mes

       ms2 = interface_H_mu%mesh2(ms)
       m2 = H_mesh%neighs(ms2)
       ms1 = interface_H_mu%mesh1(ms)
       m1 = H_mesh%neighs(ms1)

       ref = 1.d-8+SUM((H_mesh%rr(:,H_mesh%jjs(1,ms1)) - H_mesh%rr(:,H_mesh%jjs(2,ms1)))**2)
       diff =      SUM((H_mesh%rr(:,H_mesh%jjs(1,ms1)) - H_mesh%rr(:,H_mesh%jjs(1,ms2)))**2)
       IF (diff/ref .LT. 1.d-10) THEN ! 1 = 1
          w_cs = wws
       ELSE                ! 1 = 2
          DO ls = 1, l_Gs
             w_cs(1,ls)= wws(2,ls)
             w_cs(2,ls)= wws(1,ls)
             IF (n_ws1==3) w_cs(n_ws1,ls) = wws(n_ws1,ls) 
             WRITE(*,*) ' Ouaps! oder of shape functions changed?'
          END DO
       END IF

       DO ls = 1, l_Gs
          gauss2(1,ls) = SUM(H_mesh%rr(1,H_mesh%jjs(:,ms2))*H_mesh%gauss%wws(:,ls))
          gauss2(2,ls) = SUM(H_mesh%rr(2,H_mesh%jjs(:,ms2))*H_mesh%gauss%wws(:,ls))
          gauss1(1,ls) = SUM(H_mesh%rr(1,H_mesh%jjs(:,ms1))*H_mesh%gauss%wws(:,ls))
          gauss1(2,ls) = SUM(H_mesh%rr(2,H_mesh%jjs(:,ms1))*H_mesh%gauss%wws(:,ls))
       END DO

       DO ls2 = 1, l_Gs
          ref = SQRT(1.d-8+SUM(gauss2(:,ls2)**2))
          mark = .FALSE.
          DO ls1 = 1, l_Gs
             diff = SQRT(SUM((gauss2(:,ls2)-gauss1(:,ls1))**2))
             IF (diff .LT. 1.d-10) THEN
                dw_cs(:,:,ls2,ms1) =  H_mesh%gauss%dw_s(:,:,ls1,ms1)
                mark = .TRUE.
                EXIT
             END IF
          END DO
          IF (.NOT.mark) WRITE(*,*) ' BUG '
       END DO

    END DO

    DO ms = 1, interface_H_mu%mes

       ms2 = interface_H_mu%mesh2(ms)
       ms1 = interface_H_mu%mesh1(ms)
       m2 = H_mesh%neighs(ms2)
       m1 = H_mesh%neighs(ms1)
       mu_H = SUM(mu_H_field(H_mesh%jj(:,m1)))/H_mesh%gauss%n_w
       hm1 = stab_colle_H_mu*(((mu_H+mu_H)/mu_H)/SUM(rjs(:,ms2)))

       !====================================================================================
       !------------------------------------TERMES SUR LE BLOC H----------------------------
       !====================================================================================

       !-------------------------------hm1 (bi x ni) . (bj x nj)----------------------------
       !====================================================================================
       Hsij = 0.d0  
       DO ls = 1, l_Gs
          !--------On calcule le rayon du point gauss
          ray = SUM(H_mesh%rr(1,H_mesh%jjs(:,ms2))* H_mesh%gauss%wws(:,ls))
          x = hm1*rjs(ls,ms2)*ray

          !June 14 2008, muhl
          muhl1 = SUM(mu_H_field(H_mesh%jjs(:,ms1))*w_cs(:,ls))
          muhl2 = SUM(mu_H_field(H_mesh%jjs(:,ms2))* wws(:,ls))
          norm = 1.d0/MAX(muhl1,muhl2)**2
          !June 14 2008, muhl

          DO ci = 1, 2
             IF (ci==1) THEN
                normi = rnorms(:,ls,ms1)
                wwsi = w_cs(:,ls)
                n_wsi = n_ws1
                muhi = muhl1
             ELSE
                normi = rnorms(:,ls,ms2)
                wwsi = wws(:,ls)
                n_wsi = n_ws2
                muhi = muhl2
             END IF
             DO cj = 1, 2 
                IF (cj==1) THEN
                   normj = rnorms(:,ls,ms1)
                   wwsj = w_cs(:,ls)
                   n_wsj = n_ws1
                   muhj = muhl1
                ELSE
                   normj = rnorms(:,ls,ms2)
                   wwsj = wws(:,ls)
                   n_wsj = n_ws2
                   muhj = muhl2
                END IF

                DO ni = 1, n_wsi
                   DO nj = 1, n_wsj
                      y = x * wwsi(ni)*wwsj(nj)
                      ! June 14 2008, added z
                      z = y * norm * muhi * muhj
                      ! June 14 2008, added z
                      Hsij(1,ni,nj,ci,cj) = Hsij(1,ni,nj,ci,cj) + y*normi(2)*normj(2) &
                           + z*normi(1)*normj(1)
                      Hsij(4,ni,nj,ci,cj) = Hsij(4,ni,nj,ci,cj) - y*normj(1)*normi(2) &
                           + z*normi(1)*normj(2)
                      Hsij(5,ni,nj,ci,cj) = Hsij(5,ni,nj,ci,cj) + y*(normi(1)*normj(1) + normi(2)*normj(2)) 
                      Hsij(6,ni,nj,ci,cj) = Hsij(6,ni,nj,ci,cj) + y*normi(1)*normj(1) &
                           + z*normi(2)*normj(2)
                   END DO
                END DO
             END DO
          END DO
       END DO

       DO ci = 1, 2
          DO ki = 1, 3 
             DO ni = 1, n_wsi 
                IF (ci==1) THEN
                   i = interface_H_mu%jjs1(ni,ms)
                ELSE
                   i = interface_H_mu%jjs2(ni,ms)
                END IF
                ib = i + (ki-1)*H_mesh%np 
                DO cj = 1, 2
                   DO kj = 1, 3
                      DO nj = 1, n_wsj
                         IF (cj==1) THEN
                            j = interface_H_mu%jjs1(nj,ms)
                         ELSE
                            j = interface_H_mu%jjs2(nj,ms)
                         END IF
                         jb = j + (kj-1)*H_mesh%np 
                         DO p = ia(ib),  ia(ib+1) - 1
                            IF (ja(p) == jb) THEN
                               IF  ((ki == 1) .AND. (kj == 1)) THEN
                                  aa1(p) = aa1(p) + Hsij(1,ni,nj,ci,cj)
                                  aa2(p) = aa2(p) + Hsij(1,ni,nj,ci,cj)
                                  EXIT
                               ELSEIF  ((ki == 1) .AND. (kj == 3)) THEN
                                  aa1(p) = aa1(p) + Hsij(4,ni,nj,ci,cj)
                                  aa2(p) = aa2(p) + Hsij(4,ni,nj,ci,cj)
                                  EXIT                                                                                
                               ELSEIF  ((ki == 3) .AND. (kj == 1)) THEN
                                  aa1(p) = aa1(p) + Hsij(4,nj,ni,cj,ci)
                                  aa2(p) = aa2(p) + Hsij(4,nj,ni,cj,ci)
                                  EXIT
                               ELSEIF  ((ki == 2) .AND. (kj == 2)) THEN
                                  aa1(p) = aa1(p) + Hsij(5,ni,nj,ci,cj) 
                                  aa2(p) = aa2(p) + Hsij(5,ni,nj,ci,cj)
                                  EXIT                                                                                
                               ELSEIF  ((ki == 3) .AND. (kj == 3)) THEN
                                  aa1(p) = aa1(p) + Hsij(6,ni,nj,ci,cj)
                                  aa2(p) = aa2(p) + Hsij(6,ni,nj,ci,cj)
                                  EXIT 
                               ENDIF
                               !WRITE(*,*) ' BUG in mat_maxwell_mu '
                            ENDIF
                         END DO
                      END DO
                   END DO
                END DO
             END DO
          END DO
       END DO

       !====================================================================================
       !------------------------(1/sigma) (Rot bj) . (bi x ni)------------------------------
       !====================================================================================

       !terme sans derivee
       Hsij = 0.d0
       Gsij = 0.d0
       DO ls = 1, H_mesh%gauss%l_Gs

          !--------On calcule le rayon du point gauss
          ray = SUM(H_mesh%rr(1,H_mesh%jjs(:,ms2))* H_mesh%gauss%wws(:,ls))
          x = rjs(ls,ms2)*ray

          DO ci = 1, 2
             IF (ci==1) THEN
                normi = rnorms(:,ls,ms1)
                wwsi = w_cs(:,ls)
                n_wsi = n_ws1
                sigmai = sigma(m1)
             ELSE
                normi = rnorms(:,ls,ms2)
                wwsi = wws(:,ls)
                n_wsi = n_ws2
                sigmai = sigma(m2)
             END IF
             DO cj = 1, 2
                IF (cj==1) THEN
                   normj = rnorms(:,ls,ms1)
                   wwsj = w_cs(:,ls)
                   n_wsj = n_ws1
                   sigmaj = sigma(m1)
                ELSE
                   normj = rnorms(:,ls,ms2)
                   wwsj = wws(:,ls)
                   n_wsj = n_ws2
                   sigmaj = sigma(m2)
                END IF

                DO ni = 1,n_wsi  !
                   DO nj = 1, n_wsj!
                      y = x*wwsi(ni)*wwsj(nj)/(2*sigmaj)
                      Hsij(2,ni,nj,ci,cj) = Hsij(2,ni,nj,ci,cj) + y * (-mode/ray)*normi(1)
                      Hsij(3,ni,nj,ci,cj) = Hsij(3,ni,nj,ci,cj) + y *   mode/ray *normi(1)
                      Hsij(5,ni,nj,ci,cj) = Hsij(5,ni,nj,ci,cj) + y * (-1/ray)   *normi(1)
                      Hsij(8,ni,nj,ci,cj) = Hsij(8,ni,nj,ci,cj) + y * (-mode/ray)*normi(2)
                      Hsij(9,ni,nj,ci,cj) = Hsij(9,ni,nj,ci,cj) + y *   mode/ray *normi(2)
                      y = x*wwsi(ni)*wwsj(nj)/(2*sigmai)
                      Gsij(2,ni,nj,ci,cj) = Gsij(2,ni,nj,ci,cj) + y * (-mode/ray)*normj(1)
                      Gsij(3,ni,nj,ci,cj) = Gsij(3,ni,nj,ci,cj) + y * ( mode/ray)*normj(1)
                      Gsij(5,ni,nj,ci,cj) = Gsij(5,ni,nj,ci,cj) + y * (-1/ray)   *normj(1)
                      Gsij(8,ni,nj,ci,cj) = Gsij(8,ni,nj,ci,cj) + y * (-mode/ray)*normj(2)
                      Gsij(9,ni,nj,ci,cj) = Gsij(9,ni,nj,ci,cj) + y *   mode/ray *normj(2)
                   ENDDO
                ENDDO
             ENDDO
          END DO
       END DO

       !June 14 2008
       Gsij = c_sym*Gsij
       !June 14 2008

       DO ci = 1, 2
          DO ki = 1, 3
             DO ni = 1, n_wsi
                IF (ci==1) THEN
                   i = interface_H_mu%jjs1(ni,ms)
                ELSE
                   i = interface_H_mu%jjs2(ni,ms)
                END IF
                ib = i + (ki-1)*H_mesh%np
                DO cj = 1, 2
                   DO kj = 1, 3
                      DO nj = 1, n_wsj
                         IF (cj==1) THEN
                            j = interface_H_mu%jjs1(nj,ms)
                         ELSE
                            j = interface_H_mu%jjs2(nj,ms)
                         END IF
                         jb = j + (kj-1)*H_mesh%np
                         DO p = ia(ib),  ia(ib+1) - 1
                            IF (ja(p) == jb) THEN
                               IF  ((ki == 2) .AND. (kj == 1)) THEN
                                  aa1(p) = aa1(p) + Hsij(2,ni,nj,ci,cj)
                                  aa2(p) = aa2(p) + Hsij(3,ni,nj,ci,cj)
                                  EXIT
                               ELSE IF((ki == 1) .AND. (kj == 2)) THEN
                                  aa1(p) = aa1(p) + Gsij(2,ni,nj,ci,cj)
                                  aa2(p) = aa2(p) + Gsij(3,ni,nj,ci,cj)
                                  EXIT
                               ELSEIF  ((ki == 2) .AND. (kj == 2))  THEN  
                                  aa1(p) = aa1(p) + Hsij(5,ni,nj,ci,cj)+Gsij(5,ni,nj,ci,cj)
                                  aa2(p) = aa2(p) + Hsij(5,ni,nj,ci,cj)+Gsij(5,ni,nj,ci,cj)
                                  EXIT
                               ELSEIF ((ki == 2) .AND. (kj == 3)) THEN
                                  aa1(p) = aa1(p) + Hsij(8,ni,nj,ci,cj)
                                  aa2(p) = aa2(p) + Hsij(9,ni,nj,ci,cj)
                                  EXIT
                               ELSEIF ((ki == 3) .AND. (kj == 2)) THEN
                                  aa1(p) = aa1(p) + Gsij(8,ni,nj,ci,cj)
                                  aa2(p) = aa2(p) + Gsij(9,ni,nj,ci,cj)
                                  EXIT
                               ENDIF
                            ENDIF
                         END DO
                      END DO
                   END DO
                END DO
             END DO
          END DO
       END DO

       !terme avec derivees
       Hsij = 0.d0
       Gsij = 0.d0
       DO ls = 1, H_mesh%gauss%l_Gs

          !--------On calcule le rayon du point gauss
          ray = SUM(H_mesh%rr(1,H_mesh%jjs(:,ms2))* H_mesh%gauss%wws(:,ls))
          x = rjs(ls,ms2)*ray

          DO ci = 1, 2
             IF (ci==1) THEN
                normi = rnorms(:,ls,ms1)
                wwsi = w_cs(:,ls)
                dwsi = dw_cs(:,:,ls,ms1)
                n_wsi = n_ws1
                n_wi = n_w1
                sigmai = sigma(m1)
             ELSE
                normi = rnorms(:,ls,ms2)
                wwsi = wws(:,ls)
                dwsi = dw_s(:,:,ls,ms2)
                n_wsi = n_ws2
                n_wi = n_w2
                sigmai = sigma(m2)
             END IF
             DO cj = 1, 2
                IF (cj==1) THEN
                   normj = rnorms(:,ls,ms1)
                   wwsj = w_cs(:,ls)
                   dwsj = dw_cs(:,:,ls,ms1)
                   n_wsj = n_ws1
                   n_wj = n_w1
                   sigmaj = sigma(m1)
                ELSE
                   normj = rnorms(:,ls,ms2)
                   wwsj = wws(:,ls)
                   dwsj = dw_s(:,:,ls,ms2) 
                   n_wsj = n_ws2
                   n_wj = n_w2
                   sigmaj = sigma(m2)
                END IF

                !termes avec derivees
                DO ni = 1,n_wsi
                   DO nj = 1, n_wj
                      y = x*wwsi(ni)/(2*sigmaj)
                      Hsij(1,ni,nj,ci,cj) = Hsij(1,ni,nj,ci,cj) + y*(-dwsj(2,nj))*normi(2)
                      Hsij(4,ni,nj,ci,cj) = Hsij(4,ni,nj,ci,cj) + y*  dwsj(1,nj) *normi(2)
                      Hsij(5,ni,nj,ci,cj) = Hsij(5,ni,nj,ci,cj) + y*(-dwsj(2,nj) *normi(2)-dwsj(1,nj)*normi(1))
                      Hsij(6,ni,nj,ci,cj) = Hsij(6,ni,nj,ci,cj) + y*(-dwsj(1,nj))*normi(1)
                      Hsij(7,ni,nj,ci,cj) = Hsij(7,ni,nj,ci,cj) + y*  dwsj(2,nj) *normi(1)
                   ENDDO
                END DO
                DO ni = 1,n_wi
                   DO nj = 1, n_wsj
                      y = x*wwsj(nj)/(2*sigmai)
                      Gsij(1,ni,nj,ci,cj) = Gsij(1,ni,nj,ci,cj) + y*(-dwsi(2,ni))*normj(2)
                      Gsij(4,ni,nj,ci,cj) = Gsij(4,ni,nj,ci,cj) + y*  dwsi(2,ni) *normj(1)
                      Gsij(5,ni,nj,ci,cj) = Gsij(5,ni,nj,ci,cj) + y*(-dwsi(2,ni) *normj(2)-dwsi(1,ni)*normj(1))
                      Gsij(6,ni,nj,ci,cj) = Gsij(6,ni,nj,ci,cj) + y*(-dwsi(1,ni))*normj(1)
                      Gsij(7,ni,nj,ci,cj) = Gsij(7,ni,nj,ci,cj) + y*  dwsi(1,ni) *normj(2)
                   ENDDO
                END DO

             ENDDO
          ENDDO
       ENDDO

       !June 14 2008
       Gsij = c_sym*Gsij
       !June 14 2008

       DO ci = 1, 2
          DO ki = 1, 3
             DO ni = 1, n_wsi
                IF (ci==1) THEN
                   i = interface_H_mu%jjs1(ni,ms)
                ELSE
                   i = interface_H_mu%jjs2(ni,ms)
                END IF
                ib = i + (ki-1)*H_mesh%np
                DO cj = 1, 2
                   DO kj = 1, 3
                      DO nj = 1, n_wj
                         IF (cj==1) THEN
                            j = H_mesh%jj(nj,m1)
                         ELSE
                            j = H_mesh%jj(nj,m2)
                         END IF
                         jb = j + (kj-1)*H_mesh%np 
                         DO p = ia(ib),  ia(ib+1) - 1
                            IF (ja(p) == jb) THEN
                               IF      ((ki == 1) .AND. (kj == 1)) THEN
                                  aa1(p) = aa1(p) + Hsij(1,ni,nj,ci,cj)
                                  aa2(p) = aa2(p) + Hsij(1,ni,nj,ci,cj)
                                  EXIT
                               ELSEIF  ((ki == 1) .AND. (kj == 3)) THEN
                                  aa1(p) = aa1(p) + Hsij(4,ni,nj,ci,cj)
                                  aa2(p) = aa2(p) + Hsij(4,ni,nj,ci,cj)
                                  EXIT
                               ELSEIF  ((ki == 3) .AND. (kj == 1)) THEN
                                  aa1(p) = aa1(p) + Hsij(7,ni,nj,ci,cj)
                                  aa2(p) = aa2(p) + Hsij(7,ni,nj,ci,cj)
                                  EXIT
                               ELSEIF  ((ki == 2) .AND. (kj == 2)) THEN  
                                  aa1(p) = aa1(p) + Hsij(5,ni,nj,ci,cj)
                                  aa2(p) = aa2(p) + Hsij(5,ni,nj,ci,cj)
                                  EXIT
                               ELSEIF  ((ki == 3) .AND. (kj == 3)) THEN  
                                  aa1(p) = aa1(p) + Hsij(6,ni,nj,ci,cj)   
                                  aa2(p) = aa2(p) + Hsij(6,ni,nj,ci,cj)
                                  EXIT
                               ENDIF
                            ENDIF
                         END DO
                      END DO
                   END DO
                END DO
             END DO
          END DO
       END DO

       DO ci = 1, 2
          DO ki = 1, 3
             DO ni = 1, n_wi
                IF (ci==1) THEN
                   i = H_mesh%jj(ni,m1)
                ELSE
                   i = H_mesh%jj(ni,m2)
                END IF
                ib = i + (ki-1)*H_mesh%np
                DO cj = 1, 2
                   DO kj = 1, 3
                      DO nj = 1, n_wsj
                         IF (cj==1) THEN
                            j = interface_H_mu%jjs1(nj,ms)
                         ELSE
                            j = interface_H_mu%jjs2(nj,ms)
                         END IF
                         jb = j + (kj-1)*H_mesh%np
                         DO p = ia(ib),  ia(ib+1) - 1
                            IF (ja(p) == jb) THEN
                               IF      ((ki == 1) .AND. (kj == 1)) THEN
                                  aa1(p) = aa1(p) + Gsij(1,ni,nj,ci,cj)
                                  aa2(p) = aa2(p) + Gsij(1,ni,nj,ci,cj)
                                  EXIT
                               ELSEIF  ((ki == 1) .AND. (kj == 3)) THEN
                                  aa1(p) = aa1(p) + Gsij(4,ni,nj,ci,cj)
                                  aa2(p) = aa2(p) + Gsij(4,ni,nj,ci,cj)
                                  EXIT
                               ELSEIF  ((ki == 3) .AND. (kj == 1)) THEN
                                  aa1(p) = aa1(p) + Gsij(7,ni,nj,ci,cj)
                                  aa2(p) = aa2(p) + Gsij(7,ni,nj,ci,cj)
                                  EXIT
                               ELSEIF  ((ki == 2) .AND. (kj == 2)) THEN
                                  aa1(p) = aa1(p) + Gsij(5,ni,nj,ci,cj)
                                  aa2(p) = aa2(p) + Gsij(5,ni,nj,ci,cj)
                                  EXIT
                               ELSEIF  ((ki == 3) .AND. (kj == 3)) THEN
                                  aa1(p) = aa1(p) + Gsij(6,ni,nj,ci,cj)
                                  aa2(p) = aa2(p) + Gsij(6,ni,nj,ci,cj)
                                  EXIT
                               ENDIF
                            ENDIF
                         END DO
                      END DO
                   END DO
                END DO
             END DO
          END DO
       END DO

    END DO

    RETURN

  END SUBROUTINE  mat_maxwell_mu

  SUBROUTINE courant_mu(H_mesh,interface_H_mu,sigma,mu_H_field,time,mode,src_H,nl)
    !forcage faisant intervenir J, volumique et interface pour H et phi
    !pour le probleme en entier

    USE def_type_mesh
    USE gauss_points
    USE boundary

    IMPLICIT NONE
    TYPE(mesh_type),                       INTENT(IN)   :: H_mesh
    TYPE(interface_type),                  INTENT(IN)   :: interface_H_mu
    REAL(KIND=8),                          INTENT(IN)   :: time
    REAL(KIND=8), DIMENSION(H_mesh%me),    INTENT(IN)   :: sigma
    REAL(KIND=8), DIMENSION(H_mesh%np),    INTENT(IN)   :: mu_H_field
    INTEGER,                               INTENT(IN)   :: mode
    REAL(KIND=8), DIMENSION(H_mesh%np,6),  INTENT(IN)   :: nl
    REAL(KIND=8), DIMENSION(:,:),          INTENT(INOUT):: src_H

    REAL(KIND=8), DIMENSION(H_mesh%gauss%n_ws,H_mesh%gauss%l_Gs)   :: w_cs   
    REAL(KIND=8), DIMENSION(2) :: normi, gaussp1, gaussp2
    REAL(KIND=8), DIMENSION(H_mesh%gauss%n_ws) :: wwsi 
    REAL(KIND=8) :: j_courant, rot, x, ray
    INTEGER :: i, ni, n, ms, k, ls, n_ws1, n_ws2, ms1, ms2, n_w1, n_w2, m1, m2, ci, ki, n_wsi,n_wi
    REAL(KIND=8), DIMENSION(6)             :: JsolH_anal, test
    REAL(KIND=8) :: muhl1, muhl2, ref, diff
    !April 17th, 2008, JLG
    REAL(KIND=8) :: one
    DATA one/1.d0/
    !April 17th, 2008, JLG

    !**********************************************************************************
    !--------------------TERMES SUR L'INTERFACE SIGMA_MU-------------------------------
    !**********************************************************************************

    WRITE(*,*) 'Assemblage interface_mu'
    CALL gauss(H_mesh)
    n_ws1 = H_mesh%gauss%n_ws
    n_ws2 = H_mesh%gauss%n_ws
    n_w1  = H_mesh%gauss%n_w
    n_w2  = H_mesh%gauss%n_w

    DO ms = 1, interface_H_mu%mes
       ms1 = interface_H_mu%mesh1(ms)
       ms2 = interface_H_mu%mesh2(ms)
       m1 = H_mesh%neighs(ms1)
       m2 = H_mesh%neighs(ms2)

       ref = 1.d-8+SUM((H_mesh%rr(:,H_mesh%jjs(1,ms1)) - H_mesh%rr(:,H_mesh%jjs(2,ms1)))**2)
       diff =      SUM((H_mesh%rr(:,H_mesh%jjs(1,ms1)) - H_mesh%rr(:,H_mesh%jjs(1,ms2)))**2)
       IF (diff/ref .LT. 1.d-10) THEN ! 1 = 1
          w_cs = wws
       ELSE                ! 1 = 2
          DO ls = 1, l_Gs
             w_cs(1,ls)= wws(2,ls)
             w_cs(2,ls)= wws(1,ls)
             IF (n_ws1==3) w_cs(n_ws1,ls) = wws(n_ws1,ls) 
             WRITE(*,*) ' Ouaps! oder of shape functions changed?'
          END DO
       END IF
    END DO

    DO ms = 1, interface_H_mu%mes
       ms2 = interface_H_mu%mesh2(ms)
       ms1 = interface_H_mu%mesh1(ms)
       m2 = H_mesh%neighs(ms2)
       m1 = H_mesh%neighs(ms1)

       DO ls = 1, l_Gs
          !--------On calcule le rayon du point gauss
          ray = SUM(H_mesh%rr(1,H_mesh%jjs(:,ms2))* H_mesh%gauss%wws(:,ls))

          ! Cote 1
          muhl1=SUM(mu_H_field(H_mesh%jjs(:,ms1))*w_cs(:,ls))
          gaussp1(1) = SUM(H_mesh%rr(1,H_mesh%jjs(:,ms1))*w_cs(:,ls))    
          gaussp1(2) = SUM(H_mesh%rr(2,H_mesh%jjs(:,ms1))*w_cs(:,ls))    
          DO k=1, 6
             JsolH_anal(k) = Jexact_gauss(k, gaussp1, mode, one ,sigma(m1), muhl1, time)/sigma(m1) &
                  + muhl1 * SUM(NL(H_mesh%jjs(1:n_ws1,ms1),k)*w_cs(1:n_ws1,ls))
          ENDDO

          ! Cote 2
          muhl2=SUM(mu_H_field(H_mesh%jjs(:,ms2))*wws(:,ls))
          gaussp2(1) = SUM(H_mesh%rr(1,H_mesh%jjs(:,ms2))*wws(:,ls))    
          gaussp2(2) = SUM(H_mesh%rr(2,H_mesh%jjs(:,ms2))*wws(:,ls))    
          IF (MAXVAL(ABS(gaussp1-gaussp2)) > 1.d-11) THEN
             WRITE(*,*) ' BUG courant_mu '
             STOP
          END IF
          DO k=1, 6
             test(k) = Jexact_gauss(k, gaussp2, mode, one ,sigma(m2), muhl2, time)/sigma(m2) &
                  + muhl2 * SUM(NL(H_mesh%jjs(1:n_ws2,ms2),k)*wws(1:n_ws2,ls))
             !JsolH_anal(k) = JsolH_anal(k) + Jexact_gauss(k, gaussp2, mode, one,sigma(m2), muhl2, time)/sigma(m2) &
             !     + muhl2 * SUM(NL(H_mesh%jjs(1:n_ws2,ms2),k)*wws(1:n_ws2,ls))
             !WRITE(*,*) ABS(test(k) - JsolH_anal(k)), test(k), JsolH_anal(k)
             JsolH_anal(k) = JsolH_anal(k) + test(k)
          ENDDO
          !JsolH_anal = JsolH_anal!/2

          !---------forcage interface pour H            
          DO ci = 1, 2
             IF (ci==1) THEN
                normi = rnorms(:,ls,ms1)
                wwsi = w_cs(:,ls)
                n_wsi = n_ws1
             ELSE
                normi = rnorms(:,ls,ms2)
                wwsi = wws(:,ls)
                n_wsi = n_ws2
             END IF
             DO ni = 1, n_wsi
                IF (ci==1) THEN
                   i = interface_H_mu%jjs1(ni,ms)
                ELSE
                   i = interface_H_mu%jjs2(ni,ms)
                END IF
                x = rjs(ls,ms2)*ray*wwsi(ni)/2
                src_H(i,1) = src_H(i,1)+x*(-JsolH_anal(3)*normi(2))
                src_H(i,2) = src_H(i,2)+x*(-JsolH_anal(4)*normi(2))
                src_H(i,3) = src_H(i,3)+x*(JsolH_anal(1)*normi(2)-JsolH_anal(5)*normi(1))
                src_H(i,4) = src_H(i,4)+x*(JsolH_anal(2)*normi(2)-JsolH_anal(6)*normi(1))
                src_H(i,5) = src_H(i,5)+x*(JsolH_anal(3)*normi(1))
                src_H(i,6) = src_H(i,6)+x*(JsolH_anal(4)*normi(1))
             END DO
          ENDDO
       END DO
    END DO

  END SUBROUTINE courant_mu

  SUBROUTINE detect_cavities(Interface, phi_mesh, cav_index)
    ! Written June 7 2008; Jean-Luc Guermond
    ! The purpose of the subroutine is to detect cavities
    ! in the phi mesh
    ! MAXVAL(cav_index) = number of cavities
    ! cav_index(n) = index of the cavity to which the domain of rank n belongs
    USE def_type_mesh

    IMPLICIT NONE
    TYPE(interface_type),           INTENT(IN)  :: INTERFACE
    TYPE(mesh_type),                INTENT(IN)  :: phi_mesh
    INTEGER, POINTER, DIMENSION(:), INTENT(OUT) :: cav_index

    INTEGER, DIMENSION(phi_mesh%me) :: id_test
    INTEGER, DIMENSION(:), ALLOCATABLE :: virgin, dom_next, to_do_next, list_dom_phi, store_dom_phi
    LOGICAL :: cav_test
    INTEGER :: nb_cav, n_dom, nb_next, i_dom, m, i_d_neigh, n, n1, n2, next
    INTEGER, DIMENSION(1) :: loc

    ! I reconstruct the array list_dom_phi (I do not want to read it)
    ALLOCATE(store_dom_phi(MAXVAL(phi_mesh%i_d)))
    n_dom = 1
    store_dom_phi(n_dom) = phi_mesh%i_d(1)
    DO m = 1, phi_mesh%me
       i_dom = phi_mesh%i_d(m)
       IF (MINVAL(ABS(store_dom_phi(1:n_dom) - i_dom)) /= 0) THEN
          n_dom = n_dom + 1
          store_dom_phi(n_dom) = i_dom
       END IF
    END DO
    ALLOCATE(virgin(n_dom), dom_next(n_dom), to_do_next(n_dom), list_dom_phi(n_dom), cav_index(n_dom))
    list_dom_phi = store_dom_phi(1:n_dom)
    DEALLOCATE(store_dom_phi)

    IF (MINVAL(list_dom_phi)==0) THEN
       WRITE(*,*) ' BUG in detect_cavities '
       WRITE(*,*) ' MINVAL(list_dom_phi)==0 '
       WRITE(*,*) ' Do not use 0 as a domain marker '
       STOP
    END IF


    id_test = 0
    id_test(phi_mesh%neighs(INTERFACE%mesh2)) = 1

    nb_cav = 0
    cav_index = 0
    virgin = 0
    DO WHILE (MINVAL(virgin)==0) ! Stop when no more virgin
       DO n_dom = 1, SIZE(list_dom_phi)
          IF (virgin(n_dom)==0) EXIT ! Found one virgin
       END DO
       cav_test = .TRUE.
       dom_next = 0
       nb_next = 1
       dom_next(nb_next) = n_dom
       DO next = 1, SIZE(list_dom_phi)
          IF (dom_next(next)==0) EXIT ! Nothing else to be done
          n_dom = dom_next(next)
          virgin(n_dom) = 1 ! Not virgin now
          i_dom = list_dom_phi(n_dom)
          DO m = 1, phi_mesh%me
             IF (id_test(m) == 1) CYCLE ! Cell is on a H/phi interface
             IF (phi_mesh%i_d(m)/=i_dom) CYCLE ! Cell of no interest
             DO n = 1, 3
                IF (phi_mesh%neigh(n,m)==0) THEN
                   i_d_neigh = 0
                ELSE
                   i_d_neigh = phi_mesh%i_d(phi_mesh%neigh(n,m))
                END IF
                IF (i_d_neigh == i_dom) THEN
                   CYCLE  ! Cell of no interest
                ELSE IF (i_d_neigh == 0) THEN
                   n1 = MOD(n,3) + 1
                   n2 = MOD(n+1,3) + 1
                   IF (ABS(phi_mesh%rr(1,phi_mesh%jj(n1,m))) &
                     + ABS(phi_mesh%rr(1,phi_mesh%jj(n2,m))) > 1.d-10) THEN
                      cav_test = .FALSE. ! This is not a cavity
                   END IF
                ELSE IF (MINVAL(ABS(list_dom_phi - i_d_neigh)) == 0) THEN 
                   IF (MINVAL(to_do_next(1:nb_next) - i_d_neigh) /= 0) THEN !Not in the list yet
                      nb_next = nb_next + 1
                      IF (nb_next > SIZE(list_dom_phi)) THEN
                         WRITE(*,*) ' BUG in detect_cavities, nb_next>SIZE(list_dom_phi)'
                         STOP
                      END IF
                      to_do_next(nb_next) = i_d_neigh
                      loc = MINLOC(ABS(list_dom_phi - i_d_neigh)) !Localize rank in list_dom_phi
                      dom_next(nb_next) = loc(1)
                   END IF
                END IF
             END DO
          END DO
       END DO

       IF (cav_test .EQV. .FALSE.) CYCLE
       nb_cav = nb_cav + 1
       DO next = 1, nb_next
          cav_index(dom_next(next)) = nb_cav
       END DO
    END DO

    DEALLOCATE(virgin, dom_next, to_do_next, list_dom_phi)
    RETURN

  END SUBROUTINE detect_cavities


END MODULE update_maxwell
