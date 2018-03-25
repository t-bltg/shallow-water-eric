!
!Authors Jean-Luc Guermond, Raphael Laguerre, Caroline Nore, Copyrights 2005
!
MODULE subroutine_ns

CONTAINS


  SUBROUTINE three_level_ns(time, dt, Re, list_mode, pp_mesh, &
       vv_mesh, incpn_m1, incpn, pn_m1, pn, un_m1, un , nb_syst_lin, &
       vvrt_per, vvz_per, pp_per,  data_fichier, chmp_mag, second_order_ext_pressure) !, temps)
    !============================== 
    USE def_type_mesh
    USE fem_s_axi_M
    USE fem_s_axi
    USE dir_nodes
    !Jan 29 2007
    USE periodic
    !Jan 29 2007
    USE st_matrix
    USE solve_sp
    USE pardiso_solve
    USE matrix_type
    USE fem_tn_axi
    USE dyn_line
    USE boundary
    USE chaine_caractere
    USE sub_plot
    IMPLICIT NONE
    include 'mpif.h'
    REAL(KIND=8)                                   :: time, dt, Re
    INTEGER,      DIMENSION(:),     INTENT(IN)     :: list_mode   
    TYPE(mesh_type),                INTENT(IN)     :: pp_mesh, vv_mesh
    !jan 29
    TYPE(periodic_type),            INTENT(IN)     :: vvrt_per, vvz_per, pp_per
    !jan 29
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(INOUT)  :: incpn_m1, incpn
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(INOUT)  :: pn_m1, pn
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(INOUT)  :: un_m1, un
    INTEGER,      DIMENSION(:),     INTENT(IN)     :: nb_syst_lin
    CHARACTER(len=64),              INTENT(IN)     :: data_fichier
    LOGICAL,                        INTENT(IN)     :: second_order_ext_pressure

    !    REAL(KIND=8), DIMENSION(:),     OPTIONAL, INTENT(OUT)  :: temps ! Temps
    REAL(KIND=8), DIMENSION(:,:,:),           INTENT(IN)   :: chmp_mag ! champs magnetique

    INTEGER,                                        SAVE         :: m_max_c
    INTEGER,            POINTER,     DIMENSION(:)                :: ia, ja

    TYPE(matrice_bloc), ALLOCATABLE, DIMENSION(:),   SAVE        :: vv_mat                
    TYPE(dyn_real_line),ALLOCATABLE, DIMENSION(:),   SAVE        :: precond_vit    

    TYPE(matrice_bloc), ALLOCATABLE, DIMENSION(:),   SAVE        :: pp_mat
    REAL(KIND=8),       ALLOCATABLE, DIMENSION(:,:), SAVE        :: precond_press 

    TYPE(matrice_bloc),                              SAVE        :: mass_mat
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:),         SAVE        :: precond_mass   
    INTEGER,      ALLOCATABLE, DIMENSION(:,:),       SAVE        :: isave
    LOGICAL,                                         SAVE        :: direct_solver
    !------------DECLARATIONS FOR SOLVE_SPLIB--------------------------------------
    REAL(KIND=8), DIMENSION(4),                      SAVE        :: fpar_sp, fpar_sp_c
    INTEGER,      DIMENSION(4),                      SAVE        :: ipar_sp, ipar_sp_c
    !------------------------------------------------------------------------------
    TYPE(dyn_int_line),                            SAVE          :: pp_js_D ! Dirichlet nodes
    INTEGER,                                       SAVE          :: pp_j_size
    TYPE(dyn_int_line), DIMENSION(3),              SAVE          :: vv_js_D
    INTEGER,            DIMENSION(3),              SAVE          :: vv_j_size
    LOGICAL,                                       SAVE          :: once = .TRUE.
    INTEGER,                                       SAVE          :: ISAVE_MASS, nb_tot_syst_lin
    LOGICAL,                                       SAVE          :: per = .FALSE.
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:),   SAVE          :: rotv_vo
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE,     SAVE          :: rl_half 
    !----------FIN SAVE--------------------------------------------------------------------

    !----------Declaration sans save-------------------------------------------------------
    TYPE(dyn_int_line), DIMENSION(3)                             :: js_D_glob   
    INTEGER                                                      :: np_D_glob
    LOGICAL,      ALLOCATABLE, DIMENSION(:)                      :: DirP ! Type of BC
    LOGICAL,      ALLOCATABLE, DIMENSION(:,:)                    :: DirV ! Type of BC
    INTEGER,      ALLOCATABLE, DIMENSION(:)                      :: list_dirichlet_sides,bord_periodic
    !jan 29 2007
    TYPE(matrice_bloc),        DIMENSION(4)                      :: tampon
    !jan 29 2007

    INTEGER                                  :: i, k, n, n_bord, side1, side2, p, np_v, nu_mat, mode, nb_dirichlet_sides
    INTEGER                                  :: nb_procs, code
    REAL(KIND=8)                             :: moyenne
    !allocations des variables locales
    REAL(KIND=8), DIMENSION(pp_mesh%np, 2)   :: pn_p1, phi, ff_div, div
    REAL(KIND=8), DIMENSION(vv_mesh%np, 2)   :: p_p2, phalf
    REAL(KIND=8), DIMENSION(vv_mesh%np, 6)   :: un_p1, ff, src
    REAL(KIND=8), DIMENSION(2*vv_mesh%np, 3) :: ff_glob
    REAL(KIND=8), DIMENSION(2*vv_mesh%np, 2) :: w
    REAL(KIND=8), DIMENSION(vv_mesh%gauss%l_G*vv_mesh%me,6,SIZE(list_mode)) :: rotv_v, rotb_b

    REAL(KIND=8), DIMENSION(vv_mesh%np,6,SIZE(list_mode)) :: uext
    REAL(KIND=8), DIMENSION(vv_mesh%np) :: vel_loc, vel_tot
    REAL(KIND=8)   :: user_time, tps, dummy, tps_tot, tps_cumul, coeff, coeff1
    LOGICAL :: LES
    !April 17th 2008, JLG
    REAL(KIND=8) :: one, zero, three
    DATA zero, one, three/0.d0,1.d0,3.d0/

    !April 17th 2008, JLG
    !------------------------------END OF DECLARATION--------------------------------------

    ! Calcul du smb par SFT sur points de Gauss
    IF (once) THEN
       once = .FALSE.
       ALLOCATE(rl_half(vv_mesh%gauss%l_G*vv_mesh%me,6))
       ALLOCATE(rotv_vo(vv_mesh%gauss%l_G*vv_mesh%me,6,SIZE(list_mode)))
       rotv_vo = 0.d0

       !-------------DIMENSIONS-------------------------------------------------------
       m_max_c = SIZE(list_mode)
       !------------------------------------------------------------------------------
       !---------------------------------------------------------------------------------   

       OPEN (UNIT=22, FILE ='data', FORM='formatted', STATUS='unknown')
       CALL read_until(22, 'data_type_solver')
       READ(22,*) direct_solver
       CLOSE(22)

       OPEN(UNIT = 21, FILE = data_fichier, FORM = 'formatted', STATUS = 'unknown')    
       IF (.NOT.direct_solver) THEN
          CALL read_until(21, 'data_solver_vit')
          DO k=1,4;      READ(21,*) ipar_sp(k);
          END DO
          DO k=1,2;      READ(21,*) fpar_sp(k);
          END DO

          CALL read_until(21, 'data_solver_press')
          DO k=1,4;      READ(21,*) ipar_sp_c(k);
          END DO
          DO k=1,2;      READ(21,*) fpar_sp_c(k); 
          END DO
       END IF

       !---------PREPARE TYPE OF BOUNDARY CONDITIONS AND js_D ARRAY FOR PRESSURE--------
       CALL read_until(21, 'data_periodic')
       READ  (21, *)  n_bord
       IF (n_bord>0) THEN
          ALLOCATE(bord_periodic(2*n_bord))
          DO n= 1, n_bord
             READ  (21, *)  bord_periodic(2*n-1), bord_periodic(2*n)
          END DO
       END IF

       CALL read_until(21, 'data_condlim_ns')
       !Feb 20 2007
       ALLOCATE (DirP(MAXVAL(pp_mesh%sides)))
       DirP = .FALSE.
       READ(21,*) nb_dirichlet_sides
       IF (nb_dirichlet_sides .LE. 0) THEN
          READ(21,*) 
       ELSE
          ALLOCATE(list_dirichlet_sides(nb_dirichlet_sides))
          READ(21,*) list_dirichlet_sides 
          DO i = 1, nb_dirichlet_sides ! Check coherence of data 
             IF (MINVAL(ABS(pp_mesh%sides - list_dirichlet_sides(i))) /= 0) THEN
                WRITE(*,*) ' BUG list_dirichlet_sides not coherent with pressure mesh'
                STOP
             END IF
          END DO
          DirP(list_dirichlet_sides) = .TRUE.
          ! Verify coeherence with Periodic BCs
          DO n= 1, n_bord
             IF (MINVAL(ABS(list_dirichlet_sides - bord_periodic(2*n-1))) == 0 .OR.  &
                  MINVAL(ABS(list_dirichlet_sides -  bord_periodic(2*n))) == 0) THEN
                WRITE(*,*) ' Donnees Dirichlet et Periodiques incoherentes pour la pression'
                STOP
             END IF
          END DO
          ! End Verify coeherence with Periodic BCs

          DEALLOCATE(list_dirichlet_sides)
       END IF
       !Feb 20 2007
       CALL dirichlet_nodes(pp_mesh%jjs, pp_mesh%sides, DirP, pp_js_D%DIL)
       pp_j_size = SIZE(pp_js_D%DIL)
       DEALLOCATE(DirP)
       !Jan 29 2007
       IF (vvrt_per%n_bord /=0) THEN
          per = .TRUE.
       ELSE
          per = .FALSE.
       END IF
       !Jan 29 2007
       !------------------------------------------------------------------------------

       !---------PREPARE TYPE OF BOUNDARY CONDITIONS AND js_D ARRAY FOR VELOCITY--------
       !Feb 20 2007
       ALLOCATE (DirV(3,MAXVAL(vv_mesh%sides)))
       DirV = .FALSE.
       DO k = 1, 3
          READ(21,*) nb_dirichlet_sides
          IF (nb_dirichlet_sides .LE. 0) THEN
             READ(21,*)
          ELSE
             ALLOCATE(list_dirichlet_sides(nb_dirichlet_sides))
             READ(21,*) list_dirichlet_sides
             DO i = 1, nb_dirichlet_sides ! Check coherence of data
                IF (MINVAL(ABS(vv_mesh%sides - list_dirichlet_sides(i))) /= 0) THEN
                   WRITE(*,*) ' BUG list_dirichlet_sides not coherent with velocity mesh'
                   STOP
                END IF
             END DO
             DirV(k,list_dirichlet_sides) = .TRUE.
             ! Verify coeherence with Periodic BCs
             DO n= 1, n_bord
                IF (MINVAL(ABS(list_dirichlet_sides - bord_periodic(2*n-1) )) == 0 .OR.  &
                     MINVAL(ABS(list_dirichlet_sides - bord_periodic(2*n) )) == 0) THEN
                   WRITE(*,*) ' Donnees Dirichlet et Periodiques incoherentes pour la vitesse'
                   STOP
                END IF
             END DO
             ! End Verify coeherence with Periodic BCs
             DEALLOCATE(list_dirichlet_sides)
          END IF

          !remplissage des tableaux des noeuds contraints, cond Dirichlet
          CALL dirichlet_nodes(vv_mesh%jjs, vv_mesh%sides, DirV(k,:), vv_js_D(k)%DIL)
          vv_j_size(k) = SIZE(vv_js_D(k)%DIL)
       ENDDO
       DEALLOCATE(DirV)
       IF (ALLOCATED(bord_periodic)) DEALLOCATE(bord_periodic)
       !Feb 20 2007

       !calcul de js_D_glob(:)%DIL(:) 
       DO k=1, 3
          IF (k == 1 .OR. k == 2) THEN !pour les composantes V1 et V4, V2 et V3
             np_D_glob = SIZE(vv_js_D(1)%DIL(:)) + SIZE(vv_js_D(2)%DIL(:))
             ALLOCATE(js_D_glob(k)%DIL(np_D_glob))
             js_D_glob(k)%DIL(1:SIZE(vv_js_D(1)%DIL))           = vv_js_D(1)%DIL
             js_D_glob(k)%DIL(SIZE(vv_js_D(1)%DIL)+1:np_D_glob) = vv_js_D(2)%DIL + vv_mesh%np
          ELSEIF (k == 3) THEN !pour les composantes V5 et V6
             np_D_glob = SIZE(vv_js_D(3)%DIL)
             ALLOCATE(js_D_glob(3)%DIL(np_D_glob))
             js_D_glob(3)%DIL  = vv_js_D(3)%DIL
          ENDIF
       ENDDO
       !------------------------------------------------------------------------------

       !-------------MATRIX STRUCTURING-----------------------------------------------
       !structure pour 1 composante
       CALL st_csr(vv_mesh%jj, ia, ja)
       CALL st_csr(pp_mesh%jj, mass_mat%ia, mass_mat%ja)

       !Jan 29 2007
       CALL st_csr_bloc(ia, ja, tampon(1)%ia, tampon(1)%ja, 2)
       ALLOCATE(tampon(2)%ia(SIZE(ia)), tampon(2)%ja(SIZE(ja)))
       tampon(2)%ia = ia
       tampon(2)%ja = ja
       ALLOCATE(tampon(3)%ia(SIZE(mass_mat%ia)), tampon(3)%ja(SIZE(mass_mat%ja)))
       tampon(3)%ia = mass_mat%ia
       tampon(3)%ja = mass_mat%ja
       ALLOCATE(tampon(4)%ia(SIZE(mass_mat%ia)), tampon(4)%ja(SIZE(mass_mat%ja)))
       tampon(4)%ia = mass_mat%ia
       tampon(4)%ja = mass_mat%ja
       !Jan 29 2007

       ALLOCATE(vv_mat(3*m_max_c))
       ALLOCATE(pp_mat(m_max_c))   
       !structure pour 2 composantes
       DO i = 1, m_max_c
          nu_mat = 3*(i-1)+ 1
          IF (i==1) THEN
             !jan 29 2007
             IF (per) THEN
                CALL st_csr_per(vvrt_per, tampon(1)%ia, tampon(1)%ja, vv_mat(1)%ia, vv_mat(1)%ja)
                ALLOCATE(vv_mat(2)%ia(SIZE(vv_mat(1)%ia)), vv_mat(2)%ja(SIZE(vv_mat(1)%ja)))
                vv_mat(2)%ia = vv_mat(1)%ia
                vv_mat(2)%ja = vv_mat(1)%ja
                CALL st_csr_per(vvz_per,  tampon(2)%ia, tampon(2)%ja, vv_mat(3)%ia, vv_mat(3)%ja)
                CALL st_csr_per(pp_per,   tampon(3)%ia, tampon(3)%ja, pp_mat(1)%ia, pp_mat(1)%ja)
                CALL st_csr_per(pp_per,   tampon(3)%ia, tampon(3)%ja, mass_mat%ia, mass_mat%ja)
             ELSE
                ALLOCATE(vv_mat(1)%ia(SIZE(tampon(1)%ia)), vv_mat(1)%ja(SIZE(tampon(1)%ja)))
                vv_mat(1)%ia = tampon(1)%ia
                vv_mat(1)%ja = tampon(1)%ja

                ALLOCATE(vv_mat(2)%ia(SIZE(tampon(1)%ia)), vv_mat(2)%ja(SIZE(tampon(1)%ja)))
                vv_mat(2)%ia = tampon(1)%ia
                vv_mat(2)%ja = tampon(1)%ja

                !structure pour le vitesse vz
                ALLOCATE(vv_mat(3)%ia(SIZE(tampon(2)%ia)), vv_mat(3)%ja(SIZE(tampon(2)%ja)))
                vv_mat(3)%ia = tampon(2)%ia
                vv_mat(3)%ja = tampon(2)%ja

                !structure pour la pression
                ALLOCATE(pp_mat(1)%ia(SIZE(tampon(3)%ia)), pp_mat(1)%ja(SIZE(tampon(3)%ja)))
                pp_mat(1)%ia = tampon(3)%ia
                pp_mat(1)%ja = tampon(3)%ja
             END IF

          ELSE
             ALLOCATE(vv_mat(nu_mat)%ia(SIZE(vv_mat(1)%ia)),&
                  vv_mat(nu_mat)%ja(SIZE(vv_mat(1)%ja)))
             vv_mat(nu_mat)%ia = vv_mat(1)%ia
             vv_mat(nu_mat)%ja = vv_mat(1)%ja

             nu_mat = nu_mat + 1
             ALLOCATE(vv_mat(nu_mat)%ia(SIZE(vv_mat(1)%ia)),&
                  vv_mat(nu_mat)%ja(SIZE(vv_mat(1)%ja)))
             vv_mat(nu_mat)%ia = vv_mat(1)%ia
             vv_mat(nu_mat)%ja = vv_mat(1)%ja

             !structure pour le vitesse vz
             nu_mat = nu_mat + 1
             ALLOCATE(vv_mat(nu_mat)%ia(SIZE(vv_mat(3)%ia)),&
                  vv_mat(nu_mat)%ja(SIZE(vv_mat(3)%ja)))
             vv_mat(nu_mat)%ia = vv_mat(3)%ia
             vv_mat(nu_mat)%ja = vv_mat(3)%ja

             !structure pour la pression
             ALLOCATE(pp_mat(i)%ia(SIZE(pp_mat(1)%ia)),&
                  pp_mat(i)%ja(SIZE(pp_mat(1)%ja)))
             pp_mat(i)%ia = pp_mat(1)%ia
             pp_mat(i)%ja = pp_mat(1)%ja

          END IF
       END DO

       DEALLOCATE(ia,ja) 
       !---------------------------------------------------------------------------------  


       CALL read_until(21, 'data_stab_LES_NS')
       READ(21,*) LES
       !JLG+CN, Nov 1, 2010: Bug corrected. 
       IF (LES) THEN
          READ(21,*) coeff1
       ELSE
          coeff1=0.d0
       END IF
       !JLG+CN, Nov 1, 2010: Bug corrected. 
       !--------------------------------------------------------------------------------

       CLOSE(21) 

       !-----------CALCUL DES OPERATEURS --------------------------------------------
       !mass_mat : operateur de masse pour le calcul de la divergence, defini sur pp_mesh
       !vv_mat: operateurs relatifs aux calculs sur la vitesse, il y a trois operateurs
       !        par mode.
       !pp_mat: operateurs relatifs aux calculs sur la pression. Il y a un operateur par mode 
       !        le tableau pp_mat d'operateur commencant a 1. Il ya donc mod_max+1 operateurs 
       !        pour la pression.   
       !Preconditionnement : chaque operateur est preconditionne en divisant toutes les valeurs
       !                     de chaque ligne par la somme de la ligne correspondante, c'est l'objet
       !                     de la procedure  'precond_mat'
       !----------------------------------------------------------------------------------


       !-------operateur de masse pour le calcul de la divergence
       ALLOCATE (mass_mat%aa(SIZE(mass_mat%ja)))
       !jan 29 2007
       ALLOCATE (tampon(4)%aa(SIZE(tampon(4)%ja)))
       tampon(4)%aa = 0.d0
       !jan 29 2007
       CALL qs_00_M(pp_mesh, one, tampon(4)%ia, tampon(4)%ja, tampon(4)%aa)
       ALLOCATE(precond_mass(pp_mesh%np))
       IF (per) THEN
          CALL bc_per_M(pp_per, tampon(4), mass_mat)
       ELSE
          mass_mat%aa = tampon(4)%aa
       END IF

       CALL precond_mat(mass_mat%ia, mass_mat%aa, precond_mass)
       !----------------------------------------------------------------------------------


       !-------allocations des operateurs pour la vitesse et la pression
       ALLOCATE(precond_vit(3*m_max_c)) 
       !Jan 29 2007
       ALLOCATE (tampon(1)%aa(SIZE(tampon(1)%ja)))
       ALLOCATE (tampon(2)%aa(SIZE(tampon(2)%ja)))
       !Jan 29 2007
       nu_mat = 0
       DO i = 1, m_max_c 
          nu_mat = nu_mat + 1
          ALLOCATE(precond_vit(nu_mat)%DRL(2*vv_mesh%np)) 
          nu_mat = nu_mat + 1
          ALLOCATE(precond_vit(nu_mat)%DRL(2*vv_mesh%np)) 
          nu_mat = nu_mat + 1
          ALLOCATE(precond_vit(nu_mat)%DRL(vv_mesh%np)) 
       END DO
       ALLOCATE(precond_press(pp_mesh%np,m_max_c))
       ALLOCATE(ISAVE(4,m_max_c))

       !calcul des operateurs
       DO i=1, m_max_c

          mode = list_mode(i)

          !vitesse
          DO k= 1, 3
             nu_mat = 3*(i-1) + k
             ALLOCATE(vv_mat(nu_mat)%aa(SIZE(vv_mat(nu_mat)%ja)))
             !jan 29 2007
             IF (k<3) THEN
                tampon(1)%aa = 0.d0
                CALL qs_mass_diff_vect(k, vv_mesh, one/Re, three/(2*dt), coeff1, mode , tampon(1)%ia, &
                     tampon(1)%ja, tampon(1)%aa)
                !CALL qs_mass_diff_vect(k, vv_mesh, one/Re, three/(2*dt), zero, mode , tampon(1)%ia, &
                !     tampon(1)%ja, tampon(1)%aa)


                IF (per) THEN
                   CALL bc_per_M(vvrt_per, tampon(1), vv_mat(nu_mat))
                ELSE
                   vv_mat(nu_mat)%aa = tampon(1)%aa
                END IF
             ELSE
                tampon(2)%aa = 0.d0
                CALL qs_mass_diff_vect(k, vv_mesh, one/Re, three/(2*dt), coeff1, mode , tampon(2)%ia, &
                     tampon(2)%ja, tampon(2)%aa)
                ! CALL qs_mass_diff_vect(k, vv_mesh, one/Re, three/(2*dt), zero, mode , tampon(2)%ia, &
                !     tampon(2)%ja, tampon(2)%aa)

                IF (per) THEN
                   CALL bc_per_M(vvz_per, tampon(2), vv_mat(nu_mat))
                ELSE
                   vv_mat(nu_mat)%aa = tampon(2)%aa
                END IF
             END IF
             !jan 29 2007
             CALL precond_mat(vv_mat(nu_mat)%ia, vv_mat(nu_mat)%aa, precond_vit(nu_mat)%DRL) 
             CALL Dirichlet_M (js_D_glob(k)%DIL, vv_mat(nu_mat)%ia, vv_mat(nu_mat)%ja, &
                  vv_mat(nu_mat)%aa)
          ENDDO

          !pression
          ALLOCATE(pp_mat(i)%aa(SIZE(pp_mat(i)%ja)))
          !Jan 29 2007
          !pp_mat(i)%aa   = 0.d0
          !CALL qs_mass_diff_vect(3, pp_mesh, one, 0.d0, mode, pp_mat(i)%ia, &
          !     pp_mat(i)%ja, pp_mat(i)%aa(:))
          !IF (mode == 0) THEN !regularisation de la pression
          !   k = pp_mesh%np/3
          !   DO p = 1, pp_mat(i)%ia(k), pp_mat(i)%ia(k+1) - 1
          !      IF (pp_mat(i)%ja(p)==k) THEN
          !         pp_mat(i)%aa(p) = 1.d30
          !      ELSE
          !         pp_mat(i)%aa(p) = 0.d0
          !      END IF
          !   END DO
          !END IF
          !jan 29 2007
          ALLOCATE (tampon(3)%aa(SIZE(tampon(3)%ja)))
          tampon(3)%aa = 0.d0
          !jan 29 2007
          CALL qs_mass_diff_vect(3, pp_mesh, one, zero, zero, mode, tampon(3)%ia, &
               tampon(3)%ja, tampon(3)%aa(:))
          IF (per) THEN
             CALL bc_per_M(pp_per, tampon(3), pp_mat(i))
          ELSE
             pp_mat(i)%aa = tampon(3)%aa
          END IF
          IF (mode ==0  .AND. SIZE(pp_js_D%DIL)==0) THEN !regularisation de la pression
             k = pp_mesh%np/3
             DO p = 1, pp_mat(i)%ia(k), pp_mat(i)%ia(k+1) - 1
                IF (pp_mat(i)%ja(p)==k) THEN
                   pp_mat(i)%aa(p) = 1.d30
                ELSE
                   pp_mat(i)%aa(p) = 0.d0
                END IF
             END DO
          END IF
          !Jan 29 2007
          CALL precond_mat(pp_mat(i)%ia, pp_mat(i)%aa, precond_press(:,i))
          ISAVE(1,i) = -(4*(i-1) + 1)
          ISAVE(2,i) = -(4*(i-1) + 2)
          ISAVE(3,i) = -(4*(i-1) + 3)
          ISAVE(4,i) = -(4*(i-1) + 4)
          IF (SIZE(pp_js_D%DIL)/=0) THEN
             CALL Dirichlet_M (pp_js_D%DIL, pp_mat(i)%ia, pp_mat(i)%ja, &
                  pp_mat(i)%aa)
          END IF
       ENDDO
       ISAVE = ISAVE - nb_syst_lin(2)
       ISAVE_MASS = -(4*m_max_c + 1) - nb_syst_lin(2)
       nb_tot_syst_lin = nb_syst_lin(1)

       DO k = 1, 3
          DEALLOCATE(js_D_glob(k)%DIL)
       END DO
       !Jan 29 2007
       DO k = 1, 4
          DEALLOCATE(tampon(k)%ia,tampon(k)%ja,tampon(k)%aa)
       END DO
       !Jan 29 2007
    ENDIF

    tps_tot = user_time(dummy)
    tps_cumul = 0

    !Calcul du smb par SFT sur les points de Gauss
    tps = user_time(dummy)
    uext = 2*un-un_m1
    CALL smb_cross_prod_gauss_sft_par(vv_mesh,list_mode,uext,rotv_v)
    IF (SIZE(chmp_mag,1)>1) THEN
       CALL smb_cross_prod_gauss_sft_par(vv_mesh,list_mode,chmp_mag,rotb_b)
       rotv_v = rotv_v - rotb_b
    END IF
    tps = user_time(dummy) - tps; tps_cumul=tps_cumul+tps
    !WRITE(*,*) ' Tps fft vitesse', tps
    !------------DEBUT BOUCLE SUR LES MODES----------------

    vel_loc = 0.d0
    DO i = 1, m_max_c
       IF (list_mode(i)==0) THEN 
          coeff = 1.d0  
       ELSE
          coeff = .5d0
       END IF
       vel_loc = vel_loc + coeff*(un(:,1,i)**2+un(:,2,i)**2+un(:,3,i)**2+un(:,4,i)**2+un(:,5,i)**2+un(:,6,i)**2)
    END DO

    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,code)
    CALL MPI_ALLREDUCE(vel_loc,vel_tot,vv_mesh%np,MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, code)

    vel_tot = sqrt(vel_tot)
    !write(*,*) ' MAXVAL(VEL_TOT) : ' , MAXVAL(VEL_TOT)

    DO i = 1, m_max_c
       mode = list_mode(i)
       !-------------CALCUL DU SEND MEMBRE VITESSE-------------------------------------
       !Calcul de phi
       tps = user_time(dummy)
       !jan 29 2007
       IF (second_order_ext_pressure) THEN
          pn_p1(:,:) = 2*pn(:,:,i) -pn_m1(:,:,i) 
       ELSE
          pn_p1(:,:) = pn(:,:,i)
       END IF
       phi = pn_p1(:,:) + (4.d0 * incpn(:,:,i) - incpn_m1(:,:,i))/3.d0
       !jan 29 2007

       !Injection de la pression P1 -> P2  
       DO k = 1, 2
          CALL inject(pp_mesh%jj, vv_mesh%jj, phi(:,k), p_p2(:,k))
       ENDDO

       !Etape de prediction
       DO k=1,6
          src(:,k) = ff_source(k, vv_mesh%rr, mode, time, Re, 'ns')	
       END DO

       !November 6, 2008
       !We try div stabilization
       !CALL qs_ns_2006(vv_mesh, mode, src, &
       !     (4*un(:,:,i)-un_m1(:,:,i))/(2*dt), p_p2(:,:), dt, ff, rotv_v(:,:,i))
       !CALL qs_ns_stab_2008(vv_mesh, mode, src, vel_tot, &
       !     (4*un(:,:,i)-un_m1(:,:,i))/(2*dt), uext(:,:,i), p_p2(:,:), dt, ff, rotv_v(:,:,i))
       rl_half = (rotv_v(:,:,i)+ rotv_vo(:,:,i))/2
       rotv_vo(:,:,i) = rotv_v(:,:,i)
       DO k = 1, 2
          CALL inject(pp_mesh%jj, vv_mesh%jj, (pn(:,k,i) + pn_m1(:,k,i))/2, phalf(:,k))
       ENDDO

       !JLG: JAN 26 2010
       CALL qs_ns_stab_new(vv_mesh, mode, src, vel_tot, &
            (4*un(:,:,i)-un_m1(:,:,i))/(2*dt), uext(:,:,i), p_p2(:,:),(un(:,:,i)-un_m1(:,:,i))/dt-src, phalf, &
            rl_half, dt, ff, rotv_v(:,:,i))
       !JLG: JAN 26 2010
       !CALL qs_ns_stab_2010(vv_mesh, pp_mesh, mode, src, vel_tot, &
       !    (4*un(:,:,i)-un_m1(:,:,i))/(2*dt), uext(:,:,i), p_p2(:,:),(un(:,:,i)-un_m1(:,:,i))/dt-src, phalf, &
       !     rl_half, dt, ff, rotv_v(:,:,i))
       !November 6, 2008

       np_v = vv_mesh%np
       !Changement de format
       ff_glob(1:np_v, 1) = ff(:,1)
       ff_glob(np_v+1:,1) = ff(:,4)
       ff_glob(1:np_v, 2) = ff(:,2)
       ff_glob(np_v+1:,2) = ff(:,3)
       ff_glob(1:np_v, 3) = ff(:,5)
       ff_glob(np_v+1:,3) = ff(:,6)

       !Jan 29 2007
       IF (per) THEN 
          CALL bc_per(vvrt_per, ff_glob(:,1))
          CALL bc_per(vvrt_per, ff_glob(:,2))
          CALL bc_per(vvz_per,  ff_glob(1:np_v,3))
          CALL bc_per(vvz_per,  ff_glob(np_v+1:,3))
       END IF
       !Jan 29 2007

       !Preconditionnement
       nu_mat = 3*(i-1) + 1
       ff_glob(:, 1) = ff_glob(:,1) * precond_vit(nu_mat)%DRL
       nu_mat = nu_mat + 1
       ff_glob(:, 2) = ff_glob(:,2) * precond_vit(nu_mat)%DRL
       nu_mat = nu_mat + 1
       ff_glob(1:np_v, 3) = ff_glob(1:np_v, 3) * precond_vit(nu_mat)%DRL
       ff_glob(np_v+1:,3) = ff_glob(np_v+1:,3) * precond_vit(nu_mat)%DRL

       !Conditions aux limites
       ff_glob(vv_js_D(1)%DIL, 1)      = vv_exact(1, vv_mesh%rr(:,vv_js_D(1)%DIL), mode, time)
       ff_glob(vv_js_D(2)%DIL+np_v, 1) = vv_exact(4, vv_mesh%rr(:,vv_js_D(2)%DIL), mode, time)
       ff_glob(vv_js_D(1)%DIL, 2)      = vv_exact(2, vv_mesh%rr(:,vv_js_D(1)%DIL), mode, time)
       ff_glob(vv_js_D(2)%DIL+np_v, 2) = vv_exact(3, vv_mesh%rr(:,vv_js_D(2)%DIL), mode, time)
       ff_glob(vv_js_D(3)%DIL, 3)      = vv_exact(5, vv_mesh%rr(:,vv_js_D(3)%DIL), mode, time)
       ff_glob(vv_js_D(3)%DIL+np_v, 3) = vv_exact(6, vv_mesh%rr(:,vv_js_D(3)%DIL), mode, time)

       !adaptation au format pour second membre et extrapolation
       w(:np_v,  1) = uext(:,1,i)
       w(np_v+1:,1) = uext(:,4,i)
       w(:np_v,  2) = uext(:,2,i)
       w(np_v+1:,2) = uext(:,3,i)
       un_p1(:,5)   = uext(:,5,i)
       un_p1(:,6)   = uext(:,6,i)

       !Enforce boundary conditions on initial guess
       w(vv_js_D(1)%DIL, 1)          = vv_exact(1, vv_mesh%rr(:,vv_js_D(1)%DIL), mode, time)
       w(vv_js_D(2)%DIL+np_v, 1)     = vv_exact(4, vv_mesh%rr(:,vv_js_D(2)%DIL), mode, time)
       w(vv_js_D(1)%DIL, 2)          = vv_exact(2, vv_mesh%rr(:,vv_js_D(1)%DIL), mode, time)
       w(vv_js_D(2)%DIL+np_v, 2)     = vv_exact(3, vv_mesh%rr(:,vv_js_D(2)%DIL), mode, time)
       un_p1(vv_js_D(3)%DIL, 5)      = vv_exact(5, vv_mesh%rr(:,vv_js_D(3)%DIL), mode, time)
       un_p1(vv_js_D(3)%DIL, 6) = vv_exact(6, vv_mesh%rr(:,vv_js_D(3)%DIL), mode, time)

       tps = user_time(dummy) - tps; tps_cumul=tps_cumul+tps
       !WRITE(*,*) ' Tps second membre vitesse', tps
       !------------------------------------------------------------------------------------- 


       !--------------------INVERSION DE L'OPERATEUR 1 --------------
       tps = user_time(dummy)
       !Inversion operateur 1
       nu_mat = 3*(i-1) + 1
       IF (.NOT.direct_solver) THEN
          CALL solve_splib(ipar_sp, fpar_sp, vv_mat(nu_mat), &
               ff_glob(:,1), w(:,1), ISAVE(1,i), nb_tot_syst_lin)
       ELSE
          CALL solve_pardiso(vv_mat(nu_mat)%aa, vv_mat(nu_mat)%ia, vv_mat(nu_mat)%ja,  &
               ff_glob(:,1), w(:,1), ISAVE(1,i), nb_tot_syst_lin)
       END IF
       ISAVE(1,i) = ABS(ISAVE(1,i))
       !Inversion operateur 2
       nu_mat  = nu_mat + 1
       IF (.NOT.direct_solver) THEN
          CALL solve_splib(ipar_sp, fpar_sp, vv_mat(nu_mat), &
               ff_glob(:,2), w(:,2), ISAVE(2,i))
       ELSE
          CALL solve_pardiso(vv_mat(nu_mat)%aa, vv_mat(nu_mat)%ia, vv_mat(nu_mat)%ja, &
               ff_glob(:,2), w(:,2), ISAVE(2,i))
       END IF

       ISAVE(2,i) = ABS(ISAVE(2,i))

       !Inversion operateur 3
       nu_mat = nu_mat + 1
       IF (.NOT.direct_solver) THEN
          CALL solve_splib(ipar_sp, fpar_sp, vv_mat(nu_mat), &
               ff_glob(:np_v,  3), un_p1(:,5), ISAVE(3,i))
       ELSE
          CALL solve_pardiso(vv_mat(nu_mat)%aa, vv_mat(nu_mat)%ia, vv_mat(nu_mat)%ja, &
               ff_glob(:np_v,  3), un_p1(:,5), ISAVE(3,i))
       END IF
       ISAVE(3,i) = ABS(ISAVE(3,i))
       IF (.NOT.direct_solver) THEN
          CALL solve_splib(ipar_sp, fpar_sp, vv_mat(nu_mat), &
               ff_glob(np_v+1:,3), un_p1(:,6), ISAVE(3,i))
       ELSE
          CALL solve_pardiso(vv_mat(nu_mat)%aa, vv_mat(nu_mat)%ia, vv_mat(nu_mat)%ja, &
               ff_glob(np_v+1:,3), un_p1(:,6), ISAVE(3,i))
       END IF

       !ecriture de la solution ds l'autre format  
       un_p1(:,1) =  w(:np_v,  1)
       un_p1(:,4) =  w(np_v+1:,1)
       un_p1(:,2) =  w(:np_v,  2)
       un_p1(:,3) =  w(np_v+1:,2)
       tps = user_time(dummy) - tps; tps_cumul=tps_cumul+tps
       !WRITE(*,*) ' Tps solution des pb de vitesse', tps, 'for mode ', mode
       !------------------------------------------------------------------------------------- 


       !---------------RESOLUTION DE L'EQUATION DE POISSON--------------
       !------------on resout -LAP(PH3) = -3/(2*dt)*DIV(un_p1)
       !Inversion d'un seul operateur et calcul de la divergence de la vitesse predite
       tps = user_time(dummy)
       CALL qs_01_div_hybrid_2006(vv_mesh, pp_mesh, mode, un_p1 , div)

       !Jan 29 2007       
       IF (per) THEN
          CALL bc_per(pp_per, div(:,1))
          CALL bc_per(pp_per, div(:,2))
       END IF
       !Jan 29 2007       

       ff_div(:,1) = -(1.5d0/dt)*div(:,1)*precond_press(:,i)
       ff_div(:,2) = -(1.5d0/dt)*div(:,2)*precond_press(:,i)
       IF (.NOT.direct_solver) THEN
          CALL solve_splib(ipar_sp_c, fpar_sp_c, pp_mat(i), ff_div(:,1), phi(:,1), ISAVE(4,i))
       ELSE
          CALL solve_pardiso(pp_mat(i)%aa,pp_mat(i)%ia,pp_mat(i)%ja, ff_div(:,1), phi(:,1), ISAVE(4,i))
       END IF

       ISAVE(4,i) = ABS(ISAVE(4,i))
       IF (.NOT.direct_solver) THEN
          CALL solve_splib(ipar_sp_c, fpar_sp_c, pp_mat(i), ff_div(:,2), phi(:,2), ISAVE(4,i))
       ELSE
          CALL solve_pardiso(pp_mat(i)%aa,pp_mat(i)%ia,pp_mat(i)%ja, ff_div(:,2), phi(:,2), ISAVE(4,i))
       END IF
       tps = user_time(dummy) - tps; tps_cumul=tps_cumul+tps
       !WRITE(*,*) ' Tps solution des pb de pression', tps, 'for mode ', mode
       !------------------------------------------------------------------------------------- 


       !---------------CORRECTION DE LA PRESSION-----------------------
       tps = user_time(dummy)
       DO k=1, 2     !2 composantes (cos et sin)  pour la divergence
          ff_div(:,k) = precond_mass * div(:,k)
          !Jan 29 2007       
          IF (per) THEN
             CALL bc_per(pp_per, ff_div(:,k))
          END IF
          IF (SIZE(pp_js_D%DIL)/=0) THEN
             ff_div(pp_js_D%DIL,k) = 0.d0
          END IF
          !Jan 29 2007       
          IF (.NOT.direct_solver) THEN
             CALL solve_splib(ipar_sp, fpar_sp, mass_mat, ff_div(:,k), div(:,k), ISAVE_MASS)
          ELSE
             CALL solve_pardiso(mass_mat%aa,mass_mat%ia,mass_mat%ja,ff_div(:,k), div(:,k), ISAVE_MASS)
          END IF
          ISAVE_MASS = ABS(ISAVE_MASS)
          !calcul de la pression
          !jan 29 2007
          pn_p1(:,k) = pn_p1(:,k) + phi(:,k) - div(:,k)/Re
          !jan 29 2007
       ENDDO
       tps = user_time(dummy) - tps; tps_cumul=tps_cumul+tps
       !WRITE(*,*) ' Tps correction de divergence', tps
       !------------------------------------------------------------------------------------- 


       !---------------UPDATES-----------------------
       tps = user_time(dummy)
       !rehaussement par la  pression moyenne
       IF (SIZE(pp_js_D%DIL)/=0) THEN
          pn_p1(pp_js_D%DIL,1) = pp_exact(1,pp_mesh%rr(:,pp_js_D%DIL),mode,time) 
          pn_p1(pp_js_D%DIL,2) = pp_exact(2,pp_mesh%rr(:,pp_js_D%DIL),mode,time) 
       ELSE IF (mode == 0)  THEN  
          CALL Moy(pp_mesh, pn_p1(:,1),moyenne)
          pn_p1(:,1) = pn_p1(:,1)-moyenne
       ENDIF

       !JLG AR, Dec 18 2008/JLG Bug corrige Jan 23 2010
       IF (mode==0) THEN
          un_p1 (:,2) = 0.d0
          un_p1 (:,4) = 0.d0
          un_p1 (:,6) = 0.d0
          pn_p1 (:,2) = 0.d0
       END IF
       !JLG AR, Dec 18 2008/JLG Bug corrige Jan 23 2010 

       un_m1(:,:,i)    = un (:,:,i)
       un   (:,:,i)    = un_p1 

       pn_m1(:,:,i)    = pn(:,:,i) 
       pn   (:,:,i)    = pn_p1

       incpn_m1(:,:,i) = incpn(:,:,i) 
       incpn   (:,:,i) = phi

       tps = user_time(dummy) - tps; tps_cumul=tps_cumul+tps
       !WRITE(*,*) ' Tps  des updates', tps
       !-------------------------------------------------------------------------------------
    ENDDO

    tps_tot = user_time(dummy) - tps_tot
    !WRITE(*,'(A,2(f13.3,2x))') '  Tps boucle en temps Navier_stokes', tps_tot, tps_cumul
    !WRITE(*,*) ' TIME = ', time, '========================================'

  END SUBROUTINE three_level_ns
  !============================================

  SUBROUTINE inject(jj_c, jj_f, pp_c, pp_f)

    IMPLICIT NONE

    INTEGER,      DIMENSION(:,:), INTENT(IN)  :: jj_c, jj_f
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: pp_c
    REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: pp_f
    REAL(KIND=8) :: half = 0.5
    INTEGER:: m
    IF (SIZE(jj_c,1)==3) THEN
       DO m = 1, SIZE(jj_f,2)
          !         pp_f(jj_f(4,:)) = (pp_c(jj_c(2,:)) + pp_c(jj_c(3,:)))*half
          !         pp_f(jj_f(5,:)) = (pp_c(jj_c(3,:)) + pp_c(jj_c(1,:)))*half
          !         pp_f(jj_f(6,:)) = (pp_c(jj_c(1,:)) + pp_c(jj_c(2,:)))*half 
          pp_f(jj_f(1:3,m)) =  pp_c(jj_c(:,m))
          pp_f(jj_f(4,m)) = (pp_c(jj_c(2,m)) + pp_c(jj_c(3,m)))*half
          pp_f(jj_f(5,m)) = (pp_c(jj_c(3,m)) + pp_c(jj_c(1,m)))*half
          pp_f(jj_f(6,m)) = (pp_c(jj_c(1,m)) + pp_c(jj_c(2,m)))*half 
       END DO

    ELSE         
       pp_f(jj_f(1:4,m)) =  pp_c(jj_c(:,m))
       pp_f(jj_f(5,:)) = (pp_c(jj_c(3,:)) + pp_c(jj_c(4,:)))*half
       pp_f(jj_f(6,:)) = (pp_c(jj_c(4,:)) + pp_c(jj_c(2,:)))*half 
       pp_f(jj_f(7,:)) = (pp_c(jj_c(2,:)) + pp_c(jj_c(3,:)))*half 
       pp_f(jj_f(8,:)) = (pp_c(jj_c(1,:)) + pp_c(jj_c(4,:)))*half 
       pp_f(jj_f(9,:)) = (pp_c(jj_c(3,:)) + pp_c(jj_c(1,:)))*half 
       pp_f(jj_f(10,:)) = (pp_c(jj_c(1,:)) + pp_c(jj_c(2,:)))*half

    END IF

  END SUBROUTINE inject

  SUBROUTINE inject_champ(jj_c, jj_f, pp_c, pp_f)

    IMPLICIT NONE

    INTEGER,      DIMENSION(:,:), INTENT(IN)  :: jj_c, jj_f
    REAL(KIND=8), DIMENSION(:,:,:),   INTENT(IN)  :: pp_c
    REAL(KIND=8), DIMENSION(:,:,:),   INTENT(OUT) :: pp_f
    REAL(KIND=8) :: half = 0.5
    INTEGER:: m
    IF (SIZE(jj_c,1)==3) THEN
       DO m = 1, SIZE(jj_f,2) !Dimension deux
          pp_f(jj_f(1:3,m),:,:) =  pp_c(jj_c(:,m),:,:)
          pp_f(jj_f(4,m),  :,:) = (pp_c(jj_c(2,m),:,:) + pp_c(jj_c(3,m),:,:))*half
          pp_f(jj_f(5,m),  :,:) = (pp_c(jj_c(3,m),:,:) + pp_c(jj_c(1,m),:,:))*half
          pp_f(jj_f(6,m),  :,:) = (pp_c(jj_c(1,m),:,:) + pp_c(jj_c(2,m),:,:))*half 
       END DO

    ELSE  !Dimension trois
       pp_f(jj_f(1:4,m),:,:) =  pp_c(jj_c(:,m),:,:)
       pp_f(jj_f(5,:),:,:) = (pp_c(jj_c(3,:),:,:) + pp_c(jj_c(4,:),:,:))*half
       pp_f(jj_f(6,:),:,:) = (pp_c(jj_c(4,:),:,:) + pp_c(jj_c(2,:),:,:))*half 
       pp_f(jj_f(7,:),:,:) = (pp_c(jj_c(2,:),:,:) + pp_c(jj_c(3,:),:,:))*half 
       pp_f(jj_f(8,:),:,:) = (pp_c(jj_c(1,:),:,:) + pp_c(jj_c(4,:),:,:))*half 
       pp_f(jj_f(9,:),:,:) = (pp_c(jj_c(3,:),:,:) + pp_c(jj_c(1,:),:,:))*half 
       pp_f(jj_f(10,:),:,:) = (pp_c(jj_c(1,:),:,:) + pp_c(jj_c(2,:),:,:))*half

    END IF

  END SUBROUTINE inject_champ

  SUBROUTINE project_champ(pp_f, pp_c, mesh_f, mesh_c)

    USE def_type_mesh
    IMPLICIT NONE
    TYPE(mesh_type),                  INTENT(IN)  :: mesh_f, mesh_c
    REAL(KIND=8), DIMENSION(:,:,:),   INTENT(IN)  :: pp_f
    REAL(KIND=8), DIMENSION(:,:,:),   INTENT(OUT) :: pp_c
    INTEGER :: m

    DO m =1, mesh_c%me
       pp_c(mesh_c%jj(1:3,m),:,:) = pp_f(mesh_f%jj(1:3,m),:,:)
    END DO

  END SUBROUTINE project_champ


  SUBROUTINE smb_cross_prod_gauss_sft_par(mesh,list_mode,V_in,V_out)
    !=================================
    !Calcul du second membre pour le probleme de Stokes
    !sans terme non lineaire et source (ff)

    USE Gauss_points     
    USE sft_parallele

    IMPLICIT NONE

    TYPE(mesh_type), TARGET                                  :: mesh
    INTEGER,      DIMENSION(:),     INTENT(IN)  :: list_mode    
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: V_in
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT) :: V_out

    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%me,6,SIZE(list_mode)) :: RotV, W
    INTEGER,      DIMENSION(mesh%gauss%n_w)                  :: j_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w)   :: dw_loc     
    INTEGER                                                  ::  m, l , i, mode, index, k
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,6)   :: Vs
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8)   :: ray
    REAL(KIND=8)   :: user_time, tps, dummy
    REAL(KIND=8), DIMENSION(3)                  :: temps

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me


    tps = user_time(dummy)
    DO i = 1, SIZE(list_mode)
       mode = list_mode(i)
       index = 0
       DO m = 1, me
          j_loc = jj(:,m)
          DO k = 1, 6
             Vs(:,k) = V_in(j_loc,k,i)
          END DO
          DO l = 1, l_G
             index = index + 1
             dw_loc = dw(:,:,l,m)

             !--------On calcul le rayon du point gauss
             ray = SUM(mesh%rr(1,j_loc)*ww(:,l))

             !-----------------vitesse sur les points de Gauss---------------------------
             W(index,1,i) = SUM(Vs(:,1)*ww(:,l))
             W(index,3,i) = SUM(Vs(:,3)*ww(:,l))
             W(index,5,i) = SUM(Vs(:,5)*ww(:,l))

             W(index,2,i) = SUM(Vs(:,2)*ww(:,l))
             W(index,4,i) = SUM(Vs(:,4)*ww(:,l))
             W(index,6,i) = SUM(Vs(:,6)*ww(:,l))

             !-----------------rotational sur les points de Gauss---------------------------
             !coeff sur les cosinus 
             RotV(index,1,i) = mode/ray*W(index,6,i) &
                  -SUM(Vs(:,3)*dw_loc(2,:))
             RotV(index,4,i) =          SUM(Vs(:,2)*dw_loc(2,:)) &
                  -SUM(Vs(:,6)*dw_loc(1,:))
             RotV(index,5,i) =    1/ray*W(index,3,i) &
                  +SUM(Vs(:,3)*dw_loc(1,:)) &
                  -mode/ray*W(index,2,i)
             !coeff sur les sinus       
             RotV(index,2,i) =-mode/ray*W(index,5,i) &
                  -SUM(Vs(:,4)*dw_loc(2,:))
             RotV(index,3,i) =         SUM(Vs(:,1)*dw_loc(2,:)) &
                  -SUM(Vs(:,5)*dw_loc(1,:))
             RotV(index,6,i) =    1/ray*W(index,4,i) &
                  +SUM(Vs(:,4)*dw_loc(1,:))&
                  +mode/ray*W(index,1,i)
          ENDDO
       ENDDO
    END DO
    tps = user_time(dummy) - tps
    !WRITE(*,*) ' Tps dans la grande boucle', tps
    tps = user_time(dummy)
    temps = 0
    CALL FFT_PAR_CROSS_PROD(RotV, W, V_out, temps) 
    tps = user_time(dummy) - tps
    !WRITE(*,*) ' Tps dans FFT_PAR_PROD_VECT', tps
    !write(*,*) ' Temps de Comm   ', temps(1)
    !write(*,*) ' Temps de Calc   ', temps(2)
    !write(*,*) ' Temps de Chan   ', temps(3)

  END SUBROUTINE smb_cross_prod_gauss_sft_par

END MODULE subroutine_ns
