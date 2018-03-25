MODULE boundary


  LOGICAL :: test_de_convergence
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:), PRIVATE  :: val_dir

CONTAINS

  SUBROUTINE prep_VKN(vv_mesh)
    USE dir_nodes 
    USE dyn_line
    USE def_type_mesh 
    USE chaine_caractere
    USE sub_plot
    IMPLICIT NONE
    TYPE(mesh_type)                                   :: vv_mesh        
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:)         :: rr, v_test 
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:)         :: val_sup, val_cote
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:)           :: r_sup, z_cote
    !facteur de correction pour la vitesse des pales   
    REAL(KIND=8)                                      :: xi, dr, dz, d1, d2, a, b, c, x, PI, eps, S, z
    REAL(KIND=8),              DIMENSION(3)           :: va, vb, vc
    REAL(KIND=8)                                      :: t_al, vn
    INTEGER                                           :: n_c, n_dir
    INTEGER                                           :: ms, ls, i, k, m, n, side_sup, side_inf, side_cote

    INTEGER,            DIMENSION(3)                  :: vv_j_size
    TYPE(dyn_int_line), DIMENSION(3)                  :: vv_js_D
    LOGICAL,      ALLOCATABLE, DIMENSION(:,:)         :: DirV ! Type of BC
    INTEGER,      ALLOCATABLE, DIMENSION(:)           :: list_dirichlet_sides,bord_periodic
    INTEGER :: n_bord, nb_dirichlet_sides

    ! We recompute the dirichlet nodes
    OPEN(UNIT = 21, FILE = 'data', FORM = 'formatted', STATUS = 'unknown')
    CALL read_until(21, 'data_periodic')
    READ  (21, *)  n_bord
    IF (n_bord>0) THEN
       ALLOCATE(bord_periodic(2*n_bord))
       DO n= 1, n_bord
          READ  (21, *)  bord_periodic(2*n-1), bord_periodic(2*n)
       END DO
    END IF

    CALL read_until(21, 'data_condlim_ns')
    READ(21,*) 
    READ(21,*)

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
    CLOSE(21)
    ALLOCATE(rr(2,vv_j_size(1)))
    rr = vv_mesh%rr(:,vv_js_D(1)%DIL)
    !------------------------------------------------------------------------------


    PI = ACOS(-1.d0)
    n_dir = size(rr,2) !nombre de points contraints
    write(*,*) ' n_dir =============', n_dir, vv_j_size(1)

    ALLOCATE(val_sup(3,51),val_cote(3,21))
    ALLOCATE(val_dir(3,n_dir))
    !val_sup = valeur sur le vord superieur, 51 points en ordre croissant de r=0,1
    !val_cote = valeur sur le bord lateral, 21 points en ordre croissant de z=-0.75,0.75

    !abcisse des points de la paroi sup du cylindre vke
    ALLOCATE(r_sup(51),z_cote(21))
    r_sup = 0.d0     
    dr = 1.d0/50.d0  !dr maillage vke
    dz = 1.8d0/24.d0
    DO i=1, 51
       !ordre croissant en r, comme dans val_sup 
       r_sup(i) = (i-1)*dr  
       !ordre croissant en z, comme dans val_cote
       IF (i<22) THEN 
          z_cote(i) = -0.75d0+(i-1)*dz
       ENDIF
    ENDDO
    write(*,*) 'Lecture CL'
    CALL  LECT_CL_VKN(val_sup, val_cote)

    !JLG: Je pense que les points 10 et 12 sont incoherents.
    !JLG: Je les redefinis ici 
    val_cote(3,10) = (val_cote(3,9)+ val_cote(3,11))/2
    val_cote(3,12) = (val_cote(3,11)+ val_cote(3,13))/2

    write(*,*) 'projection conditions limites'
    DO n=1, n_dir 
       IF (rr(1,n) < 1.d0) THEN    !points d'un des disque
          i = INT(rr(1,n)/dr)+1 !La division euclidienne est beaucoup plus simple (JLG).
          a = r_sup(i);   va = val_sup(1:3,i)
          b = r_sup(i+1); vb = val_sup(1:3,i+1)
          IF (i==1) THEN
             c = r_sup(i+2); vc = val_sup(1:3,i+2)
          ELSE
             c = r_sup(i-1); vc = val_sup(1:3,i-1)
          END IF
          x = rr(1,n)
          val_dir(1:3,n) = (b-x)*(c-x)/((b-a)*(c-a))* va &
               +(c-x)*(a-x)/((c-b)*(a-b))* vb &
               +(a-x)*(b-x)/((a-c)*(b-c))* vc !Second order interpolation
          IF (rr(2,n) == -0.75d0) THEN !disque inferieur
             val_dir(2:3,n) = -val_dir(2:3,n)
          ENDIF
       ELSE !point du bord lateral
          i = INT((rr(2,n)+0.75)/dz)+1 !La division euclidienne est beaucoup plus simple (JLG).
          a = z_cote(i);   va = val_cote(1:3,i)
          IF (i==SIZE(z_cote)) THEN
             b = z_cote(i-1); vb = val_cote(1:3,i-1)
             c = z_cote(i-2); vc = val_cote(1:3,i-2)
          ELSE IF (i==1) THEN
             b = z_cote(i+1); vb = val_cote(1:3,i+1)
             c = z_cote(i+2); vc = val_cote(1:3,i+2)
          ELSE
             b = z_cote(i-1); vb = val_cote(1:3,i-1)
             c = z_cote(i+1); vc = val_cote(1:3,i+1)
          END IF
          x = rr(2,n)
          val_dir(1:3,n) = (b-x)*(c-x)/((b-a)*(c-a))* va &
               +(c-x)*(a-x)/((c-b)*(a-b))* vb &
               +(a-x)*(b-x)/((a-c)*(b-c))* vc !Second order interpolation
       ENDIF
    ENDDO

    !VERIFICATION DES PROFILS DE VITESSE OBTENUS

    DO i = 1, 21
       WRITE(14,*) i, val_cote(3,i)
       WRITE(15,*) i, val_cote(2,i)
    END DO
    DO i = 1, 51
       WRITE(16,*) i, val_sup(3,i) !*i*dr
    END DO

    DEALLOCATE(val_sup,val_cote)
    DEALLOCATE(r_sup,z_cote)
100 FORMAT(5(e11.5,3x))

    ALLOCATE (v_test(2,vv_mesh%np))

    v_test(2,vv_js_D(3)%DIL) = val_dir(2,:)
    DO ms = 1,  vv_mesh%mes
       i = vv_mesh%jjs(1,ms)
       IF (vv_mesh%rr(1,i)>0.999d0) write(12,100) vv_mesh%rr(2,i),v_test(2,i) 
       i = vv_mesh%jjs(2,ms)
       IF (vv_mesh%rr(1,i)>0.999d0) write(12,100) vv_mesh%rr(2,i),v_test(2,i) 
    END DO

    v_test = 0
    v_test(1,vv_js_D(1)%DIL) = val_dir(1,:)
    v_test(2,vv_js_D(3)%DIL) = val_dir(3,:)
    DO ms = 1,  vv_mesh%mes
       i = vv_mesh%jjs(1,ms)
       IF (vv_mesh%rr(1,i)>0.999d0) write(13,100) vv_mesh%rr(2,i),v_test(2,i)
       i = vv_mesh%jjs(2,ms)
       IF (vv_mesh%rr(1,i)>0.999d0) write(13,100) vv_mesh%rr(2,i),v_test(2,i) 
    END DO

    s = 0
    x = 0
    DO ms =1 , vv_mesh%mes
       z = SUM(vv_mesh%rr(2,vv_mesh%jjs(:,ms))*vv_mesh%gauss%wws(:,1))
       IF (ABS(z)<0.7499999d0) CYCLE
       DO ls = 1, vv_mesh%gauss%l_Gs
          vn = SUM(vv_mesh%gauss%wws(:,ls)*v_test(1,vv_mesh%jjs(:,ms)))*vv_mesh%gauss%rnorms(1,ls,ms) &
               +SUM(vv_mesh%gauss%wws(:,ls)*v_test(2,vv_mesh%jjs(:,ms)))*vv_mesh%gauss%rnorms(2,ls,ms)
          s = s + vv_mesh%gauss%rjs(ls,ms)*vn*SUM(vv_mesh%rr(1,vv_mesh%jjs(:,ms))*vv_mesh%gauss%wws(:,ls)) 
          x = x + vv_mesh%gauss%rjs(ls,ms)*SUM(vv_mesh%rr(1,vv_mesh%jjs(:,ms))*vv_mesh%gauss%wws(:,ls))
       END DO
    END DO

    write(*,*) ' flux avant correction', s, x
    DO i = 1, n_dir
       IF ((rr(2,i) > 0.7499999d0)) THEN
          val_dir(3,i) = val_dir(3,i) - s/x
       ELSE IF ((rr(2,i) < -0.7499999d0)) THEN
          val_dir(3,i) = val_dir(3,i) + s/x
       END IF
    END DO

    v_test = 0
    v_test(1,vv_js_D(1)%DIL) = val_dir(1,:)
    v_test(2,vv_js_D(3)%DIL) = val_dir(3,:)
    CALL plot_arrow_label(vv_mesh%jj, vv_mesh%rr, v_test,'vect.plt')
    s = 0
    x = 0
    DO ms =1 , vv_mesh%mes
       z = SUM(vv_mesh%rr(2,vv_mesh%jjs(:,ms))*vv_mesh%gauss%wws(:,1))
       IF (ABS(z)<0.7499999d0) CYCLE
       DO ls = 1, vv_mesh%gauss%l_Gs
          vn = SUM(vv_mesh%gauss%wws(:,ls)*v_test(1,vv_mesh%jjs(:,ms)))*vv_mesh%gauss%rnorms(1,ls,ms) &
               +SUM(vv_mesh%gauss%wws(:,ls)*v_test(2,vv_mesh%jjs(:,ms)))*vv_mesh%gauss%rnorms(2,ls,ms)
          s = s + vv_mesh%gauss%rjs(ls,ms)*vn*SUM(vv_mesh%rr(1,vv_mesh%jjs(:,ms))*vv_mesh%gauss%wws(:,ls))
          x = x + vv_mesh%gauss%rjs(ls,ms)
       END DO
    END DO

    write(*,*) ' flux apres corection ', s
    DEALLOCATE(v_test, rr)
  END SUBROUTINE prep_VKN

!============================================================================================
!                       CONDITIONS LIMITES POUR NAVIER_STOKES
!============================================================================================
  SUBROUTINE init_up(mesh_f,mesh_c,time,dt,list_mode,un_m1,un,pn_m1,pn,phin_m1,phin)

    USE boundary_anal
    USE def_type_mesh

    IMPLICIT NONE
    TYPE(mesh_type)                            :: mesh_f, mesh_c        
    REAL(KIND=8),                   INTENT(OUT):: time
    REAL(KIND=8),                   INTENT(IN) :: dt
    INTEGER,      DIMENSION(:),     INTENT(IN) :: list_mode   
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT):: un_m1, un
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT):: pn_m1, pn, phin_m1, phin

    IF (test_de_convergence) THEN
       CALL init_up_analytique(mesh_f,mesh_c,time,dt,list_mode,un_m1,un,pn_m1,pn,phin_m1,phin)
       RETURN
    END IF

    un_m1   = 0.d0
    un      = 0.d0
    pn_m1   = 0.d0
    pn      = 0.d0
    phin_m1 = 0.d0
    phin    = 0.d0   
    time    = 0.d0

  END SUBROUTINE init_up

  FUNCTION ff_source(type,rr,m,t, Re, ty) RESULT(vv)
    !=================================== 
    USE boundary_anal

    IMPLICIT NONE

    INTEGER     ,                        INTENT(IN)   :: type
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv
    INTEGER     ,                        INTENT(IN)   :: m  !mode  
    REAL(KIND=8),                        INTENT(IN)   :: t   
    INTEGER                                           :: n
    CHARACTER(LEN=2),                    INTENT(IN)   :: ty
    REAL(KIND=8),                        INTENT(IN)   :: Re


    IF (test_de_convergence) THEN
       vv = ff_anal(type,rr,m,t,Re,ty)
       RETURN
    END IF

    vv = 0.d0
    RETURN

  END FUNCTION ff_source

  FUNCTION pp_exact(type,rr,m,t) RESULT (vv)
    
    USE boundary_anal

    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: type
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv 
    INTEGER     ,                        INTENT(IN)   :: m  !mode
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8)                                      :: PI, k
    
    IF (test_de_convergence) THEN
       vv = pp_exact_anal(type,rr,m,t)
       RETURN
    END IF

    vv=0.d0

   END FUNCTION pp_exact

  ! Extension of the velocity field in the solid
  FUNCTION extension_vel(type,rr,m,t) RESULT(vv)

    USE boundary_anal
    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: type
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER,                             INTENT(IN)   :: m  !mode 
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8)                                      :: Re 
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv, r, z 
    INTEGER :: n
    REAL(KIND=8)                                      :: PI, k, Ri

    vv = 0.d0
!   Ri = 1.d0
!   IF(m==0 .and. type == 3) THEN
!     DO n = 1, SIZE(rr,2)
!        vv(n) =  t**2/(0.1d0+t**2)*rr(1,n)/Ri !           1.d0
!     END DO
!   END IF
    RETURN

  END FUNCTION extension_vel 



!============================================================================================
!                       CONDITIONS LIMITES POUR MAXWELL
!============================================================================================

  FUNCTION Vexact(m, H_mesh) RESULT(vv)  !Set uniquement a l'induction

    USE boundary_anal
    USE def_type_mesh
    IMPLICIT NONE

    TYPE(mesh_type),                       INTENT(IN) :: H_mesh     !type de maillage
    INTEGER,                               INTENT(IN) :: m
    REAL(KIND=8), DIMENSION(H_mesh%np,6)              :: vv
    REAL(KIND=8), DIMENSION(:), POINTER               :: r, z
    REAL(KIND=8):: rho, phi,Ri
    INTEGER :: i, n
    REAL(Kind=8) :: pol_sur_tor,  PI, rr, zz, eps

    IF (test_de_convergence) THEN
       vv = Vexact_anal(m, H_mesh)
       RETURN
    END IF

  !Taylor couette avec cylindre interieur anime
    !Ri = 1.d0
    !vv = 0.d0
    !IF(m==0 .and. type == 3) THEN
    !vv(n) =  t**2/(0.1+t**2)*rr(1,n)/Ri !           1.d0   
    !END IF
    !END DO
    !END IF
    !RETURN


     vv = 0.d0
! VK experimental
     IF (m==0) THEN
       CALL init_flow_VKE(H_mesh, vv, 1)
     END IF

  END FUNCTION Vexact

  FUNCTION Hexact_init(type, rr, m, mu_H_field, t) RESULT(vv) 

    USE boundary_anal
    USE bessel
    IMPLICIT NONE

    INTEGER     ,                        INTENT(IN)   :: type
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m  !mode 
    REAL(KIND=8),                        INTENT(IN)   :: t 
    REAL(KIND=8), DIMENSION(:),          INTENT(IN)   :: mu_H_field
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv
    REAL(KIND=8)                                      :: PI, Lz, ri, ro, k0
    REAL(KIND=8),DIMENSION(SIZE(rr,2))                :: r, z, f, fp
    INTEGER                                           :: n
    
    IF (test_de_convergence) THEN
       vv = Hexact_anal(type, rr, m, mu_H_field, t)
       RETURN
    END IF

    vv(:) = 0.d0
!
!! CHAMP J. LEORAT pour TC non lineaire; COMPATIBLE avec phi=0 pr r>ro
!   PI = ACOS(-1.d0)
!   r = rr(1,:)
!   z = rr(2,:)
!
!    Lz = 2.d0*PI
!    ri = 1.d0                         !rayon du cylindre interieur
!    ro = 2.d0                         !rayon du cylindre exterieur
!    k0 = 2*2.d0*PI/Lz
!    f  = (ri-r)**2*(ro-r)**2
!    fp = -2.d0*((ri-r)**2*(ro-r)+(ri-r)*(ro-r)**2)
!    amp_h = 1d-1
!
!    IF (m == 1) THEN
!       IF (type == 1) THEN
!          vv = f/r*SIN(k0*z)
!       ELSEIF (type == 4) THEN
!          vv = -fp*SIN(k0*z)
!       ENDIF
!    ENDIF
!    vv = amp_h * vv
!
!! champ horizontal Bx (m=1) en extinction et champ vertical Bz (m=0) en extinction
! H_theta(m=0) peut etre qcq
!
    PI = ACOS(-1.d0)
    Lz=1.8d0
    ro=1.d0
    IF (m == 1) THEN
       IF (type == 1) THEN
          vv = 0.05d0          !chp Hx constant; VOIR avec Phiexact
       ELSEIF (type == 4) THEN
          vv = -0.05d0         !chp Hx constant; VOIR avec Phiexact
       ENDIF
    ELSE IF ((m==0).AND.(type==3)) THEN
      vv = 0.05d0*sin(PI*rr(1,:)/ro)
      !vv = 0.05d0*cos(PI*rr(2,:)/Lz)*sin(PI*rr(1,:)/ro)
    ELSE IF ((m==0).AND.(type==5)) THEN ! champ Hz en extinction; VOIR Phiexact
      vv = 0.05d0
    ENDIF

!! cas alpha EsR = champ J LEORAT entre e=ri <r< R=1 pr dependance en z
!!
!    PI = ACOS(-1.d0)
!    r = rr(1,:)
!    z = rr(2,:)
!
!    Lz=2.0d0
!    ro=1.d0
!    ri=0.5d0 ! a changer avec l'epaisseur pour le cas alpha EsR
!    k0 = PI
!    f  = (ri-r)**2*(ro-r)**2
!    fp = -2.d0*((ri-r)**2*(ro-r)+(ri-r)*(ro-r)**2)
!    print*,'Lz,ro,ri,k0=', Lz,ro,ri,k0
!
!  DO n=1, SIZE(rr,2)
!   IF ((ri.LE.r(n)).AND.(r(n).LE.ro)) THEN
!    ! print*,'r(n), z(n)= ',r(n), z(n)
!    IF (m == 1) THEN
!       IF (type == 1) THEN
!         !vv = 0.05d0 !chp Hx constant; VOIR avec Phiexact
!         !vv = 0.05d0*(1.d0 + f/r*SIN(k0*z))! chp Hx constant + chp variant avec z
!          vv(n) = 1.d0*(f(n)*SIN(k0*z(n))/r(n))
!       ELSEIF (type == 4) THEN
!         !vv = -0.05d0! chp Hx constant; VOIR avec Phiexact
!         !vv = -0.05d0*(1.d0 +fp*SIN(k0*z))! chp Hx constant + chp variant avec z
!          vv(n) = -1.d0*(fp(n)*SIN(k0*z(n)))
!       ENDIF
!    ELSE IF ((m==0).AND.(type==3)) THEN
!      vv(n) = 0.05d0*sin(PI*r(n)/ro)*SIN(k0*z(n))
!      !vv = 0.05d0*cos(PI*rr(2,:)/Lz)*sin(PI*rr(1,:)/ro)
!    ENDIF
!   ENDIF
!  ENDDO

 ! avant le 11/03/08
!   ELSE IF ((m==0).AND.(type==3)) THEN
!     vv = 0.05d0 ! pb de H_theta constant en r=0
!   ELSE IF ((m==0).AND.(type==5)) THEN ! champ Hz en extinction; VOIR Phiexact
!     vv = 0.05d0
!   ENDIF


    
  END FUNCTION Hexact_init

  FUNCTION Hexact(type, rr, m, mu_H_field, t) RESULT(vv) 

    USE boundary_anal
    USE bessel
    IMPLICIT NONE

    INTEGER     ,                        INTENT(IN)   :: type
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m  !mode 
    REAL(KIND=8),                        INTENT(IN)   :: t 
    REAL(KIND=8), DIMENSION(:),          INTENT(IN)   :: mu_H_field
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv
    REAL(KIND=8)                                      :: PI, Lz, ri, ro, k0
    REAL(KIND=8),DIMENSION(SIZE(rr,2))                :: r, z, f, fp, amp_h


    IF (test_de_convergence) THEN
       vv = Hexact_anal(type, rr, m, mu_H_field, t)
       RETURN
    END IF

    vv(:) = 0.d0
!
!CHAMP J. LEORAT
!   PI = ACOS(-1.d0)
!   r = rr(1,:)
!   z = rr(2,:)
!
!    Lz = 2.d0*PI
!    ri = 1.d0                         !rayon du cylindre interieur
!    ro = 2.d0                         !rayon du cylindre exterieur
!    k0 = 2*2.d0*PI/Lz
!    f  = (ri-r)**2*(ro-r)**2
!    fp = -2.d0*((ri-r)**2*(ro-r)+(ri-r)*(ro-r)**2)
!    amp_h = 1d-1
!
!    IF (m == 1) THEN
!       IF (type == 1) THEN
!          vv = f/r*SIN(k0*z)
!       ELSEIF (type == 4) THEN
!          vv = -fp*SIN(k0*z)
!       ENDIF
!    ENDIF
!    vv = amp_h * vv
!
!! champ horizontal Bx en extinction et champ vertical Bz constant
!
!    IF (m == 1) THEN
!       IF (type == 1) THEN
!          vv = 0.05d0
!       ELSEIF (type == 4) THEN
!          vv = -0.05d0
!       ENDIF
!    ELSE IF ((m==0).AND.(type==5)) THEN
!      vv = 0.05d0
!    ENDIF
!

    
  END FUNCTION Hexact

  FUNCTION Phiexact(type, rr, m, mu_phi,t) RESULT(vv) 
     
    USE boundary_anal
    USE bessel
    IMPLICIT NONE

    INTEGER     ,                        INTENT(IN)   :: type
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m  !mode 
    REAL(KIND=8),                        INTENT(IN)   :: mu_phi, t
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv   
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: r, z
    INTEGER                                           :: n
    REAL(KIND=8)                                      :: t_ext

    IF (test_de_convergence) THEN
       vv = Phiexact_anal(type, rr, m, mu_phi,t)
       RETURN
    END IF
 
    vv(:) = 0.d0
!
!! champ horizontal Bx (m=1) en extinction et champ vertical Bz (m=0) en extinction
! H_theta(m=0) peut etre qcq

   t_ext = .1d0 ! 1.d-2 
   
      IF ((m ==1).AND.(type==1)) THEN !champ horizontal Bx en extinction 
            vv = 0.05d0*rr(1,:)/(1.d0+(t/t_ext)**3)
      ELSE IF ((m ==0).AND.(type==1)) THEN !champ vertical Bz en extinction
            vv = 0.05d0*rr(2,:)/(1.d0+(t/t_ext)**3)
      ENDIF 



  END FUNCTION Phiexact


  FUNCTION Jexact_gauss(type, rr, m, mu_phi, sigma, mu_H, t) RESULT(vv) 


    USE boundary_anal
    USE bessel
    IMPLICIT NONE

    INTEGER     ,                        INTENT(IN)   :: type
    REAL(KIND=8), DIMENSION(:),          INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m  !mode 
    REAL(KIND=8),                        INTENT(IN)   :: mu_phi, sigma, mu_H, t 
    REAL(KIND=8)                                      :: vv
    REAL(KIND=8)                                      :: r, z
    REAL(KIND=8)                                      :: f1, f2, f3, f4
    REAL(KIND=8)                                      :: df2, df3, d2f2, df1, d2f3, df4
    REAL(KIND=8)                                      :: f1_r0, f2_r0, f3_r0
    REAL(KIND=8)                                      :: df2_r0, df3_r0
    REAL(KIND=8)                                      :: A, B, z0

    REAL(KIND=8)                                      :: H1, H2, H3, dH3, dH2

    IF (test_de_convergence) THEN
       vv = Jexact_gauss_anal(type, rr, m, mu_phi, sigma, mu_H, t) 
       RETURN
    END IF

    vv = 0.d0


  END FUNCTION Jexact_gauss

  FUNCTION Eexact_gauss(type, rr, m, mu_phi, sigma, mu_H, t) RESULT(vv) 

    USE boundary_anal
    USE bessel
    IMPLICIT NONE

    INTEGER,                             INTENT(IN)   :: type
    REAL(KIND=8), DIMENSION(:),          INTENT(IN)   :: rr
    INTEGER,                             INTENT(IN)   :: m  !mode 
    REAL(KIND=8),                        INTENT(IN)   :: mu_phi, sigma, mu_H, t
    REAL(KIND=8)                                      :: vv 
    REAL(KIND=8)                                      :: r, z
    REAL(KIND=8)                                      :: f1, f2, f3, f4
    REAL(KIND=8)                                      :: f1_r0, f2_r0, f3_r0
    REAL(KIND=8)                                      :: df2_r0, df3_r0
    REAL(KIND=8)                                      :: A, B, z0

    IF (test_de_convergence) THEN
       vv = Eexact_gauss_anal(type, rr, m, mu_phi, sigma, mu_H, t) 
       RETURN
    END IF

   vv = 0.d0
 

  END FUNCTION Eexact_gauss


  SUBROUTINE init_maxwell(H_mesh,phi_mesh,time,dt,mu_H_field,mu_phi,list_mode,Hn1,Hn,phin1,phin)

    USE boundary_anal
    USE def_type_mesh


    IMPLICIT NONE
    TYPE(mesh_type)                            :: H_mesh, phi_mesh     
    REAL(KIND=8),                   INTENT(OUT):: time
    REAL(KIND=8),                   INTENT(IN) :: dt
    REAL(KIND=8),                   INTENT(IN) :: mu_phi
    REAL(KIND=8), DIMENSION(:),     INTENT(IN) :: mu_H_field
    INTEGER,      DIMENSION(:),     INTENT(IN) :: list_mode    
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT):: Hn, Hn1
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT):: phin, phin1
    INTEGER                                    :: i, k, mode, m_max_c
    REAL(KIND=8), DIMENSION(:),   ALLOCATABLE  :: x
 
    ALLOCATE(x(H_mesh%np))
    
    IF (test_de_convergence) THEN
       CALL init_maxwell_analytique(H_mesh,phi_mesh,time,dt,mu_H_field,mu_phi,list_mode,Hn1,Hn,phin1,phin) 
       RETURN
    END IF

    time = -dt
    m_max_c = size(list_mode)
    DO k=1,6
      DO i=1, m_max_c
         mode = list_mode(i)
         Hn1(:,k,i) = Hexact_init(k, H_mesh%rr, mode, mu_H_field, time)
!        Hn1(:,k,i) = Hexact(k, H_mesh%rr, mode, mu_H_field, time)
         IF (k<3) THEN
            phin1(:,k,i) = Phiexact(k, phi_mesh%rr, mode, mu_phi, time)
         ENDIF
      ENDDO
    ENDDO

    time = time + dt
    DO k=1,6
      DO i=1, m_max_c
         mode = list_mode(i)
         Hn(:,k,i) = Hexact_init(k, H_mesh%rr, mode, mu_H_field, time)
!        Hn(:,k,i) = Hexact(k, H_mesh%rr, mode, mu_H_field, time)
         IF (k<3) THEN
            phin(:,k,i) = Phiexact(k, phi_mesh%rr, mode, mu_phi, time)
         ENDIF
      ENDDO
    ENDDO

   DO i = 1, SIZE(list_mode)
      IF (list_mode(i)==0) THEN
         CALL RANDOM_SEED
         CALL RANDOM_NUMBER(x)
         Hn(:,3,i)= Hn(:,3,i)*(x - 0.5d0)
      END IF
   END DO


  END SUBROUTINE init_maxwell


FUNCTION vv_exact(type, rr, mode, t) RESULT(vv) 
  
   USE def_type_mesh
   USE sub_plot
   IMPLICIT NONE
   
    INTEGER     ,                        INTENT(IN)   :: type
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER,                             INTENT(IN)   :: mode  !mode
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8), DIMENSION(size(rr,2))    :: vv
    REAL(KIND=8) :: t0

    t0 = .1d0
    vv = 0.d0
    IF (mode == 0) THEN
       IF (type == 1) THEN
          vv(:) = (t/t0)**2/(1.+(t/t0)**2)*val_dir(1,:)
       ELSEIF (type == 3) THEN
          vv(:) = (t/t0)**2/(1.+(t/t0)**2)*val_dir(2,:)
       ELSEIF (type == 5) THEN
          vv(:) = (t/t0)**2/(1.+(t/t0)**2)*val_dir(3,:)
       ENDIF
    ENDIF

END FUNCTION vv_exact


 SUBROUTINE LECT_CL_VKN(val_sup, val_cote)
!initialisation du champ de vitesse a partir du champ
!mesure dans VKEau.
!la lecture des donnees se fait ici ainsi que la projection sur le 
!maillage FEM
!num_dom est le numero du domain contenant le fluide em mvt
!z \in [-0.9,0.9]
!r \in [0,1]

 USE gauss_points
 USE def_type_mesh
 USE chaine_caractere

 IMPLICIT NONE

   TYPE(mesh_type)                                     :: mesh_vke
   !composantes 1 -> r , 2 -> theta, 3 -> z
   !valeurs des composantes de la vitesse
   !sur le disque superieur du champ vke sur le maillage expe
   !c'est pour cette raison qu'il y a 51 points
   REAL(KIND=8), DIMENSION(3,51)        , INTENT(OUT)   :: val_sup 
   !valeurs des composantes de la vitesse
   !sur la paroi laterale du champ vke sur le maillage expe, pour z\in[-0.7,0.7]
   !c'est pour cette raison qu'il y a 21 points
   REAL(KIND=8), DIMENSION(3,21)        , INTENT(OUT)   :: val_cote    
   !champ de vitesse lu :  r, th, z
   REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:)       :: v_vke
   INTEGER                                         :: i, j, k, m, n 
   REAL(KIND=8)                                    :: dr, dz, r, z, min, dist, eps
   REAL(KIND=8)                                    :: r_vke, z_vke, r_fem, z_fem
   REAL(KIND=8)                                    :: r_ref, z_ref
   INTEGER                                         :: num_ligne, me_ligne, n_min
   INTEGER                                         :: np_ligne, n_inf, n_sup
   INTEGER                                         :: n1, m1
   LOGICAL                                         :: test, mark
   REAL(KIND=8)                                    :: Vmax
!---------------END OF DECLARATION---------------------------------

   !initialisation du maillage de donnees (vke)
   !-------------------------------------
   write(*,*) 'maillage des donnees'
   dr = 1.d0/50.d0  !car r \in [0,1]
   dz = 1.8d0/24.d0 !car z \in [-0.9,0.9]

   mesh_vke%me = 1200      !=nel_r*nel_z=50*24
   mesh_vke%np = 1275      !=np_r*np_z  =51*25  
   mesh_vke%gauss%n_w = 4  !4 sommets par carre
   me_ligne = 50
   np_ligne = 51
   
   ALLOCATE(mesh_vke%jj(4,mesh_vke%me))
   ALLOCATE(mesh_vke%rr(2,mesh_vke%np))
   ALLOCATE(v_vke(3,mesh_vke%np))

   !tableau de connectivite
   !numerotation globale de gauche a droite et de haut en bas
   write(*,*) 'tableau de connectivite'
   num_ligne = 1
   DO m=1, mesh_vke%me
      IF (m > me_ligne*num_ligne) THEN
         num_ligne = num_ligne + 1
      ENDIF    
      mesh_vke%jj(1,m) = m + num_ligne - 1             !coin sup gauche
      mesh_vke%jj(2,m) = mesh_vke%jj(1,m) + 1          !coin sup droite
      mesh_vke%jj(3,m) = m + me_ligne + num_ligne + 1  !coin inf droit
      mesh_vke%jj(4,m) = mesh_vke%jj(3,m) - 1          !coin inf gauche
   ENDDO   
   !calcul des coordonnees
   write(*,*) 'tableau des coordonnees'
   num_ligne = 1
   eps=1.d-4

   DO i=1, mesh_vke%np
      IF (i == 1) THEN
          mesh_vke%rr(1,i) = 0.d0
          mesh_vke%rr(2,i) = -0.9d0
      ELSEIF(i > np_ligne*num_ligne) THEN  !alors changement de ligne et incre de z
         num_ligne = num_ligne + 1
         mesh_vke%rr(1,i) = 0.d0
         mesh_vke%rr(2,i) = mesh_vke%rr(2,i-1) + dz
      ELSE   !alors on incremente en r cae onreste sur la meme ligne
         mesh_vke%rr(1,i) = mesh_vke%rr(1,i-1) + dr
         mesh_vke%rr(2,i) = mesh_vke%rr(2,i-1)
      ENDIF
      IF ((abs(mesh_vke%rr(2,i)-0.9) .le. eps).OR. &
          (abs(mesh_vke%rr(2,i)-0.9) .le. 1.d-16)) THEN
! CN 29/09/08      IF (abs(mesh_vke%rr(2,i)-0.9) .le. eps) THEN
         mesh_vke%rr(2,i) = 0.9d0
      ENDIF   
   ENDDO

   !lecture du champ de vitesse expe
   !------------------------------------
   write(*,*) 'lecture des du champ VKEau'
    OPEN (UNIT = 31, FILE ='tm732.div0', & 
      FORM='formatted', STATUS='unknown')
    write(*,*) 'composante radiale'
    DO i=1, 25   !25 lignes
       n_inf = 1 + np_ligne*(i-1)
       n_sup = n_inf + (np_ligne-1)
       READ (31, *) v_vke(1,n_inf:n_sup) 
    ENDDO  
    write(*,*) 'composante azimuthale'
    DO i=1, 25  
       n_inf = 1 + np_ligne*(i-1)
       n_sup = n_inf + (np_ligne-1)
       READ (31, *) v_vke(2,n_inf:n_sup)
    ENDDO
    write(*,*) 'composante axiale'
    DO i=1, 25
       n_inf = 1 + np_ligne*(i-1)
       n_sup = n_inf + (np_ligne-1)
       READ (31, *) v_vke(3,n_inf:n_sup)
    ENDDO
    CLOSE(31)
    !normalisation du champ de vitesse par ||Vmax||
    Vmax = 0.d0
    DO i=1, mesh_vke%np
       !write(*,*) sqrt(v_vke(1,i)**2+v_vke(2,i)**2+v_vke(3,i)**2)
       IF (Vmax .lt. sqrt(v_vke(1,i)**2+v_vke(2,i)**2+v_vke(3,i)**2)) THEN
          Vmax = sqrt(v_vke(1,i)**2+v_vke(2,i)**2+v_vke(3,i)**2)
       ENDIF
    ENDDO   
    v_vke = v_vke / Vmax  
    !avant avant derniere ligne (z=0.75) dans val_sup(z=0.75) et val_(z=-0,75) 
      DO i=1, 51
         val_sup(1,i) = v_vke(1,mesh_vke%np-153+i) 
         val_sup(2,i) = v_vke(2,mesh_vke%np-153+i)
         val_sup(3,i) = v_vke(3,mesh_vke%np-153+i)
         !write(10,*) val_sup(1,i), val_sup(2,i),val_sup(3,i)
      ENDDO   
      DO i=1, 21
         !on saute les deux premieres lignes
         val_cote(1,i) = v_vke(1,np_ligne*(i+2))
         val_cote(2,i) = v_vke(2,np_ligne*(i+2))
         val_cote(3,i) = v_vke(3,np_ligne*(i+2))
         !write(11,*) val_cote(1,i), val_cote(2,i),val_cote(3,i)
      ENDDO
  
END SUBROUTINE  LECT_CL_VKN

SUBROUTINE init_flow_VKE(mesh_fem, v_fem, num_dom)
!initialisation du champ de vitesse a partir du champ
!mesure dans VKEau.
!la lecture des donnees se fait ici ainsi que la projection sur le 
!maillage FEM
!num_dom est le numero du domain contenant le fluide em mvt
!z \in [-0.9,0.9]
!r \in [0,1]

 USE gauss_points
 USE def_type_mesh
 USE chaine_caractere

 IMPLICIT NONE

   TYPE(mesh_type), TARGET                             :: mesh_fem
   !champ de vitesse pour les trois composantes pour le mode 0
   !les comp 2,4 et 6(correspondant a des sinus sont nulles car m=0)
   !champ de vitesse VKEau axisymetrie
   TYPE(mesh_type)                                     :: mesh_vke
   REAL(KIND=8), DIMENSION(mesh_fem%np,6), INTENT(OUT) :: v_fem
   INTEGER,                                INTENT(IN)  :: num_dom
   !champ de vitesse lu, r, th, z
   REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:)       :: v_vke
   INTEGER                                         :: i, j, k, m, n 
   REAL(KIND=8)                                    :: dr, dz, r, z, min, dist
   REAL(KIND=8)                                    :: r_vke, z_vke, r_fem, z_fem
   REAL(KIND=8)                                    :: r_ref, z_ref
   INTEGER                                         :: num_ligne, me_ligne, n_min
   INTEGER                                         :: np_ligne, n_inf, n_sup
   REAL(KIND=8)                                    :: T1, T2, T3, T4, eps
   INTEGER                                         :: j1, j2, j3, j4, m_loc
   INTEGER                                         :: n1, m1
   LOGICAL                                         :: test, mark
   REAL(KIND=8)                                    :: Vmax
!---------------END OF DECLARATION---------------------------------

   v_fem = 0.d0

   !initialisation du maillage de donnees (vke)
   !-------------------------------------
   write(*,*) 'maillage des donnees'
   dr = 1.d0/50.d0  !car r \in [0,1]
   dz = 1.8d0/24.d0 !car z \in [-0.9,0.9]

   mesh_vke%me = 1200      !=nel_r*nel_z=50*24
   mesh_vke%np = 1275      !=np_r*np_z  =51*25  
   mesh_vke%gauss%n_w = 4  !4 sommets par carre
   me_ligne = 50
   np_ligne = 51
   
   ALLOCATE(mesh_vke%jj(4,mesh_vke%me))
   ALLOCATE(mesh_vke%rr(2,mesh_vke%np))
   ALLOCATE(v_vke(3,mesh_vke%np))

   !tableau de connectivite
   !numerotation globale de gauche a droite et de haut en bas
   write(*,*) 'tableau de connectivite'
   num_ligne = 1
   DO m=1, mesh_vke%me
      IF (m > me_ligne*num_ligne) THEN
         num_ligne = num_ligne + 1
      ENDIF    
      mesh_vke%jj(1,m) = m + num_ligne - 1             !coin sup gauche
      mesh_vke%jj(2,m) = mesh_vke%jj(1,m) + 1          !coin sup droite
      mesh_vke%jj(3,m) = m + me_ligne + num_ligne + 1  !coin inf droit
      mesh_vke%jj(4,m) = mesh_vke%jj(3,m) - 1          !coin inf gauche
   ENDDO   
   !calcul des coordonnees
   write(*,*) 'tableau des coordonnees'
   num_ligne = 1
   eps=1.d-4
   DO i=1, mesh_vke%np
      IF (i == 1) THEN
          mesh_vke%rr(1,i) = 0.d0
          mesh_vke%rr(2,i) = -0.9d0
      ELSEIF(i > np_ligne*num_ligne) THEN  !alors changement de ligne et incre de z
         num_ligne = num_ligne + 1
         mesh_vke%rr(1,i) = 0.d0
         mesh_vke%rr(2,i) = mesh_vke%rr(2,i-1) + dz
      ELSE   !alors on incremente en r cae onreste sur la meme ligne
         mesh_vke%rr(1,i) = mesh_vke%rr(1,i-1) + dr
         mesh_vke%rr(2,i) = mesh_vke%rr(2,i-1)
      ENDIF
      IF (abs(mesh_vke%rr(2,i)-0.9) .le. eps) THEN
         mesh_vke%rr(2,i) = 0.9d0
      ENDIF   
   ENDDO

!TEST
   ! DO i=1, mesh_vke%np
   !    write(70,*) i ,mesh_vke%rr(1,i),mesh_vke%rr(2,i)
   ! ENDDO
   ! DO i=1, mesh_vke%me
   !    write(80,*) i, mesh_vke%jj(1,i), mesh_vke%jj(2,i), mesh_vke%jj(3,i), mesh_vke%jj(4,i)
   ! ENDDO
   ! STOP
!TEST

   !lecture du champ de vitesse expe
   !------------------------------------
   write(*,*) 'lecture du champ VKEau'
    OPEN (UNIT = 31, FILE ='tm732.div0', & 
      FORM='formatted', STATUS='unknown')
    write(*,*) 'composante radiale'
    DO i=1, 25   !25 lignes
       n_inf = 1 + np_ligne*(i-1)
       n_sup = n_inf + (np_ligne-1)
       READ (31, *) v_vke(1,n_inf:n_sup) 
    ENDDO  
    write(*,*) 'composante azimuthale'
    DO i=1, 25  
       n_inf = 1 + np_ligne*(i-1)
       n_sup = n_inf + (np_ligne-1)
       READ (31, *) v_vke(2,n_inf:n_sup)
    ENDDO
    write(*,*) 'composante axiale'
    DO i=1, 25
       n_inf = 1 + np_ligne*(i-1)
       n_sup = n_inf + (np_ligne-1)
       READ (31, *) v_vke(3,n_inf:n_sup)
    ENDDO
    !normalisation du champ de vitesse par ||Vmax||
    Vmax = 0.d0
    DO i=1, mesh_vke%np
       !write(*,*) sqrt(v_vke(1,i)**2+v_vke(2,i)**2+v_vke(3,i)**2)
       IF (Vmax .lt. sqrt(v_vke(1,i)**2+v_vke(2,i)**2+v_vke(3,i)**2)) THEN
          Vmax = sqrt(v_vke(1,i)**2+v_vke(2,i)**2+v_vke(3,i)**2)
       ENDIF
    ENDDO   
    v_vke = v_vke / Vmax  
    write(*,*) 'Vmax(champ vitesse)=',Vmax  
    !STOP  

!TEST genere fichier visualisable avec gnuplot
!   DO i=1, mesh_vke%np
       !write(80,*)  mesh_vke%rr(1,i), mesh_vke%rr(2,i), v_vke(1,i)
       !write(81,*)  mesh_vke%rr(1,i), mesh_vke%rr(2,i), v_vke(2,i)
       !write(82,*)  mesh_vke%rr(1,i), mesh_vke%rr(2,i), v_vke(3,i)
!   ENDDO
!   TEST
      
    !projection du champ sur le maillage FEM
    !-----------------------------------------
    write(*,*) 'projection du champ sur le maillage FEM'
    
    DO i=1, mesh_fem%np
       r_fem = mesh_fem%rr(1,i)
       z_fem = mesh_fem%rr(2,i)
       min = sqrt(dr**2+dz**2) +eps
       !write(*,*) 'i=',i,'r_fem=',r_fem,'z_fem=',z_fem,'min=',min
       !recherche d l'element contenant le point i
       test = .false.
       DO m=1, mesh_vke%me
          !write(*,*) '**********m=',m
          j1 = mesh_vke%jj(1,m) ; j2 = mesh_vke%jj(2,m)
          j3 = mesh_vke%jj(3,m) ; j4 = mesh_vke%jj(4,m)
          IF ((r_fem .ge. mesh_vke%rr(1,j1)) .and.(r_fem .le. mesh_vke%rr(1,j2)) &
          .and. (z_fem .ge. mesh_vke%rr(2,j1)) &
          .and. (z_fem .le. mesh_vke%rr(2,j3))) THEN
             test = .true.                   
             m_loc = m
             EXIT
          ENDIF
       ENDDO

       IF (test == .false.) THEN
          !recherche de l'element de mesh_fem d'appartenance du point
          mark = .false.
          DO m1 = 1, mesh_fem%me
             DO n1 = 1, mesh_fem%gauss%n_w
                IF ((i == mesh_fem%jj(n1,m1)) .and. (mesh_fem%i_d(m1) /= num_dom)) THEN
                    mark = .true.
                    EXIT
                ENDIF    
             ENDDO
             IF (mark) THEN
                EXIT
             ENDIF   
          ENDDO
          IF (.not.(mark)) THEN
             write(*,*) 'aucun element de &
             mesh_vke trouve contenant le point du maillage FEM : '
             write(*,*) 'i=',i,'r=',r_fem,'z=',z_fem
             STOP
          ENDIF    
       
       ELSE

!TEST POUR SAVOIR SI L'element contient le point!
     
      !write(*,*) 'i=',i,'r_fem=',r_fem,'z_fem=',z_fem
      !write(*,*) 'element vke = ',m_loc
      !write(*,*) 'rj1=',mesh_vke%rr(1,mesh_vke%jj(1,m_loc))
      !write(*,*) 'zj1=',mesh_vke%rr(2,mesh_vke%jj(1,m_loc))
      !write(*,*) 'rj2=',mesh_vke%rr(1,mesh_vke%jj(2,m_loc))
      !write(*,*) 'zj2=',mesh_vke%rr(2,mesh_vke%jj(2,m_loc))
      !write(*,*) 'rj3=',mesh_vke%rr(1,mesh_vke%jj(3,m_loc))
      !write(*,*) 'zj3=',mesh_vke%rr(2,mesh_vke%jj(3,m_loc))
      !write(*,*) 'rj4=',mesh_vke%rr(1,mesh_vke%jj(4,m_loc))
      !write(*,*) 'zj4=',mesh_vke%rr(2,mesh_vke%jj(4,m_loc))
!TEST

       !calcul des valeurs des fonctions de forme T_k au point du maillage FEM
       !on se ramene au carre de reference
          j1 = mesh_vke%jj(1,m_loc) ; j2 = mesh_vke%jj(2,m_loc)
          j3 = mesh_vke%jj(3,m_loc) ; j4 = mesh_vke%jj(4,m_loc)
          r_ref = 1.d0/dr * (r_fem - mesh_vke%rr(1,j1))
          z_ref = 1.d0/dz * (z_fem - mesh_vke%rr(2,j1))
          !write(*,*) 'r_fem=',r_fem,'mesh_vke%rr(1,j1)=',mesh_vke%rr(1,j1)
          !write(*,*) 'element ',m_loc,'r_ref=',r_ref,'z_ref=',z_ref
          T2 = r_ref        * (1.d0-z_ref)
          T3 = r_ref        * z_ref      
          T4 = (1.d0-r_ref) * z_ref     
          T1 = (1.d0-r_ref) * (1.d0-z_ref)
       !interpolation lineaire par element Q1
       !vitesse = SUM_j[v_vke(j)*T_j]
       !composante radiale
          v_fem(i,1) = T1*v_vke(1,j1) + T2*v_vke(1,j2) + T3*v_vke(1,j3) &
          + T4*v_vke(1,j4)
       !composante azimuthale
          v_fem(i,3) = T1*v_vke(2,j1) + T2*v_vke(2,j2) + T3*v_vke(2,j3) &
          + T4*v_vke(2,j4)
       !composante axiale
          v_fem(i,5) = T1*v_vke(3,j1) + T2*v_vke(3,j2) + T3*v_vke(3,j3) &
          + T4*v_vke(3,j4)
       ENDIF
    ENDDO
    
!TEST : genere des fichiers gnuplot pour le champ de vitesse FEM
   DO i=1, mesh_fem%np
      write(90,*)  mesh_fem%rr(1,i), mesh_fem%rr(2,i), v_fem(i,1)
      write(91,*)  mesh_fem%rr(1,i), mesh_fem%rr(2,i), v_fem(i,3)
      write(92,*)  mesh_fem%rr(1,i), mesh_fem%rr(2,i), v_fem(i,5)
   ENDDO

!   CALL comp_champ_p1_2d_vtk(mesh_fem, v_fem(1,:), 10, 'V_fem_r.vtk')
!   CALL comp_champ_p1_2d_vtk(mesh_fem, v_fem(3,:), 11, 'Vth_VKE.vtk')
!   CALL comp_champ_p1_2d_vtk(mesh_fem, v_fem(5,:), 12, 'V_fem_z.vtk')
!   CALL champ_vect_p1_2D_vtk(mesh_fem, v_fem(1,:), v_fem(5,:) &
!         , 13, 'V', 'Vmer_VKE.vtk')
!!
   !CALL comp_champ_p1_2d_vtk(mesh_vke, v_vke(1,:), 13, 'V_vke_r.vtk')
   !CALL comp_champ_p1_2d_vtk(mesh_vke, v_vke(2,:), 14, 'V_vke_th.vtk')
   !CALL comp_champ_p1_2d_vtk(mesh_vke, v_vke(3,:), 15, 'V_vke_z.vtk')   
!TEST
    
    !STOP

END SUBROUTINE init_flow_VKE

  FUNCTION alpha(rr) result(vv)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(2),          INTENT(IN)   :: rr
    REAL(KIND=8)                                      :: vv
    REAL(KIND=8)                                      :: z_top,z_bot,z_del,z,z_mil
    REAL(KIND=8)                                      :: r, PI
! Forcage alpha localise
    z_top =  0.8d0 !0.7d0 epais!0.8d0 alpha jusqu au 09/04/08
    z_bot = -0.8d0 !-0.7d0 epais !-0.8d0 alpha jusqu au 09/04/08
    z_del =  0.05d0
    z_mil =  0.1d0
    PI    =  ACOS(-1.d0)

         r     = rr(1)
         z     = rr(2)
         IF (r.le.1.d0) THEN
         ! 09/04/08 simple alpha avec z_top=0.8d0=z_bot OU 0.7 pour epais
          vv    = (2.d0 + tanh((z-z_top)/z_del) - tanh((z-z_bot)/z_del))/2.d0
         ! 10/04/08 double alpha avec z_top=0.8d0=z_bot,z_del=0.05d0,z_mil=0.1d0
         !vv    = (2.d0 + tanh((z-z_top)/z_del) - tanh((z-z_bot)/z_del))/2.d0 + &
         !        (tanh((z-z_mil)/z_del) - tanh((z+z_mil)/z_del))/2.d0
         ! 09/04/08 simple alpha avec z_top=0.8d0=z_bot OU 0.7 pour epais
         !vv  = (2.d0 + tanh((z-z_top)/z_del) - tanh((z-z_bot)/z_del))/2.d0
         !vv  = vv*r**2*SIN(PI*r)
         ! 24/09/08 Nice cas alpha constante par morceaux alpha_c pr
         ! comparer avec Andre Giesecke
         !IF((z.GE.z_top).OR.(z.LE.z_bot)) THEN
         !   vv = 1.d0
         !END IF
         ENDIF

    RETURN

  END FUNCTION alpha

  FUNCTION alpha_K(rr) result(vv)
    IMPLICIT NONE
    REAL(KIND=8)                                      :: vv
    REAL(KIND=8), DIMENSION(2),          INTENT(IN)   :: rr
    REAL(KIND=8)                                      :: r
! Forcage alpha_K Karlsruhe alpha=constante sur H_r(1:2) et H_theta(3:4)

         r     = rr(1)
         IF (r.le.1.d0) THEN
            vv    = 1.d0
         ENDIF

    RETURN

  END FUNCTION alpha_K

  FUNCTION alpha_EsR(rr) result(vv)
    IMPLICIT NONE
    REAL(KIND=8)                                      :: vv
    REAL(KIND=8), DIMENSION(2),          INTENT(IN)   :: rr
    REAL(KIND=8)                                      :: e, r
!
! Forcage alpha EsR alpha=constante sur H_r(1:2) et H_theta(3:4) entre r=e et R=1
!

         e = 0.5d0
         r     = rr(1)
         IF ((e.LE.r).AND.(r.LE.1.d0)) THEN
            vv    = 1.d0
         ENDIF
    RETURN

  END FUNCTION alpha_EsR

END MODULE boundary
