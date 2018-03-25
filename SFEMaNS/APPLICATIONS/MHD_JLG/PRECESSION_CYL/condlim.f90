MODULE boundary


  LOGICAL :: test_de_convergence
  !LOGICAL :: second_order_ext_pressure

CONTAINS

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
    INTEGER                                    :: i
    REAL(KIND=8), DIMENSION(:),   ALLOCATABLE  :: x
    REAL(KIND=8)                               :: ampl

    ALLOCATE(x(mesh_f%np))


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


! cas PRECESSION 28/07/09
    !initilisation rotation rigide
    DO i=1, SIZE(list_mode)
       IF (list_mode(i) == 0) THEN
          un_m1(:,3,i) = mesh_f%rr(1,:)
          un(:,3,i) = mesh_f%rr(1,:)
          pn_m1(:,1,i) = (mesh_c%rr(1,:))**2/2.d0
          pn(:,1,i) =  (mesh_c%rr(1,:))**2/2.d0
       ENDIF
    ENDDO

!    !initialisation onde inertielle
!    !CALL  init_up_oi(mesh_f,mesh_c,time,dt,list_mode,un_m1,un,pn_m1,pn,phin_m1,phin)
!
     ampl=1.d-4
    DO i = 1, SIZE(list_mode)
       IF ((list_mode(i).LE.5).AND.(list_mode(i).GE.1)) THEN
          CALL RANDOM_SEED
          CALL RANDOM_NUMBER(x)
          un(:,3,i)= ampl*(x - 0.5d0)
       END IF
    END DO
    RETURN

  END SUBROUTINE init_up


  SUBROUTINE init_up_oi(mesh_f,mesh_c,time,dt,list_mode,un_m1,un,pn_m1,pn,phin_m1,phin)

    USE boundary_anal
    USE def_type_mesh
    USE bessel

    IMPLICIT NONE
    TYPE(mesh_type)                            :: mesh_f, mesh_c
    REAL(KIND=8),                   INTENT(OUT):: time
    REAL(KIND=8),                   INTENT(IN) :: dt
    INTEGER,      DIMENSION(:),     INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT):: un_m1, un
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT):: pn_m1, pn, phin_m1, phin
    REAL(KIND=8), DIMENSION(size(mesh_f%rr,2))   :: rf, zf
    REAL(KIND=8), DIMENSION(size(mesh_c%rr,2))   :: rc, zc
    INTEGER :: i

    REAL(KIND=8)  :: kn, lambda, pi, d, amp 
    
    un_m1   = 0.d0
    un      = 0.d0
    pn_m1   = 0.d0
    pn      = 0.d0
    phin_m1 = 0.d0
    phin    = 0.d0
    time    = 0.d0
  
    !cas KB94 pour m=0, l=1 et n=1
    rf(:) = mesh_f%rr(1,:)
    rc(:) = mesh_c%rr(1,:)
    zf(:) = mesh_f%rr(2,:)
    zc(:) = mesh_c%rr(2,:)

    kn = 3.83171  !premiere racine de J1
    pi = acos(-1.d0)
    d  = 1.35045  !hauteur du cylindrepour KB94
    lambda = 2.d0/sqrt(1.d0+(kn**2*d**2)/(pi**2))  
    
    !initialisation vitesse
    DO i=1, mesh_f%np
        un_m1(i,1,1) = 2.d0*lambda*BESSJ1(kn*rf(i))*COS(pi*zf(i)/d)*SIN(lambda)
        un_m1(i,3,1) = 4.d0*BESSJ1(kn*rf(i))*COS(pi*zf(i)/d)*COS(lambda)
        un_m1(i,5,1) = -2.d0*lambda*d/pi*BESSJ0(kn*rf(i))*SIN(pi*zf(i)/d)*SIN(lambda)
        un(i,1,1) = 2.d0*lambda*BESSJ1(kn*rf(i))*COS(pi*zf(i)/d)*SIN(lambda*(1.d0+dt))
        un(i,3,1) = 4.d0*BESSJ1(kn*rf(i))*COS(pi*zf(i)/d)*COS(lambda*(1.d0+dt))
        un(i,5,1) = -2.d0*lambda*d/pi*BESSJ0(kn*rf(i))*SIN(pi*zf(i)/d)*SIN(lambda*(1.d0+dt))
    ENDDO
    !initialisation pression
    DO i=1, mesh_c%np
        pn_m1(i,1,1) = -1.d0/kn*BESSJ0(kn*rc(i))*COS(pi*zc(i)/d)*COS(lambda)
        pn(i,1,1) = -1.d0/kn*BESSJ0(kn*rc(i))*COS(pi*zc(i)/d)*COS(lambda*(1.d0+dt))
    ENDDO

    amp = 1.d-2
    pn = amp * pn
    pn_m1 = amp * pn_m1
    un = amp * un
    un_m1 = amp * un_m1
    

  END SUBROUTINE init_up_oi


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

  FUNCTION vv_exact(type,rr,m,t) RESULT(vv)
    
    USE boundary_anal

    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: type
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER,                             INTENT(IN)   :: m  !mode 
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8)                                      :: Re 
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv, r, z 
    INTEGER                                           :: n
    REAL(KIND=8)                                      :: PI, k, ri, ro, eps
    REAL(KIND=8), SAVE                                :: h
    LOGICAL, SAVE                                     :: once=.true.

    IF (test_de_convergence) THEN
       vv = vv_exact_anal(type,rr,m,t) 
       RETURN
    END IF

    vv = 0.d0
! cas PRECESSION 28/07/09
    IF ((m == 0) .and. (type == 3)) THEN
       vv(:) = rr(1,:)
    ENDIF
    RETURN

!  !Taylor couette avec cylindre interieur anime
!    Ri = 1.d0
!    Ro = 2.d0
!    eps = 1.d-4
!     !determination de h=LZ/2
!    IF (once) THEN
!       once = .false.
!       h = 0.d0
!       DO n=1, SIZE(rr,2)
!          IF (abs(rr(2,n)) > h) THEN
!              h=abs(rr(2,n))
!          ENDIF
!       ENDDO
!     write(*,*) 'h=',h
!    ENDIF
!
!    vv = 0.d0
!    IF(m==0 .and. type == 3) THEN
!      DO n = 1, SIZE(rr,2)
!         IF(rr(1,n) > Ro - eps ) THEN
!           vv(n) = 0
!         ELSE IF(rr(1,n) < Ri + eps) THEN
!           vv(n) = t**2/(0.1+t**2) !           1.d0
!         ELSEIF (abs(rr(2,n))>h-eps) THEN
!!        ELSEIF ((rr(2,n)>h-eps) .or. (rr(2,n)<-1.d0*eps+h)) THEN
!!          write(*,*) 'oki'
!           vv(n) = rr(1,n)* t**2/(0.1+t**2) !           1.d0
!         END IF
!      END DO
!    END IF
!    RETURN
!
!    !Von Karmann
!    vv = 0.d0
!    IF(m==0 .and. type == 3) THEN
!      DO n = 1, SIZE(rr,2)
!         IF(rr(1,n) > 0.49999d0) THEN
!           vv(n) = 0
!         ELSE IF(rr(2,n) < 0.001d0) THEN
!           vv(n) = rr(1,n)
!         ELSE
!           vv(n) = -rr(1,n)
!         END IF
!      END DO
!    END IF

  END FUNCTION vv_exact

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
! cas PRECESSION 28/07/09
    IF ((m == 0) .and. (type == 1)) THEN
       vv(:) = (rr(1,:))**2/2.d0
    ENDIF
    RETURN
   END FUNCTION pp_exact
!
  FUNCTION vv_exact_rot_axez(type,rr,m,t) RESULT(vv)
    
    USE boundary_anal

    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: type
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER,                             INTENT(IN)   :: m  !mode 
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8)                                      :: Re 
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv, r, z 
    INTEGER                                           :: n
    REAL(KIND=8)                                      :: PI, k, ri, ro, eps
    REAL(KIND=8), SAVE                                :: h
    LOGICAL, SAVE                                     :: once=.true.

    IF (test_de_convergence) THEN
       vv = vv_exact_anal(type,rr,m,t) 
       RETURN
    END IF

    vv = 0.d0
! cas PRECESSION 28/07/09
    IF ((m == 0) .and. (type == 3)) THEN
       vv(:) = rr(1,:)
    ENDIF
    RETURN
  END FUNCTION vv_exact_rot_axez

  ! Extension of the velocity field in the solid
  FUNCTION extension_vel(type, H_mesh, mode, t, n_start) RESULT(vv)
    USE def_type_mesh
    USE boundary_anal
    IMPLICIT NONE
    TYPE(mesh_type),                     INTENT(IN)   :: H_mesh
    INTEGER     ,                        INTENT(IN)   :: type
    INTEGER,                             INTENT(IN)   :: mode, n_start
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8), DIMENSION(SIZE(H_mesh%rr,2))        :: vv
    INTEGER                                           :: n
    REAL(KIND=8)                                      :: Ri

    Ri = 1.d0
    vv = 0.d0
!    IF(mode==0 .and. type == 3) THEN
!      DO n = n_start, SIZE(H_mesh%rr,2)
!         vv(n-n_start+1) = H_mesh%rr(1,n) ! t**2/(0.1+t**2)*H_mesh%rr(1,n)/Ri 
!      END DO
!    END IF
!    RETURN

  END FUNCTION extension_vel 



!============================================================================================
!                       CONDITIONS LIMITES POUR MAXWELL
!============================================================================================

  FUNCTION Vexact(m, H_mesh, time) RESULT(vv)  !Set uniquement a l'induction

    USE boundary_anal
    USE def_type_mesh
    USE chaine_caractere
    IMPLICIT NONE

    TYPE(mesh_type),                       INTENT(IN) :: H_mesh     !type de maillage
    INTEGER,                               INTENT(IN) :: m
! cas PRECESSION 28/07/09
    REAL(KIND=8), OPTIONAL,                INTENT(IN) :: time
    REAL(KIND=8), DIMENSION(H_mesh%np,6)              :: vv_int, vv

    REAL(KIND=8), DIMENSION(:), POINTER               :: r, z
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:)           :: r_i
    LOGICAL, SAVE                                     :: once = .true.
    REAL(KIND=8), ALLOCATABLE,DIMENSION(:,:,:), SAVE  :: vv_f
    INTEGER, SAVE                                     :: n, nr, mcut
    REAL(KIND=8), SAVE                                :: puls, amp
    REAL(KIND=8)                                      :: dth, mdth, x, Delta
    INTEGER                                           :: i, j
    IF (test_de_convergence) THEN
       vv = Vexact_anal(m, H_mesh)
       RETURN
    END IF

    vv = 0.d0
    RETURN

   !Cas Waleed - Jacques
   IF(once) THEN
      once = .false.
      nr = 21
      mcut = 16
      ALLOCATE(vv_f( 0:mcut-1,6,nr))
      !Lecture du champ de vitesse dans les fichiers
      !Composante radiale    
      OPEN(UNIT = 40, FILE = 'Crad.txt', FORM = 'formatted', STATUS = 'unknown')
WRITE(*,*) 'coucou1'
      READ(40,*)
      READ(40,*)
      DO i=1, 21
         READ(40,'(16(e13.7,2x))') vv_f(0:mcut-1,1,i)
      ENDDO  
      CLOSE(40)
      OPEN(UNIT = 40, FILE = 'Srad.txt', FORM = 'formatted', STATUS = 'unknown')
WRITE(*,*) 'coucou1'
      READ(40,*)
      READ(40,*)
      DO i=1, 21
         READ(40,'(16(e13.7,2x))') vv_f(0:mcut-1,2,i)
      ENDDO
      CLOSE(40)
      !composante azimutale
      OPEN(UNIT = 40, FILE = 'Caz.txt', FORM = 'formatted', STATUS = 'unknown')
WRITE(*,*) 'coucou1'
      READ(40,*)
      READ(40,*)
      DO i=1, 21
         READ(40,'(16(e13.7,2x))') vv_f(0:mcut-1,3,i)
      ENDDO
      CLOSE(40)
      OPEN(UNIT = 40, FILE = 'Saz.txt', FORM = 'formatted', STATUS = 'unknown')
WRITE(*,*) 'coucou1'
      READ(40,*)
      READ(40,*)
      DO i=1, 21
         READ(40,'(16(e13.7,2x))') vv_f(0:mcut-1,4,i)
      ENDDO
      CLOSE(40)

       OPEN(UNIT = 21, FILE = 'data', FORM = 'formatted', STATUS = 'unknown')
       CALL read_until(21, 'cyl_prec')
       READ (21,*) puls, amp
WRITE(*,*) 'puls, amp =', puls, amp
       CLOSE(21)
   ENDIF  ! fin du once

   !Interpolation   
   vv = 0.d0
   vv_int = 0.d0
   Delta = 1.d0/(nr-1) 
   ALLOCATE(r_i(nr))
   DO i=1, nr
      r_i(i) = (i-1.d0)*Delta
   ENDDO

   DO i=1, H_mesh%np
      x = H_mesh%rr(1,i)
      !recherche de sa position
      DO j=1, nr-1
         if ((r_i(j) .le. x).and.(x < r_i(j+1))) then
            exit
         endif
      ENDDO
      !composantes radiales
      vv(i,1) = vv_f(m,1,j)*(1.d0-(r_i(j)-x)/(r_i(j)-r_i(j+1))) &
               +vv_f(m,1,j+1)*(r_i(j)-x)/(r_i(j)-r_i(j+1))
      vv(i,2) = vv_f(m,2,j)*(1.d0-(r_i(j)-x)/(r_i(j)-r_i(j+1))) &
               +vv_f(m,2,j+1)*(r_i(j)-x)/(r_i(j)-r_i(j+1))
      !composantes azimutales
      vv(i,3) = vv_f(m,3,j)*(1.d0-(r_i(j)-x)/(r_i(j)-r_i(j+1))) &
               +vv_f(m,3,j+1)*(r_i(j)-x)/(r_i(j)-r_i(j+1))
      vv(i,4) = vv_f(m,4,j)*(1.d0-(r_i(j)-x)/(r_i(j)-r_i(j+1))) &
               +vv_f(m,4,j+1)*(r_i(j)-x)/(r_i(j)-r_i(j+1))
      !write(*,*) 'dans Vexact j=',j,vv(i,1), vv(i,2),vv(i,3),vv(i,4)
      

       IF (PRESENT(time)) THEN
         !write(*,*) 'calcul vitesse , time=', time
         dth  = puls*time
!        mdth = m*puls*time
         vv_int(i,1) = (vv(i,1)*cos(dth)-vv(i,3)*sin(dth))
         vv_int(i,2) = (vv(i,2)*cos(dth)-vv(i,4)*sin(dth))
         vv_int(i,3) = (vv(i,1)*sin(dth)+vv(i,3)*cos(dth))
         vv_int(i,4) = (vv(i,2)*sin(dth)+vv(i,4)*cos(dth))
       ENDIF

   ENDDO
       IF (PRESENT(time)) THEN
         vv = amp*vv_int
       ELSE
         vv = amp*vv
       ENDIF

   RETURN
  END FUNCTION Vexact

  FUNCTION Hexact_init(TYPE, rr, m, mu_H_field, t) RESULT(vv) 

    USE boundary_anal
    USE bessel
    IMPLICIT NONE

    INTEGER     ,                        INTENT(IN)   :: TYPE
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m  !mode 
    REAL(KIND=8),                        INTENT(IN)   :: t 
    REAL(KIND=8), DIMENSION(:),          INTENT(IN)   :: mu_H_field
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv
    REAL(KIND=8)                                      :: PI, Lz, ri, ro, k0
    REAL(KIND=8),DIMENSION(SIZE(rr,2))                :: r, z, f, fp, amp_h

    IF (test_de_convergence) THEN
       vv = Hexact_anal(TYPE, rr, m, mu_H_field, t)
       RETURN
    END IF

    vv(:) = 0.d0
     IF (m.ge.2) THEN
        vv = 1.d-5
     ENDIF
     RETURN
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
    !!
    !    PI = ACOS(-1.d0)
    !    ro=1.d0
    !    IF (m == 1) THEN
    !       IF (type == 1) THEN
    !          vv = 0.05d0
    !       ELSEIF (type == 4) THEN
    !          vv = -0.05d0
    !       ENDIF
    !    ELSE IF ((m==0).AND.(type==3)) THEN
    !      vv = 0.05d0*sin(PI*rr(1,:)/ro)
    !      !vv = 0.05d0*cos(PI*rr(2,:)/Lz)*sin(PI*rr(1,:)/ro)
    !    ENDIF
    ! ! avant le 11/03/08
    !!   ELSE IF ((m==0).AND.(type==3)) THEN
    !!     vv = 0.05d0
    !!   ELSE IF ((m==0).AND.(type==5)) THEN
    !!     vv = 0.05d0
    !!   ENDIF
    !


  END FUNCTION Hexact_init

  FUNCTION Hexact(TYPE, rr, m, mu_H_field, t) RESULT(vv) 

    USE boundary_anal
    USE bessel
    IMPLICIT NONE

    INTEGER     ,                        INTENT(IN)   :: TYPE
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m  !mode 
    REAL(KIND=8),                        INTENT(IN)   :: t 
    REAL(KIND=8), DIMENSION(:),          INTENT(IN)   :: mu_H_field
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv
    REAL(KIND=8)                                      :: PI, Lz, ri, ro, k0
    REAL(KIND=8),DIMENSION(SIZE(rr,2))                :: r, z, f, fp, amp_h


    IF (test_de_convergence) THEN
       vv = Hexact_anal(TYPE, rr, m, mu_H_field, t)
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

  FUNCTION Phiexact(TYPE, rr, m, mu_phi,t) RESULT(vv) 

    USE boundary_anal
    USE bessel
    IMPLICIT NONE

    INTEGER     ,                        INTENT(IN)   :: TYPE
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m  !mode 
    REAL(KIND=8),                        INTENT(IN)   :: mu_phi, t
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv   
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: r, z
    INTEGER                                           :: n
    REAL(KIND=8)                                      :: t_ext, t_al

    IF (test_de_convergence) THEN
       vv = Phiexact_anal(TYPE, rr, m, mu_phi,t)
       RETURN
    END IF

    vv(:) = 0.d0
    !
    !! champ horizontal Bx en extinction et champ vertical Bz constant
    !
    t_al = 4.d-1 ! 1.d-2 
    t_ext = 2.d0 ! 1.d-2 

    IF ((m ==1).AND.(TYPE==1)) THEN
       vv = 0.05d0*rr(1,:)*&
            (t/t_al)**3/(1.d0+(t/t_al)**3)*&
            (1.d0-(t/t_ext)**4/(1.d0+(t/t_ext)**4))
             ELSE IF ((m ==0).AND.(type==1)) THEN
       !           vv = 0.05d0*rr(2,:)/(1.d0+(t/t_ext)**3)
                   vv = 0.05d0*rr(2,:)*&
            (t/t_al)**3/(1.d0+(t/t_al)**3)*&
            (1.d0-(t/t_ext)**4/(1.d0+(t/t_ext)**4))
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
    m_max_c = SIZE(list_mode)
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
       IF (list_mode(i)==0) THEN
          Hn(:,3,i)=  1.d-5
       END IF
       ENDDO
    ENDDO


  END SUBROUTINE init_maxwell

END MODULE boundary
