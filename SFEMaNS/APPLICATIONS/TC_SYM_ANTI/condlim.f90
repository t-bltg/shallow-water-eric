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

  FUNCTION vv_exact(type,rr,m,t) RESULT(vv)
    
    USE boundary_anal

    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: type
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER,                             INTENT(IN)   :: m  !mode 
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8)                                      :: Re 
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv, r, z 
    INTEGER :: n
    REAL(KIND=8)                                      :: PI, k, ri, ro, eps

    IF (test_de_convergence) THEN
       vv = vv_exact_anal(type,rr,m,t) 
       RETURN
    END IF

  !Taylor couette avec cylindre interieur anime
    Ri = 1.d0
    Ro = 2.d0
    eps = 1.d-4
    vv = 0.d0
    IF(m==0 .and. type == 3) THEN
      DO n = 1, SIZE(rr,2)
         IF(rr(1,n) > Ro - eps ) THEN
           vv(n) = 0
         ELSE IF(rr(1,n) < Ri + eps) THEN
           vv(n) =  t**2/(0.1+t**2) !           1.d0
         END IF
      END DO
    END IF
    RETURN

    !Von Karmann
    vv = 0.d0
    IF(m==0 .and. type == 3) THEN
      DO n = 1, SIZE(rr,2)
         IF(rr(1,n) > 0.49999d0) THEN
           vv(n) = 0
         ELSE IF(rr(2,n) < 0.001d0) THEN
           vv(n) = rr(1,n)
         ELSE
           vv(n) = -rr(1,n)
         END IF
      END DO
    END IF

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

    vv=0
   END FUNCTION pp_exact

  ! Extension of the velocity field in the solid
  FUNCTION extension_vel(type,rr,m,t) RESULT(vv)

    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: type
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER,                             INTENT(IN)   :: m  !mode 
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8)                                      :: Re 
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv, r, z 
    INTEGER :: n
    REAL(KIND=8)                                      :: PI, k, Ri

    Ri = 1.d0
    vv = 0.d0
    IF(m==0 .and. type == 3) THEN
      DO n = 1, SIZE(rr,2)
         vv(n) =  t**2/(0.1+t**2)*rr(1,n)/Ri !           1.d0
      END DO
    END IF
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

  END FUNCTION Vexact

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
   !CHAMP J. LEORAT
   PI = ACOS(-1.d0)
   r = rr(1,:)
   z = rr(2,:)

    Lz = 2.d0*PI
    ri = 1.d0                         !rayon du cylindre interieur
    ro = 2.d0                         !rayon du cylindre exterieur
    k0 = 2*2.d0*PI/Lz
    f  = (ri-r)**2*(ro-r)**2
    fp = -2.d0*((ri-r)**2*(ro-r)+(ri-r)*(ro-r)**2)
    amp_h = 1d-1

    IF (m == 1) THEN
       IF (type == 1) THEN
          vv = f/r*SIN(k0*z)      !f/r*COS(k0*z)                     !f/r*SIN(k0*z)
       ELSEIF (type == 4) THEN
          vv = -fp*SIN(k0*z)      !-fp*COS(k0*z)                     !-fp*SIN(k0*z)
       ENDIF
    ENDIF
    vv = amp_h * vv

    
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

    IF (test_de_convergence) THEN
       vv = Phiexact_anal(type, rr, m, mu_phi,t)
       RETURN
    END IF
 
    vv(:) = 0.d0

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
    INTEGER :: i, k, mode, m_max_c
    
    IF (test_de_convergence) THEN
       CALL init_maxwell_analytique(H_mesh,phi_mesh,time,dt,mu_H_field,mu_phi,list_mode,Hn1,Hn,phin1,phin) 
       RETURN
    END IF

    time = -dt
    m_max_c = size(list_mode)
    DO k=1,6
      DO i=1, m_max_c
         mode = list_mode(i)
         Hn1(:,k,i) = Hexact(k, H_mesh%rr, mode, mu_H_field, time)
         IF (k<3) THEN
            phin1(:,k,i) = Phiexact(k, phi_mesh%rr, mode, mu_phi, time)
         ENDIF
      ENDDO
    ENDDO

    time = time + dt
    DO k=1,6
      DO i=1, m_max_c
         mode = list_mode(i)
         Hn(:,k,i) = Hexact(k, H_mesh%rr, mode, mu_H_field, time)
         IF (k<3) THEN
            phin(:,k,i) = Phiexact(k, phi_mesh%rr, mode, mu_phi, time)
         ENDIF
      ENDDO
    ENDDO

!   Hn1= 0.d0
!   DO i = 1, SIZE(list_mode)
!      IF (list_mode(i)==0) THEN
!         CALL RANDOM_NUMBER(Hn1(:,3,i))
!         Hn1(:,3,i)= 1.d-4*Hn1(:,3,i)
!      END IF
!   END DO
!   Hn = 0.d0
!   phin  = 0.d0
!   phin1  = 0.d0
!   time  = 0.d0


  END SUBROUTINE init_maxwell

END MODULE boundary
