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
    INTEGER                                    :: i,k,m_max_c,mode    

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

  
    time = -dt
    m_max_c = size(list_mode)
    DO k=1,6
      DO i=1, m_max_c
         mode = list_mode(i)
         un_m1(:,k,i) = vv_exact(k,mesh_f%rr, mode, time) ! vv_exact_init(k,mesh_f%rr, mode, time)
         IF (k<3) THEN
            pn_m1(:,k,i) = pp_exact(k, mesh_c%rr, mode, time)
         ENDIF
      ENDDO
    ENDDO

    time = time + dt
    DO k=1,6
      DO i=1, m_max_c
         mode = list_mode(i)
         un(:,k,i) = vv_exact(k,mesh_f%rr, mode, time) ! vv_exact_init(k,mesh_f%rr, mode, time)
         IF (k<3) THEN
            pn(:,k,i) = pp_exact(k, mesh_c%rr, mode, time)
         ENDIF
      ENDDO
    ENDDO

         
  END SUBROUTINE init_up

  FUNCTION ff_source(type,rr,m,t, Re, ty) RESULT(vv)
    !=================================== 
    USE boundary_anal
    USE chaine_caractere
    
    IMPLICIT NONE

    INTEGER     ,                        INTENT(IN)   :: type
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv
    INTEGER     ,                        INTENT(IN)   :: m  !mode  
    REAL(KIND=8),                        INTENT(IN)   :: t   
    INTEGER                                           :: n
    CHARACTER(LEN=2),                    INTENT(IN)   :: ty
    REAL(KIND=8),                        INTENT(IN)   :: Re
    REAL(KIND=8), SAVE                                      :: H0,S0,R0
    REAL(KIND=8)                                      :: omg0, t_al
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: X
    INTEGER                                           :: i
    LOGICAL, SAVE                                     :: once=.true.


    IF (test_de_convergence) THEN
       vv = ff_anal(type,rr,m,t,Re,ty)
       RETURN
    END IF

    vv = 0.d0
   IF (once) THEN
   once=.FALSE.
    OPEN(UNIT = 21, FILE = 'data', FORM = 'formatted', STATUS = 'unknown')
          CALL read_until(21, 'data_MRI_KM')
          READ(21,*) H0,S0,R0
         !READ(21,*) S0,H0,R0 CN 10/02/08 BUG
    CLOSE(21)
   ENDIF

   
   X=(rr(1,:)/s0)
   omg0=1.d0
   t_al=5.d-1

   IF (type==3 .AND. m==0) THEN
      ! vv = omg0*rr(1,:) ! SPIN_UP
     ! vv=(t/t_al)**3/(1.d0 + (t/t_al)**3)* &
      !   (omg0/(s0*Re))*7.5d0*X**2*(1.d0+ 1.d-1*X**3)*((1.d0+X**3)**(-2.5))
      !KM source + allumage
      vv = (omg0/(s0*Re))*7.5d0*X**2*(1.d0+ 1.d-1*X**3)*((1.d0+X**3)**(-2.5))
   END IF
  
  RETURN

  END FUNCTION ff_source

   FUNCTION vv_exact(type, rr, mode, t) RESULT(vv)

   USE boundary_anal
   USE def_type_mesh
   USE sub_plot
   USE chaine_caractere
   IMPLICIT NONE

   INTEGER     ,                        INTENT(IN)   :: type
   REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
   INTEGER,                             INTENT(IN)   :: mode  !mode
   REAL(KIND=8),                        INTENT(IN)   :: t
   REAL(KIND=8), DIMENSION(size(rr,2))               :: vv

   REAL(KIND=8)                                      :: t_al, omg0
   REAL(KIND=8), SAVE                                :: H0,S0,R0
   LOGICAL, SAVE                                     :: once=.true.

   IF (test_de_convergence) THEN
      vv = vv_exact_anal(type,rr,mode,t)
      RETURN
   END IF
   
   vv = 0.d0 

   IF (once) THEN
   once=.FALSE.
    OPEN(UNIT = 21, FILE = 'data', FORM = 'formatted', STATUS = 'unknown')
          CALL read_until(21, 'data_MRI_KM')
          READ(21,*) H0,S0,R0
         !READ(21,*) S0,H0,R0 CN 10/02/08 BUG
    CLOSE(21)
   ENDIF

   omg0=1.d0
   t_al=5.d-1

   IF (type==3 .AND. mode==0) THEN
         ! vv(:)=omg0*rr(1,:)*(t/t_al)**3/(1.d0 + (t/t_al)**3) ! SPIN_UP+allumage
         ! vv(:)=omg0*(rr(1,:)/(1.d0+(rr(1,:)/S0)**3)**(0.5))*&
         !       (t/t_al)**3/(1.d0 + (t/t_al)**3) ! KM97 + allumage
           vv(:)=omg0*(rr(1,:)/(1.d0+(rr(1,:)/S0)**3)**(0.5)) ! KM97
         ! vv(:)=omg0*rr(1,:) ! Rot Sol
   END IF

END FUNCTION vv_exact

   FUNCTION vv_exact_init(type, rr, mode, t) RESULT(vv)

   USE boundary_anal
   USE def_type_mesh
   USE sub_plot
   USE chaine_caractere
   IMPLICIT NONE

   INTEGER     ,                        INTENT(IN)   :: type
   REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
   INTEGER,                             INTENT(IN)   :: mode  !mode
   REAL(KIND=8),                        INTENT(IN)   :: t
   REAL(KIND=8), DIMENSION(size(rr,2))               :: vv

   REAL(KIND=8)                                      :: t_al, omg0, PI
   REAL(KIND=8), SAVE                                :: H0,S0,R0
   LOGICAL, SAVE                                     :: once=.true.

   IF (test_de_convergence) THEN
      vv = vv_exact_anal(type,rr,mode,t)
      RETURN
   END IF

   vv = 0.d0

   IF (once) THEN
   once=.FALSE.
    OPEN(UNIT = 21, FILE = 'data', FORM = 'formatted', STATUS = 'unknown')
          CALL read_until(21, 'data_MRI_KM')
          READ(21,*) H0,S0,R0
         !READ(21,*) S0,H0,R0 CN 10/02/08 BUG
    CLOSE(21)
   ENDIF


   PI = ACOS(-1.d0)
   omg0=1.d0
   t_al=5.d-1

   IF (type==3 .AND. mode==0) THEN
           vv(:)=0.0d0*omg0*cos(3*PI*rr(2,:))*sin(2*PI*rr(1,:)/R0)+ &
                 omg0*(rr(1,:)/(1.d0+(rr(1,:)/S0)**3)**(0.5)) ! KM97
   END IF

END FUNCTION vv_exact_init

  FUNCTION pp_exact(type,rr,m,t) RESULT (vv)
    
    USE boundary_anal
    USE chaine_caractere

    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: type
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv 
    INTEGER     ,                        INTENT(IN)   :: m  !mode
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8)                                      :: omg0
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: X
    INTEGER                                           :: i
    REAL(KIND=8), SAVE                                :: H0,S0,R0
    LOGICAL, SAVE                                     :: once=.true.

    IF (test_de_convergence) THEN
       vv = pp_exact_anal(type,rr,m,t)
       RETURN
    END IF

    vv=0.d0


   IF (once) THEN
   once=.FALSE.
    OPEN(UNIT = 21, FILE = 'data', FORM = 'formatted', STATUS = 'unknown')
          CALL read_until(21, 'data_MRI_KM')
          READ(21,*) H0,S0,R0
         !READ(21,*) S0,H0,R0 CN 10/02/08 BUG
    CLOSE(21)
   ENDIF

    omg0=1.d0

X=rr(1,:)/S0

 IF (type==1 .AND. m==0) THEN
  !      +0.5d0*(omg0*(rr(1,:)/(1.d0+(rr(1,:)/S0)**3)**(0.5)))**2 
!pression dynamique KM97 02/04/08
        vv(:)=(omg0*S0)**2*(atan((2*X(:)-1.d0)/sqrt(3.d0))/sqrt(3.d0) &
        +(1.d0/6.d0)*log(( X(:)**2 -X(:)+ 1.d0)/((X(:) + 1.d0)**2))) &
        +0.5d0*(omg0*(rr(1,:)/(1.d0+(rr(1,:)/S0)**3)**(0.5)))**2 
!pression dynamique KM97 02/04/08

 END IF



   END FUNCTION pp_exact


  FUNCTION extension_vel(type,rr,m,t) RESULT(vv) ! Extension of the velocity field in the solid

    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: type
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER,                             INTENT(IN)   :: m  !mode 
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8)                                      :: Re 
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv, r, z 
    INTEGER :: n
    REAL(KIND=8)                                      :: PI, k, Ri

    vv=0.d0

  END FUNCTION extension_vel 



!============================================================================================
!                       CONDITIONS LIMITES POUR MAXWELL
!============================================================================================

  FUNCTION Vexact(m, H_mesh) RESULT(vv)  !Set uniquement a l'induction Terme non lineaire

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


    vv = 0.d0

  END FUNCTION Vexact

  FUNCTION Hexact_init(type, rr, m, mu_H_field, t) RESULT(vv) 
    
    USE chaine_caractere
    USE boundary_anal
    USE bessel
    IMPLICIT NONE

    INTEGER     ,                        INTENT(IN)   :: type
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m  !mode 
    REAL(KIND=8),                        INTENT(IN)   :: t 
    REAL(KIND=8), DIMENSION(:),          INTENT(IN)   :: mu_H_field
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv
    REaL(KIND=8)                                      :: PI
    REAL(KIND=8),DIMENSION(SIZE(rr,2))                :: r, z, f, fp, amp_h
    REAL(KIND=8), SAVE                                :: H0,S0,R0
    LOGICAL, SAVE                                     :: once=.true.
    
    IF (test_de_convergence) THEN
       vv = Hexact_anal(type, rr, m, mu_H_field, t)
       RETURN
    END IF

   IF (once) THEN
   once=.FALSE.
    OPEN(UNIT = 21, FILE = 'data', FORM = 'formatted', STATUS = 'unknown')
          CALL read_until(21, 'data_MRI_KM')
          READ(21,*) H0,S0,R0
         !READ(21,*) S0,H0,R0 CN 10/02/08 BUG
    CLOSE(21)
   ENDIF

   PI=acos(-1.d0)

    vv(:) = 0.d0
       
   IF(type==5 .AND. m==0) THEN
    vv(:)=H0
   END IF
        IF ((m==0).AND.(type==3)) THEN
       vv(:)=0.01*H0*cos(PI*rr(2,:))*sin(2*PI*rr(1,:)/R0)
      !vv = 0.01d0*H0*sin(PI*rr(1,:)/ro)*sin(PI*rr(1,:)/ro)
      !vv = 0.05d0*cos(PI*rr(2,:)/Lz)*sin(PI*rr(1,:)/ro)
    ENDIF
    
  END FUNCTION Hexact_init

  FUNCTION Hexact(type, rr, m, mu_H_field, t) RESULT(vv) 
    
    USE chaine_caractere
    USE boundary_anal
    USE bessel
    IMPLICIT NONE

    INTEGER     ,                        INTENT(IN)   :: type
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m  !mode 
    REAL(KIND=8),                        INTENT(IN)   :: t 
    REAL(KIND=8), DIMENSION(:),          INTENT(IN)   :: mu_H_field
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv
    REAL(KIND=8)                                      :: PI, Lz, ri, k0, eps
    REAL(KIND=8),DIMENSION(SIZE(rr,2))                :: r, z, fsr, fsrp
    INTEGER                                           :: i
    REAL(KIND=8), SAVE                                :: H0,S0,R0
    LOGICAL, SAVE                                     :: once=.true.

    IF (test_de_convergence) THEN
       vv = Hexact_anal(type, rr, m, mu_H_field, t)
       RETURN
    END IF


   IF (once) THEN
   once=.FALSE.
    OPEN(UNIT = 21, FILE = 'data', FORM = 'formatted', STATUS = 'unknown')
          CALL read_until(21, 'data_MRI_KM')
          READ(21,*) H0,S0,R0
         !READ(21,*) S0,H0,R0 CN 10/02/08 BUG
    CLOSE(21)
   ENDIF
   
   !CHAMP Uniforme Hz 

   PI=acos(-1.d0)
   vv(:)=0.d0
  
   IF(type==5 .AND. m==0) THEN
    vv(:)=H0
   END IF
  
   !CHAMP de pertubation en r
   
   !eps=1.d-3
   !fsr=rr(1,:)**2*(ro-rr(1,:))**2
   !fsrp=rr(1,:)*(ro-rr(1,:))*(3*ro-5*rr(1,:))
   
    ! IF (type==1 .AND. m==0) THEN
    !   vv=amp*fsr*(2*PI*sin(2*PI*rr(2,:)))
    ! END IF

     !IF (type==5 .AND. m==0) THEN
     !  vv=amp*fsrp*(1+cos(2*PI*rr(2,:)))
     !END IF
   
  !CHAMP de pertubation en Hth
 
!   IF (type==3 .AND. m==0) THEN
!     vv(:)=0.01*H0*cos(PI*rr(2,:))*sin(2*PI*rr(1,:)/R0)

!   END IF
 
  END FUNCTION Hexact

  FUNCTION Phiexact(type, rr, m, mu_phi,t) RESULT(vv) 
     
    USE boundary_anal
    USE bessel
    USE chaine_caractere
    IMPLICIT NONE

    INTEGER     ,                        INTENT(IN)   :: type
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m  !mode 
    REAL(KIND=8),                        INTENT(IN)   :: mu_phi, t
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv   
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: r, z
    INTEGER                                           :: ni
    REAL(KIND=8)                                      :: PI, t_loc
    REAL(KIND=8), SAVE                                :: H0,S0,R0,t_al,t_ext,t_restart
    LOGICAL, SAVE                                     :: once=.true.

    IF (test_de_convergence) THEN
       vv = Phiexact_anal(type, rr, m, mu_phi,t)
       RETURN
    END IF
   
    vv(:)=0.d0

   IF (once) THEN
   once=.FALSE.
    OPEN(UNIT = 21, FILE = 'data', FORM = 'formatted', STATUS = 'unknown')
           CALL read_until(21, 'data_MRI_KM')
           READ(21,*) H0,S0,R0,t_al,t_ext,t_restart
    CLOSE(21)
   ENDIF
!allumage d un champ vertical

   PI=acos(-1.d0)
   IF (type==1 .AND. m==0) THEN
     vv(:)=H0*(rr(2,:))
   END IF
!!allumage d un champ oblique incline 85 degres par phi (ds data, il faut mettre .t. pour
! !interface 3)
!    IF ((m ==1).AND.(TYPE==1)) THEN
!          t_loc  = t - t_restart
!          vv =(H0/10.d0)*rr(1,:)* &
!            (t_loc/t_al)**3/(1.d0+(t_loc/t_al)**3)*&
!            (1.d0-(t_loc/t_ext)**5/(1.d0+(t_loc/t_ext)**5))
!    ENDIF



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
    REAL(KIND=8), DIMENSION(:),   ALLOCATABLE  :: x
   

    ALLOCATE(x(H_mesh%np))
 
    IF (test_de_convergence) THEN
       CALL init_maxwell_analytique(H_mesh,phi_mesh,time,dt,mu_H_field,mu_phi,list_mode,Hn1,Hn,phin1,phin) 
       RETURN
    END IF

     time = -dt
    m_max_c = size(list_mode)
    Hn1= 0.d0
    phin1=0.d0
       DO i = 1, m_max_c

          IF(list_mode(i)==0) THEN
                 Hn1(:,5,i) = Hexact_init(5, H_mesh%rr, list_mode(i), mu_H_field, time)
                 Hn1(:,3,i) = Hexact_init(3, H_mesh%rr, list_mode(i), mu_H_field, time)
                 phin1(:,1,i) = Phiexact(1, phi_mesh%rr,list_mode(i), mu_phi, time)
          END IF

          
        IF(list_mode(i)>0 .AND. list_mode(i)< 3) THEN
             CALL RANDOM_NUMBER(Hn1(:,1,i))
             CALL RANDOM_NUMBER(Hn1(:,2,i))
             CALL RANDOM_NUMBER(Hn1(:,3,i))
             CALL RANDOM_NUMBER(Hn1(:,4,i))
             CALL RANDOM_NUMBER(Hn1(:,5,i))
             CALL RANDOM_NUMBER(Hn1(:,6,i))
              Hn1(:,:,i)= 1.d-5*(Hn1(:,:,i)-0.5d0)
          END IF

     END DO


    Hn = Hn1
    phin  =  phin1
    time  = 0.d0



  END SUBROUTINE init_maxwell


END MODULE boundary
