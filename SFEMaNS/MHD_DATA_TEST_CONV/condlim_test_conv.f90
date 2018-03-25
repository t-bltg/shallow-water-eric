MODULE boundary_anal
!TEST June 17 2008
  USE def_type_mesh
  TYPE(mesh_type), POINTER   :: H_mesh_analytic
!TEST June 17 2008
  LOGICAL       :: test_conv_ns, test_conv_maxwell
!jan 29 2007
  LOGICAL       :: dom_H_larger_dom_ns
  LOGICAL       :: ns_periodic, mxw_periodic 
  INTEGER       :: no_test ! no_test = 1 -> non periodique,  no_test = 2 -> periodique
!jan 29 2007
  REAL (KIND=8), PRIVATE :: alpha=1.d0, beta=1.d0, gamma = 1.d0, x0=3, y0=0
  
!jan 29 2007
   !parametre solution perio maxwell
   REAL(KIND=8), PRIVATE  :: r0 = 0.5d0
   REAL(KIND=8), PRIVATE  :: k0=6.28318530717958d0
   INTEGER     , PRIVATE  :: m0=1
!jan 29 2007
CONTAINS

!============================================================================================
!                       CONDITIONS LIMITES POUR NAVIER_STOKES
!============================================================================================
  SUBROUTINE init_up_analytique(mesh_f,mesh_c,time,dt,list_mode,un_m1,un,pn_m1,pn,phin_m1,phin)

    USE def_type_mesh

    IMPLICIT NONE
    TYPE(mesh_type)                              :: mesh_f, mesh_c     
    REAL(KIND=8),                   INTENT(OUT)  :: time
    REAL(KIND=8),                   INTENT(IN)   :: dt
    INTEGER,      DIMENSION(:),     INTENT(IN)   :: list_mode    
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT)  :: un_m1, un
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT)  :: pn_m1, pn, phin_m1, phin
    INTEGER                                      :: i,j, mode  
    REAL(KIND=8), DIMENSION(mesh_c%np)           :: pn_m2
    !April 17th 2008, JLG
    REAL(KIND=8) zero
    DATA zero/0.d0/
    !April 17th 2008, JLG

    DO i= 1, SIZE(list_mode)
       mode = list_mode(i) 
       DO j = 1,6 

          !vitesse
          un_m1(:,j,i) = vv_exact_anal(j,mesh_f%rr,mode,-dt)  
          un   (:,j,i) = vv_exact_anal(j,mesh_f%rr,mode,zero)
          !pression
          IF (j == 1 .or. j == 2) THEN
             pn_m2(:)       = pp_exact_anal(j,mesh_c%rr,mode ,-2*dt)
             pn_m1  (:,j,i) = pp_exact_anal(j,mesh_c%rr,mode,zero)
             pn     (:,j,i) = pp_exact_anal(j,mesh_c%rr,mode,dt)
             phin_m1(:,j,i) = pn_m1(:,j,i) - pn_m2(:)
             phin   (:,j,i) = Pn   (:,j,i) - pn_m1(:,j,i)
          ENDIF

       ENDDO
    ENDDO

    time = 0.d0

  END SUBROUTINE init_up_analytique

  FUNCTION ff_anal(type,rr,m,t,Re,ty) RESULT(vv)
    !=================================== 
    IMPLICIT NONE

    INTEGER     ,                        INTENT(IN)   :: type
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv
    INTEGER     ,                        INTENT(IN)   :: m  !mode  
    REAL(KIND=8),                        INTENT(IN)   :: t,Re
    CHARACTER(LEN=2),                    INTENT(IN)   :: ty   

    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: LAP, TEMP,GPRESS
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: r, z
    REAL(KIND=8), DIMENSION(SIZE(rr,2),6)             :: Rot, V, vt  
    INTEGER :: i, j, m1
    REAL(KIND=8)                                      :: PI, k
    REAL(KIND=8)                                      :: one
    DATA one/1.d0/
    IF (test_conv_ns .AND. no_test==1) THEN
     
         IF (m==0 .OR. m==2) THEN
          IF(ty=='st') THEN !
             vv =0
             RETURN
          END IF

          !Contribution nonlineaire
          m1 = 1
          Rot(:,1) = m1*rr(1,:)*rr(2,:)**3*(4+m1)  &
               -6*rr(1,:)**3*rr(2,:)+9*(rr(1,:)*rr(2,:))**2
          Rot(:,2) =-m1*rr(1,:)*rr(2,:)**3*(4-m1)&
               -6*rr(1,:)**3*rr(2,:)+9*(rr(1,:)*rr(2,:))**2          
          Rot(:,3) = 3*m1*(rr(1,:)*rr(2,:))**2-6*rr(1,:)**3*rr(2,:)  &
               -2*rr(1,:)*rr(2,:)**3*(4-m1)
          Rot(:,4) =-3*m1*(rr(1,:)*rr(2,:))**2-6*rr(1,:)**3*rr(2,:)  &
               -2*rr(1,:)*rr(2,:)**3*(4+m1)
          Rot(:,5) = 12*(rr(1,:)*rr(2,:))**2-9*rr(1,:)*rr(2,:)**3 + &
               m1*(m1*rr(1,:)*rr(2,:)**3+3*(rr(1,:)*rr(2,:))**2)
          Rot(:,6) = 12*(rr(1,:)*rr(2,:))**2-9*rr(1,:)*rr(2,:)**3 + &
               m1*(m1*rr(1,:)*rr(2,:)**3-3*(rr(1,:)*rr(2,:))**2)          
          Rot = Rot*COS(t)

          Vt(:,1) = m1*rr(2,:)**3 - 3*rr(1,:)*rr(2,:)**2
          Vt(:,4) =3*(rr(1,:)-rr(2,:))*rr(2,:)**2
          Vt(:,5) =(4-m1)*rr(2,:)**3
          V(:,1) = vt(:,1)*rr(1,:)**2
          V(:,4) = vt(:,4)*rr(1,:)**2
          V(:,5) = vt(:,5)*rr(1,:)**2

          Vt(:,2) =-m1*rr(2,:)**3 - 3*rr(1,:)*rr(2,:)**2
          Vt(:,3) =3*(rr(1,:)-rr(2,:))*rr(2,:)**2
          Vt(:,6) =(4+m1)*rr(2,:)**3
          V(:,2) = vt(:,2)*rr(1,:)**2
          V(:,3) = vt(:,3)*rr(1,:)**2
          V(:,6) = vt(:,6)*rr(1,:)**2
          V = V*COS(t)

          IF (m==0) THEN ! vv = -VxRot(v) pour mode 0
             IF (MOD(type,2)==0) THEN
                vv = 0.d0
             ELSEIF (type == 1) THEN
                vv = -0.5d0*(V(:,3)*Rot(:,5)-V(:,5)*Rot(:,3)  &
                     + V(:,4)*Rot(:,6)-V(:,6)*Rot(:,4))
             ELSEIF (type == 3) THEN
                vv = -0.5d0*(V(:,5)*Rot(:,1)-V(:,1)*Rot(:,5)  &
                     + V(:,6)*Rot(:,2)-V(:,2)*Rot(:,6))
             ELSEIF (type == 5) THEN
                vv = -0.5d0*(V(:,1)*Rot(:,3)-V(:,3)*Rot(:,1)  &
                     + V(:,2)*Rot(:,4)-V(:,4)*Rot(:,2))        
             ENDIF
             RETURN
          END IF

          IF (m==2) THEN ! vv = -VxRot(v) pour mode 2
             IF (type == 1) THEN
                vv = -0.5d0*(V(:,3)*Rot(:,5)-V(:,5)*Rot(:,3)  &
                     - (V(:,4)*Rot(:,6)-V(:,6)*Rot(:,4)))
             ELSE IF (type == 2) THEN
                vv = -0.5d0*(V(:,3)*Rot(:,6)-V(:,5)*Rot(:,4)  &
                     + (V(:,4)*Rot(:,5)-V(:,6)*Rot(:,3)))
             ELSEIF (type == 3) THEN
                vv = -0.5d0*(V(:,5)*Rot(:,1)-V(:,1)*Rot(:,5)  &
                     - (V(:,6)*Rot(:,2)-V(:,2)*Rot(:,6)))
             ELSEIF (type == 4) THEN
                vv = -0.5d0*(V(:,5)*Rot(:,2)-V(:,1)*Rot(:,6)  &
                     + (V(:,6)*Rot(:,1)-V(:,2)*Rot(:,5)))
             ELSEIF (type == 5) THEN
                vv = -0.5d0*(V(:,1)*Rot(:,3)-V(:,3)*Rot(:,1)  &
                     - (V(:,2)*Rot(:,4)-V(:,4)*Rot(:,2)))
             ELSEIF (type == 6) THEN
                vv = -0.5d0*(V(:,1)*Rot(:,4)-V(:,3)*Rot(:,2)  &
                     + (V(:,2)*Rot(:,3)-V(:,4)*Rot(:,1)))
             END IF

             RETURN
          END IF

       END IF
       IF (m/=1) THEN
          vv = 0 
          RETURN
       END IF

       !Calcul pour le mode 1
       Vt(:,1) = m*rr(2,:)**3 - 3*rr(1,:)*rr(2,:)**2
       Vt(:,4) =3*(rr(1,:)-rr(2,:))*rr(2,:)**2
       Vt(:,5) =(4-m)*rr(2,:)**3
       V(:,1) = vt(:,1)*rr(1,:)**2
       V(:,4) = vt(:,4)*rr(1,:)**2
       V(:,5) = vt(:,5)*rr(1,:)**2

       Vt(:,2) =-m*rr(2,:)**3 - 3*rr(1,:)*rr(2,:)**2
       Vt(:,3) =3*(rr(1,:)-rr(2,:))*rr(2,:)**2
       Vt(:,6) =(4+m)*rr(2,:)**3
       V(:,2) = vt(:,2)*rr(1,:)**2
       V(:,3) = vt(:,3)*rr(1,:)**2
       V(:,6) = vt(:,6)*rr(1,:)**2

       IF (type == 1) THEN
          LAP(:) = m*4*rr(2,:)**3 - 27*rr(1,:)*rr(2,:)**2 + m*6*rr(1,:)**2*rr(2,:) - 6*rr(1,:)**3 &
               -m**2*vt(:,1) - 2*m*vt(:,4) - vt(:,1)
          TEMP(:) = V(:,1)
          GPRESS(:) = 2*rr(1,:)*rr(2,:)**3*COS(t)
       ELSEIF (type == 2) THEN
          LAP(:) =-m*4*rr(2,:)**3 - 27*rr(1,:)*rr(2,:)**2 - m*6*rr(1,:)**2*rr(2,:) - 6*rr(1,:)**3 &
               -m**2*vt(:,2) + 2*m*vt(:,3) - vt(:,2)
          TEMP(:) = V(:,2)
          GPRESS(:) = 2*rr(1,:)*rr(2,:)**3*COS(t)
       ELSEIF (type == 3) THEN
          LAP(:) = 3*(9*rr(1,:)*rr(2,:)**2 - 4*rr(2,:)**3  + 2*rr(1,:)**3- 6*rr(1,:)**2*rr(2,:)) &
               -m**2*vt(:,3) + 2*m*vt(:,2) -vt(:,3)
          TEMP(:) = V(:,3)
          GPRESS(:) = m*rr(1,:)*rr(2,:)**3*COS(t)
       ELSEIF (type == 4) THEN
          LAP(:) = 3*(9*rr(1,:)*rr(2,:)**2 - 4*rr(2,:)**3  + 2*rr(1,:)**3- 6*rr(1,:)**2*rr(2,:)) &
               -m**2*vt(:,4) -2*m*vt(:,1) -vt(:,4)
          TEMP(:) = V(:,4)
          GPRESS(:) = -m*rr(1,:)*rr(2,:)**3*COS(t)
       ELSEIF (type == 5) THEN
          LAP(:) = (4-m)*(4*rr(2,:)**3  + 6*rr(1,:)**2*rr(2,:)) -m**2*vt(:,5)
          TEMP(:) = V(:,5)
          GPRESS(:) = 3*rr(1,:)**2*rr(2,:)**2*COS(t)
       ELSEIF (type == 6) THEN
          LAP(:) = (4+m)*(4*rr(2,:)**3  + 6*rr(1,:)**2*rr(2,:)) -m**2*vt(:,6)
          TEMP(:) = V(:,6)
          GPRESS(:) = 3*rr(1,:)**2*rr(2,:)**2*COS(t)
       ENDIF
       vv(:) = -SIN(t)*TEMP  + COS(t)*(-LAP/Re) + GPRESS
       RETURN
    END IF

    IF (test_conv_ns .AND. no_test==2) THEN 
       PI = ACOS(-one) 
       k  = 2*PI 
       r = rr(1,:)
       z = rr(2,:) 

       IF (m==0 .OR. m==2) THEN
          IF(ty=='st') THEN
             vv =0
             RETURN
          END IF
          !Contribution nonlineaire
          Rot(:,1) = r*(4.d0*COS(k*z)+1) 
          Rot(:,2) = 0.d0
          Rot(:,3) = 0.d0
          Rot(:,4) = r*(k**2*r**2*COS(k*z)-8.d0*COS(k*z)-2.d0)
          Rot(:,5) = r*(-8.d0-k*r*SIN(k*z))
          Rot(:,6) = 0.d0
          Rot = Rot*COS(t)

          Vt(:,1) = 0.d0
          Vt(:,4) = 0.d0
          Vt(:,5) = 0.d0
          Vt(:,2) = -r**2*(1.d0-k*r*SIN(k*z))
          Vt(:,3) = -3.d0*r**2
          Vt(:,6) = r**2*(4.d0*COS(k*z)+1.d0)
          V = Vt*COS(t)

          IF (m==0) THEN ! vv = -VxRot(v) pour mode 0
             IF (MOD(type,2)==0) THEN
                vv = 0.d0
             ELSEIF (type == 1) THEN
                vv = -0.5d0*(V(:,3)*Rot(:,5)-V(:,5)*Rot(:,3)  &
                     + V(:,4)*Rot(:,6)-V(:,6)*Rot(:,4))
             ELSEIF (type == 3) THEN
                vv = -0.5d0*(V(:,5)*Rot(:,1)-V(:,1)*Rot(:,5)  &
                     + V(:,6)*Rot(:,2)-V(:,2)*Rot(:,6))
             ELSEIF (type == 5) THEN
                vv = -0.5d0*(V(:,1)*Rot(:,3)-V(:,3)*Rot(:,1)  &
                     + V(:,2)*Rot(:,4)-V(:,4)*Rot(:,2))        
             ENDIF
             RETURN
          END IF

          IF (m==2) THEN ! vv = -VxRot(v) pour mode 2
             IF (type == 1) THEN
                vv = -0.5d0*(V(:,3)*Rot(:,5)-V(:,5)*Rot(:,3)  &
                     - (V(:,4)*Rot(:,6)-V(:,6)*Rot(:,4)))
             ELSE IF (type == 2) THEN
                vv = -0.5d0*(V(:,3)*Rot(:,6)-V(:,5)*Rot(:,4)  &
                     + (V(:,4)*Rot(:,5)-V(:,6)*Rot(:,3)))
             ELSEIF (type == 3) THEN
                vv = -0.5d0*(V(:,5)*Rot(:,1)-V(:,1)*Rot(:,5)  &
                     - (V(:,6)*Rot(:,2)-V(:,2)*Rot(:,6)))
             ELSEIF (type == 4) THEN
                vv = -0.5d0*(V(:,5)*Rot(:,2)-V(:,1)*Rot(:,6)  &
                     + (V(:,6)*Rot(:,1)-V(:,2)*Rot(:,5)))
             ELSEIF (type == 5) THEN
                vv = -0.5d0*(V(:,1)*Rot(:,3)-V(:,3)*Rot(:,1)  &
                     - (V(:,2)*Rot(:,4)-V(:,4)*Rot(:,2)))
             ELSEIF (type == 6) THEN
                vv = -0.5d0*(V(:,1)*Rot(:,4)-V(:,3)*Rot(:,2)  &
                     + (V(:,2)*Rot(:,3)-V(:,4)*Rot(:,1)))
             END IF

             RETURN
          END IF

       END IF
       IF (m/=1) THEN
          vv = 0 
          RETURN
       END IF

       Vt(:,1) = 0.d0
       Vt(:,4) = 0.d0
       Vt(:,5) = 0.d0
       Vt(:,2) = -r**2*(1.d0-k*r*SIN(k*z))
       Vt(:,3) = -3.d0*r**2
       Vt(:,6) = r**2*(4.d0*COS(k*z)+1.d0)

       IF (type == 1) THEN
          LAP(:) =  0.d0
          TEMP(:) = Vt(:,1)
          GPRESS(:) = 2*r*COS(k*z)
       ELSEIF (type == 2) THEN
          LAP(:) =  (-8+7*k*r*SIN(k*z)-k**3*r**3*SIN(k*z))
          TEMP(:) = Vt(:,2)
          GPRESS(:) = 0.d0
       ELSEIF (type == 3) THEN
          LAP(:) = (-8+2*k*r*SIN(k*z))
          TEMP(:) = Vt(:,3)
          GPRESS(:) = 0.d0
       ELSEIF (type == 4) THEN
          LAP(:) =  0.d0
          TEMP(:) = Vt(:,4)
          GPRESS(:) = -r*COS(k*z)
       ELSEIF (type == 5) THEN
          LAP(:) = 0.d0
          TEMP(:) = Vt(:,5)
          GPRESS(:) = -r**2*k*SIN(k*z)
       ELSEIF (type == 6) THEN
          LAP(:) = (3+12*COS(k*z)-4*k**2*r**2*COS(k*z))
          TEMP(:) = Vt(:,6)
          GPRESS(:) = 0.d0
       ENDIF

       vv(:) = -SIN(t)*TEMP  + COS(t)*(-LAP/Re) + GPRESS*cos(t)
       RETURN
   ENDIF

   WRITE(*,*) ' BUG in ff_anal, we shouldn''t arrive here'
   STOP

  END FUNCTION ff_anal

  FUNCTION vv_exact_anal(type,rr,m,t) RESULT(vv)

    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: type
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER,                             INTENT(IN)   :: m  !mode 
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8)                                      :: Re 
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv, r, z 
    INTEGER :: n
    REAL(KIND=8)                                      :: PI, k, ri, ro, eps
    REAL(KIND=8)                                      :: one
    DATA one/1.d0/

    IF (test_conv_ns .AND. no_test==1) THEN

       IF (m/=1) THEN
          vv = 0 
          RETURN
       END IF

       IF (type == 1) THEN
          vv(:) =  m*rr(1,:)**2*rr(2,:)**3 - 3*rr(1,:)**3*rr(2,:)**2
       ELSEIF (type == 2 .and. m /= 0) THEN
          vv(:) = -m*rr(1,:)**2*rr(2,:)**3 - 3*rr(1,:)**3*rr(2,:)**2
       ELSEIF (type == 3) THEN
          vv(:) = 3*rr(1,:)**3*rr(2,:)**2 - 3 *rr(1,:)**2*rr(2,:)**3
       ELSEIF (type == 4 .and. m /= 0) THEN
          vv(:) = 3*rr(1,:)**3*rr(2,:)**2 - 3 *rr(1,:)**2*rr(2,:)**3
       ELSEIF (type == 5) THEN
          vv(:) = rr(1,:)**2*rr(2,:)**3*(4-m)
       ELSEIF (type == 6 .and. m /= 0) THEN
          vv(:) = rr(1,:)**2*rr(2,:)**3*(4+m)
       ENDIF

       vv(:) = vv(:) * COS(t)
       RETURN
    END IF

    IF (test_conv_ns .AND. no_test==2) THEN
       IF (m/=1) THEN
          vv = 0 
          RETURN
       END IF

       PI= ACOS(-one)
       k = 2*PI
       r = rr(1,:)
       z = rr(2,:) 
       IF (type == 1) THEN
          vv(:) = 0.d0
       ELSEIF (type == 2 .and. m /= 0) THEN
          vv(:) = -r**2*(1-k*r*SIN(k*z))
       ELSEIF (type == 3) THEN
          vv(:) = -3*r**2 
       ELSEIF (type == 4 .and. m /= 0) THEN
          vv(:) = 0.d0
       ELSEIF (type == 5) THEN
          vv(:) = 0.d0
       ELSEIF (type == 6 .and. m /= 0) THEN
          vv(:) = r**2*(4*COS(k*z)+1)
       ENDIF

       vv(:) = vv(:) * COS(t)
       RETURN
    ENDIF

    WRITE(*,*) ' BUG in vv_exact_anal, we shouldn''t arrive here'
    STOP

  END FUNCTION vv_exact_anal

  FUNCTION pp_exact_anal(type,rr,m,t) RESULT (vv)
    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: type
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv 
    INTEGER     ,                        INTENT(IN)   :: m  !mode
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8)                                      :: PI, k
    REAL(KIND=8)                                      :: one
    DATA one/1.d0/
    

    IF (test_conv_ns .AND. no_test==1) THEN
       IF (m/=1) THEN
          vv = 0 
          RETURN
       END IF

       vv(:) = rr(1,:)**2*rr(2,:)**3*COS(t)
       RETURN
    END IF

    IF (test_conv_ns .AND. no_test==2)  THEN
       PI=ACOS(-one)
       k = 2*PI
       IF (m/=1) THEN
          vv = 0 
          RETURN
       END IF

       IF (type==1) THEN
          vv(:)= rr(1,:)**2*COS(k*rr(2,:))*COS(t)
       ELSE
          vv(:) = 0.d0
       END IF

       RETURN
    ENDIF

    WRITE(*,*) ' BUG in pp_exact_anal, we shouldn''t arrive here'
    STOP

  END FUNCTION pp_exact_anal


!============================================================================================
!                       CONDITIONS LIMITES POUR MAXWELL
!============================================================================================

  FUNCTION Vexact_anal(m, H_mesh) RESULT(vv)  !Set uniquement a l'induction
    USE def_type_mesh
    IMPLICIT NONE

    TYPE(mesh_type),                       INTENT(IN) :: H_mesh     !type de maillage
    INTEGER,                               INTENT(IN) :: m
    REAL(KIND=8), DIMENSION(H_mesh%np,6)              :: vv
    REAL(KIND=8), DIMENSION(:), POINTER               :: r, z
    REAL(KIND=8):: rho, phi
    INTEGER :: i, n
    REAL(Kind=8) :: pol_sur_tor,  PI, rr, zz, eps

    IF (test_conv_maxwell .AND. no_test==1) THEN


       IF (m==0) THEN
          vv = 0
          RETURN
       END IF
       r => H_mesh%rr(1,:)
       z => H_mesh%rr(2,:)
       vv(:,1) = alpha*z*(r**(m-1))*m     !-alpha*z*gamma/r**(m+1)*m
       vv(:,2) = beta *z*(r**(m-1))*m     !-beta *z*gamma/r**(m+1)*m
       vv(:,3) = beta *z*(r**(m-1))*m     !+beta *z*gamma/r**(m+1)*m
       vv(:,4) =-alpha*z*(r**(m-1))*m     !-alpha*z*gamma/r**(m+1)*m
       vv(:,5) = alpha*(r**m)             !+  alpha*gamma/r**m)
       vv(:,6) = beta *(r**m)             !+   beta*gamma/r**m)
       VV = vv/m**3
       RETURN
    ELSE
       vv  =0 
       RETURN
    END IF

    WRITE(*,*) ' BUG in Vexact_anal: Why did we arrive here?'
    STOP

  END FUNCTION Vexact_anal


  FUNCTION Hexact_anal(type, rr, m, mu_H_field, t) RESULT(vv) 

    USE bessel

    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: type
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m  !mode 
    REAL(KIND=8),                        INTENT(IN)   :: t 
    REAL(KIND=8), DIMENSION(:),          INTENT(IN)   :: mu_H_field
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: r, z 
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: f1, f2, f3, f4
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: df2, df3
    REAL(KIND=8)                                      :: f1_r0, f2_r0, f3_r0
    REAL(KIND=8)                                      :: df2_r0, df3_r0
    REAL(KIND=8)                                      :: A, B, z0, muH

    REAL(KIND=8)                                      :: aa=0.5d0, bb=1.d0, mu, mu0=1.d0, &
         capA,capB, capC, capD
    REAL(KIND=8),DIMENSION(SIZE(rr,2))                :: theta, rho
    INTEGER                                           :: n, ms, ns

    IF (test_conv_maxwell .AND. no_test==1) THEN
       IF (MAXVAL(mu_H_field) /= MINVAL(mu_H_field)) THEN
          WRITE(*,*) ' BUG in condlim, mu not constant'
          WRITE(*,*) MAXVAL(mu_H_field), MINVAL(mu_H_field)
          STOP
       END IF
       muH=mu_H_field(1) 
       r = rr(1,:)
       z = rr(2,:)
       IF (m==0) THEN
          vv = 0
          RETURN
       END IF
       IF (type == 1) THEN
          vv = alpha*z*(r**(m-1))*m !-alpha*z*gamma/r**(m+1)*m
       ELSEIF (type == 2) THEN
          vv = beta *z*(r**(m-1))*m !-beta *z*gamma/r**(m+1)*m
       ELSEIF (type ==3) THEN
          vv = beta *z*(r**(m-1))*m !+beta *z*gamma/r**(m+1)*m
       ELSEIF (type == 4)  THEN
          vv =-alpha*z*(r**(m-1))*m !+-alpha*z*gamma/r**(m+1)*m
       ELSEIF (type == 5) THEN
          vv = alpha*(r**m) ! +  alpha*(gamma/r**m)
       ELSEIF (type == 6) THEN
          vv = beta *(r**m) ! +  beta *(gamma/r**m)
       ENDIF

       vv = (vv/muH)*cos(t)/m**3
       RETURN
    END IF

    IF (test_conv_maxwell .AND. no_test==2) THEN
       r = rr(1,:)
       z = rr(2,:)
       IF (m/=m0) THEN
          vv = 0
          RETURN
       END IF
       !valeur en r0
       z0 = r0*k0
       f3_r0  = 1.d0
       df3_r0 = 2.d0/r0  
       f1_r0  = 1.d0/k0*(df3_r0-m0/r0*BESSK(m0,z0))
       !f2_r0  = m0/(k0*r0)*f3_r0 - (BESSI(m0-1,z0)-m0/(k0*r0)*BESSI(m0,z0))
       f2_r0  = m0/(k0*r0)*f3_r0 - (-BESSK(m0-1,z0)-m0/(k0*r0)*BESSK(m0,z0))
       df2_r0 = (m0/r0*f1_r0-f2_r0/r0-k0*BESSK(m0,z0))
       !pour evaluer f2
       A = 3*f2_r0 - r0*df2_r0
       B = r0*df2_r0-2*f2_r0    
       !fonction pour tous les points
       f1  = (r/r0)**2*f1_r0
       f2  = (r/r0)**2*(A+B*r/r0)
       f3  = (r/r0)**2
       df3 = 2*r/r0**2
       df2 = r/r0**2*(2*A+3*B*r/r0)
       f4  = r/r0**2*(A+B*r/r0) + df2 - m0*r/r0**2*f1_r0 

       IF (type == 1) THEN
          vv = (m0*r/r0**2-k0*f2)*COS(k0*z) 
       ELSEIF (type == 2) THEN
          vv = 0.d0
       ELSEIF (type ==3) THEN
          vv = 0.d0
       ELSEIF (type == 4)  THEN
          vv = (k0*f1-df3)*COS(k0*z)
       ELSEIF (type == 5) THEN
          vv = f4*SIN(k0*z)
       ELSEIF (type == 6) THEN
          vv = 0.d0
       ENDIF
       vv = vv * COS(t)
       RETURN
    END IF

    IF (test_conv_maxwell .AND. no_test==3) THEN
       IF (m/=0) THEN
          vv = 0.d0
          RETURN
       END IF
       IF (SIZE(rr,2)==0) RETURN
       IF (SIZE(rr,2)/=H_mesh_analytic%np) THEN
          WRITE(*,*) ' BUG in Hexact_anal '
          STOP
       END IF
       !CHAMP DURAND, H/H/phi configuration, H multi-valued at interface.
       ! mu1 = mu3 = mu_0 = 1
       ! mu2 = mu
       mu = MAXVAL(mu_H_field)
       r = H_mesh_analytic%rr(1,:)
       z = H_mesh_analytic%rr(2,:)
       theta = ATAN2(r,z)
       rho = SQRT(r**2+z**2)

       capA = -9*mu*mu0/((2*mu+mu0)*(mu+2*mu0) - 2*(mu-mu0)**2*(aa/bb)**3)
       capD = (2*mu+mu0)*(mu-mu0)*((bb/aa)**3-1.d0)/((2*mu+mu0)*(mu+2*mu0) - 2*(mu-mu0)**2*(aa/bb)**3)
       capC = (1.d0 - mu0/mu)*capA/3
       capB = (2.d0 + mu0/mu)*capA/3

       DO n = 1, H_mesh_analytic%np
          IF (rho(n) .LE. aa) THEN
             IF (type==1) THEN
                vv(n) =0.d0
             ELSE IF(type==5) THEN
                vv(n) = -capA
             ELSE
                vv(n) = 0.d0
             END IF
          ELSE
             IF (type==1) THEN
                vv(n) = 3*capC*(aa/rho(n))**3*COS(theta(n))*SIN(theta(n))
             ELSE IF(type==5) THEN
                vv(n) = -capB + capC*(aa/rho(n))**3*(3.d0*COS(theta(n))**2-1.d0)
             ELSE
                vv(n) = 0.d0
             END IF
          END IF
       END DO

       DO ms = 1, H_mesh_analytic%mes !Do loop on interface, since H multi-valued
          DO ns = 1, H_mesh_analytic%gauss%n_ws
             n = H_mesh_analytic%jjs(ns,ms)
             IF (H_mesh_analytic%i_d(H_mesh_analytic%neighs(ms)) == 1) THEN 
                IF (type==1) THEN
                   vv(n) = 0.d0
                ELSE IF(type==5) THEN
                   vv(n) = -capA
                ELSE
                   vv(n) = 0.d0
                END IF
             ELSE
                IF (type==1) THEN
                   vv(n) = 3*capC*(aa/rho(n))**3*COS(theta(n))*SIN(theta(n))
                ELSE IF(type==5) THEN
                   vv(n) = -capB + capC*(aa/rho(n))**3*(3.d0*COS(theta(n))**2-1.d0)
                ELSE
                   vv(n) = 0.d0
                END IF
             END IF
          END DO
       END DO
       RETURN
    END IF

    IF (test_conv_maxwell .AND. no_test==4) THEN
       IF (m/=0) THEN
          vv = 0.d0
          RETURN
       END IF

       !CHAMP DURAND, configuration phi/H/phi (H single-valued)
       mu = MAXVAL(mu_H_field)
       r = rr(1,:)
       z = rr(2,:)
       theta = ATAN2(r,z)
       rho = SQRT(r**2+z**2)

       capA = -9*mu*mu0/((2*mu+mu0)*(mu+2*mu0) - 2*(mu-mu0)**2*(aa/bb)**3)
       capD = (2*mu+mu0)*(mu-mu0)*((bb/aa)**3-1.d0)/((2*mu+mu0)*(mu+2*mu0) - 2*(mu-mu0)**2*(aa/bb)**3)
       capC = (1.d0 - mu0/mu)*capA/3
       capB = (2.d0 + mu0/mu)*capA/3

       IF (type==1) THEN
          vv = 3*capC*(aa/rho)**3*COS(theta)*SIN(theta) 
       ELSE IF(type==5) THEN
          vv = -capB + capC*(aa/rho)**3*(3.d0*COS(theta)**2-1.d0)
       ELSE
          vv = 0.d0
       END IF
       RETURN
    END IF

    WRITE(*,*) ' BUG in Hexact_anal: Why did we arrive here?'
    STOP

  END FUNCTION Hexact_anal

  FUNCTION Phiexact_anal(type, rr, m, mu_phi,t) RESULT(vv) 
    USE bessel
    IMPLICIT NONE

    INTEGER     ,                        INTENT(IN)   :: type
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m  !mode 
    REAL(KIND=8),                        INTENT(IN)   :: mu_phi, t
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv   
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: r, z, theta, rho
    INTEGER                                           :: n
    REAL(KIND=8)                                      :: a=0.5d0, b=1.d0, mu, mu0=1.d0, &
         capA,capB, capC, capD, t_al=5.d-2

    IF (test_conv_maxwell .AND. no_test==1) THEN
       r = rr(1,:)
       z = rr(2,:)
       vv(:) = 0.d0
       IF (m==0) RETURN
       IF  (type == 1) THEN
          vv = alpha*z*(r**m) ! + alpha*z*(gamma/r**m)
       ELSEIF (type == 2) THEN
          vv = beta *z*(r**m) ! + beta *z*(gamma/r**m)
       ENDIF

       vv = (vv/mu_phi)*cos(t)/m**3
       RETURN
    END IF

    IF (test_conv_maxwell .AND. no_test==2) THEN
       r = rr(1,:)
       z = rr(2,:)
       vv(:) = 0.d0
       IF (m/=m0) RETURN
       IF  (type == 1) THEN
          DO n= 1, SIZE(rr,2)
             vv(n) = BESSK(m0,k0*r(n))*COS(k0*z(n))
          END DO
       ELSEIF (type == 2) THEN
          vv = 0.d0
       ENDIF
       vv = vv*COS(t)
       RETURN
    END IF

    IF (test_conv_maxwell .AND. (no_test==3 .OR. no_test==4)) THEN
       IF (m/=0) THEN
          vv = 0.d0
          RETURN
       END IF

       !CHAMP DURAND
       !TESTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
       mu = 200.d0
       !TESTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
       r = rr(1,:)
       z = rr(2,:)   
       rho = SQRT(r**2+z**2)

       DO n = 1, SIZE(rho)
          IF (rho(n).LE.1.d-10) THEN
             theta(n) = 0.d0
          ELSE
             theta(n) = ATAN2(r(n),z(n))
          END IF
       END DO

       capA = -9*mu*mu0/((2*mu+mu0)*(mu+2*mu0) - 2*(mu-mu0)**2*(a/b)**3)
       capD = (2*mu+mu0)*(mu-mu0)*((b/a)**3-1.d0)/((2*mu+mu0)*(mu+2*mu0) - 2*(mu-mu0)**2*(a/b)**3)


       DO n = 1, SIZE(rho)
          IF (type==1 .AND. rho(n).LE. (a+1.d-1)) THEN
             vv(n) =  -capA*rho(n)*COS(theta(n)) 
          ELSE IF (type==1 .AND. rho(n) .GE. (b-1.d-1)) THEN
             vv(n) = (rho(n)*COS(theta(n)) - capD*COS(theta(n))*a**3/rho(n)**2) !*(t/t_al)**3/(1.d0+(t/t_al)**3)
          ELSE
             vv(n) = 0.d0
          END IF
       END DO

       RETURN
    END IF

    WRITE(*,*) ' BUG in Phiexact_anal: Why did we arrive here?'
    STOP

  END FUNCTION Phiexact_anal


  FUNCTION Jexact_gauss_anal(type, rr, m, mu_phi, sigma, mu_H, t) RESULT(vv) 

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

    IF (test_conv_maxwell .AND. no_test==1) THEN
       vv = -sigma* Eexact_gauss_anal(type, rr, m, mu_phi, sigma, mu_H, t)
       RETURN
    END IF
    
    IF (test_conv_maxwell .AND. no_test==2) THEN
      r = rr(1)
      z = rr(2)
      IF (m/=m0) THEN
         vv = 0
         RETURN
      END IF
      !calcul du rotationnel de H
      !valeur en r0
      z0 = r0*k0
      f3_r0  = 1.d0
      df3_r0 = 2.d0/r0
      f1_r0  = 1.d0/k0*(df3_r0-m0/r0*BESSK(m0,z0))
      !f2_r0  = m0/(k0*r0)*f3_r0 - (BESSI(m0-1,z0)-m0/(k0*r0)*BESSI(m0,z0))
      f2_r0  = m0/(k0*r0)*f3_r0 - (-BESSK(m0-1,z0)-m0/(k0*r0)*BESSK(m0,z0))
      df2_r0 = (m0/r0*f1_r0-f2_r0/r0-k0*BESSK(m0,z0))
      !pour evaluer f2
      A = 3*f2_r0 - r0*df2_r0
      B = r0*df2_r0-2*f2_r0
      !fonction pour tous les points
      f1  = (r/r0)**2*f1_r0
      f2  = (r/r0)**2*(A+B*r/r0)
      f3  = (r/r0)**2
      df3 = 2*r/r0**2
      df2 = r/r0**2*(2*A+3*B*r/r0)
      f4  = r/r0**2*(A+B*r/r0) + df2 - m0*r/r0**2*f1_r0
      df1 = 2*r/r0**2*f1_r0
      d2f3= 2/r0**2
      d2f2= (1/r0**2)*(2*A+6*B*r/r0)
      df4 = (1/r0**2)*(A+2*B*r/r0) + d2f2 -m0/r0**2*f1_r0

      !champ mag en r
      H1 = (m0*r/r0**2-k0*f2)
      H2 = (k0*f1-df3)
      H3 = f4
      dH2= k0*df1 - d2f3
      dH3= df4

      IF  (type == 1) THEN
         vv = 0.d0
      ELSEIF (type == 2) THEN
         vv = ((-m0/r)*H3+k0*H2)*SIN(k0*z)
      ELSEIF (type ==3) THEN
         vv = (-k0*H1-dH3)*SIN(k0*z)
      ELSEIF (type == 4)  THEN
         vv = 0.d0
      ELSEIF (type == 5) THEN
         vv = 0.d0
      ELSEIF (type == 6) THEN
         vv = 1/r*(H2 + r*dH2 + m0*H1)*COS(k0*z)
      ENDIF
      vv = vv*COS(t) - sigma* Eexact_gauss_anal(type, rr, m, mu_phi, sigma, mu_H, t)
     ! vv =  - sigma* Eexact_gauss(type, rr, m, mu_phi, sigma, mu_H, t)
       
      RETURN
    END IF

    IF (test_conv_maxwell .AND. (no_test==3 .OR. no_test==4)) THEN
       vv = 0.d0
       RETURN
    END IF

   WRITE(*,*) ' BUG Jexact_gauss_anal: Why did we arrive here?'
   STOP

  END FUNCTION Jexact_gauss_anal

  FUNCTION Eexact_gauss_anal(type, rr, m, mu_phi, sigma, mu_H, t) RESULT(vv) 

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

    IF (test_conv_maxwell .AND. no_test==1) THEN
       r = rr(1)
       z = rr(2)
       vv = 0.d0
      
       IF (m == 0) RETURN
       IF  (type == 1) THEN
          vv = 0.d0
       ELSEIF (type == 2) THEN
          vv = 0.d0
       ELSEIF (type ==3) THEN
          vv = alpha*(-1.d0/(m+2)*r**(m+1)) ! + alpha*(1.d0/(m-2)*gamma/r**(m-1))
       ELSEIF (type == 4)  THEN
          vv = beta *(-1.d0/(m+2)*r**(m+1))  !+  beta *(1.d0/(m-2)*gamma/r**(m-1))
       ELSEIF (type == 5) THEN
          vv =  beta*z*(r**m) !   beta*z*(-gamma/r**m)
       ELSEIF (type == 6) THEN
          vv =-alpha*z*(r**m) ! -alpha*z*(-gamma/r**m)
       ENDIF

       vv = -vv*sin(t)/m**3
       RETURN
    END IF

   IF (test_conv_maxwell .AND. no_test==2) THEN
      r = rr(1)
      z = rr(2)
      IF (m/=m0) THEN
         vv = 0
         RETURN
      END IF
      !valeur en r0
      z0 = r0*k0
      f3_r0  = 1.d0
      df3_r0 = 2.d0/r0
      f1_r0  = 1.d0/k0*(df3_r0-m0/r0*BESSK(m0,z0))
      !f2_r0  = m0/(k0*r0)*f3_r0 - (BESSI(m0-1,z0)-m0/(k0*r0)*BESSI(m0,z0))
      f2_r0  = m0/(k0*r0)*f3_r0 - (-BESSK(m0-1,z0)-m0/(k0*r0)*BESSK(m0,z0))
      df2_r0 = (m0/r0*f1_r0-f2_r0/r0-k0*BESSK(m0,z0))
      !pour evaluer f2
      A = 3*f2_r0 - r0*df2_r0
      B = r0*df2_r0-2*f2_r0
      !fonction pour tous les points
      f1  = (r/r0)**2*f1_r0
      f2  = (r/r0)**2*(A+B*r/r0)
      f3  = (r/r0)**2

      IF  (type == 1) THEN
         vv = 0.d0
      ELSEIF (type == 2) THEN
         vv = f1*SIN(k0*z)
      ELSEIF (type ==3) THEN
         vv = f2*SIN(k0*z)
      ELSEIF (type == 4)  THEN
         vv = 0.d0
      ELSEIF (type == 5) THEN
         vv = 0.d0
      ELSEIF (type == 6) THEN
         vv = f3*COS(k0*z)
      ENDIF
      vv = vv*sin(t)/mu_H
      RETURN
   END IF

    IF (test_conv_maxwell .AND. (no_test==3 .OR. no_test==4)) THEN
       vv = 0.d0
       RETURN
    END IF


   WRITE(*,*) ' BUG Eexact_gauss_anal: Why did we arrive here?'
   STOP

  END FUNCTION Eexact_gauss_anal


  SUBROUTINE init_maxwell_analytique(H_mesh,phi_mesh,time,dt,mu_H_field,mu_phi,list_mode,Hn1,Hn,phin1,phin)

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
    INTEGER :: i, k
    

    IF (test_conv_maxwell .AND. (no_test==3 .OR. no_test==4)) THEN
       Hn1 = 0.d0
       Hn = 0.d0
       phin1 = 0.d0
       phin = 0.d0
       RETURN
    END IF

    time = -dt

    DO k=1,6
       DO i=1, SIZE(list_mode)
          Hn1(:,k,i) = Hexact_anal(k, H_mesh%rr, list_mode(i), mu_H_field, time)
          IF (k<3) THEN
             phin1(:,k,i) = Phiexact_anal(k, phi_mesh%rr, list_mode(i) , mu_phi, time)
          ENDIF
       ENDDO
    ENDDO
    
    time = time + dt 
    DO k=1,6
       DO i=1, SIZE(list_mode)
          Hn(:,k,i) = Hexact_anal(k, H_mesh%rr, list_mode(i), mu_H_field, time)
          IF (k<3) THEN
             phin(:,k,i) = Phiexact_anal(k, phi_mesh%rr, list_mode(i), mu_phi, time)
          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE init_maxwell_analytique

END MODULE boundary_anal
