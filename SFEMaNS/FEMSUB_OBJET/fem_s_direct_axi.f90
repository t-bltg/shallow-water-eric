!
!
!
!Authors: Jean-Luc Guermond, Raphael Laguerre, Copyrights 2000, 2004
!
MODULE  fem_s_axi

CONTAINS


  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  !Procedure pour le cas scalaire
  !----------------------------------------------------------------------
  !----------------------------------------------------------------------

  SUBROUTINE Dirichlet (jjs, us, ff)
    !=================================

    IMPLICIT NONE
    INTEGER, DIMENSION(:),        INTENT(IN)   :: jjs
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)   :: us
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT)  :: ff

    INTEGER :: n, i

    DO n = 1, SIZE(jjs);  i = jjs(n)
       ff(i) = us(n)
    END DO

  END SUBROUTINE Dirichlet


  SUBROUTINE Moy(mesh,p,RESULT)
    !===========================
    !moyenne

    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:)  ,   INTENT(IN)  :: p
    REAL(KIND=8)                ,   INTENT(OUT) :: RESULT
    REAL(KIND=8)                                :: vol   
    INTEGER ::  m, l , i , ni,k


    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8)   :: ray

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me


    RESULT = 0.d0
    vol = 0.d0   

    DO m = 1, me
       DO l = 1, l_G

          !--------On calcule le rayon du point gauss
          ray = 0
          DO ni = 1, n_w;  i = jj(ni,m)
             ray = ray + mesh%rr(1,i)*ww(ni,l)
          END DO

          RESULT = RESULT +  SUM(p(jj(:,m)) * ww(:,l))* ray* rj(l,m)
          vol = vol + ray* rj(l,m)

       ENDDO
    ENDDO

    RESULT = RESULT / vol

  END SUBROUTINE Moy



  SUBROUTINE qs_00_ssr(mesh, ff,  u0)   !sans  le rayon dans l'integration
    !=================================

    !second membre simple 
    !  < w, f >   ===>   u0

    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
    REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0

    INTEGER ::  m, l , i , ni
    REAL(KIND=8) :: fl

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8)   :: ray

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    u0 = 0

    DO m = 1, me
       DO l = 1, l_G

          !--------On calcule le rayon du point gauss
          ray = 0
          DO ni = 1, n_w;  i = jj(ni,m)
             ray = ray + mesh%rr(1,i)*ww(ni,l)
          END DO

          fl = SUM(ff(jj(:,m)) * ww(:,l)) * rj(l,m)

          DO ni = 1, n_w 
             u0(jj(ni,m)) = u0(jj(ni,m)) +  ww(ni,l) * fl
          ENDDO

       ENDDO
    ENDDO

  END SUBROUTINE qs_00_ssr

  SUBROUTINE qs_00(mesh, ff,  u0)   !avec  le rayon dans l'integration
    !=================================
    !second membre simple
    !  < w, f >   ===>   u0

    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
    REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0

    INTEGER ::  m, l , i , ni
    REAL(KIND=8) :: fl

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8)   :: ray

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    u0 = 0

    DO m = 1, me
       DO l = 1, l_G

          !--------On calcule le rayon du point gauss
          ray = 0
          DO ni = 1, n_w;  i = jj(ni,m)
             ray = ray + mesh%rr(1,i)*ww(ni,l)
          END DO

          fl = SUM(ff(jj(:,m)) * ww(:,l)) * rj(l,m)

          DO ni = 1, n_w 
             u0(jj(ni,m)) = u0(jj(ni,m)) +  ww(ni,l) * fl *ray
          ENDDO

       ENDDO
    ENDDO

  END SUBROUTINE qs_00


  SUBROUTINE qs_00_inst_3d_ssr(mesh,ff,T1,T2,dt,u0) !sans le rayon
    !=================================
    !second membre pour l'equation de la diffusion avec BDF2
    !  < w, f > + <w,(4T(n)-T(n-1))/2dt>     ===>   u0 pour les cos(1) et sin(2)

    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff   !analytique (noeud)
    REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: T1   !T(noeud) ,Tn
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: T2   !T(noeud) ,Tn-1
    REAL(KIND=8),                 INTENT(IN)  :: dt
    INTEGER ::  m, l , i , ni
    REAL(KIND=8) :: fl   


    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8)   :: ray

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    u0 = 0

    DO m = 1, me
       DO l = 1, l_G

          !--------On calcule le rayon du point gauss
          ray = 0
          DO ni = 1, n_w;  i = jj(ni,m)
             ray = ray + mesh%rr(1,i)*ww(ni,l)
          END DO

          fl = SUM(ff(jj(:,m)) * ww(:,l)) * rj(l,m) &
               + SUM(4 * T1(jj(:,m)) * ww(:,l) - T2(jj(:,m))*ww(:,l))*rj(l,m) /(2* dt)*ray

          DO ni = 1, n_w 
             u0(jj(ni,m)) = u0(jj(ni,m)) +  ww(ni,l)  * fl 
          ENDDO

       ENDDO
    ENDDO

  END SUBROUTINE qs_00_inst_3d_ssr

  SUBROUTINE qs_00_inst_3d(mesh,ff,T1,T2,dt,u0) !avec  le rayon
    !=================================

    !second membre pour l'equation de la diffusion avec BDF2
    !  < w, f > + <w,(4T(n)-T(n-1))/2dt>     ===>   u0 pour les cos(1) et sin(2)

    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff   !analytique (noeud)
    REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: T1   !T(noeud) ,Tn
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: T2   !T(noeud) ,Tn-1
    REAL(KIND=8),                 INTENT(IN)  :: dt
    INTEGER ::  m, l , i , ni
    REAL(KIND=8) :: fl   


    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8)   :: ray

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    u0 = 0

    DO m = 1, me
       DO l = 1, l_G

          !--------On calcul le rayon du point gauss
          ray = 0
          DO ni = 1, n_w;  i = jj(ni,m)
             ray = ray + mesh%rr(1,i)*ww(ni,l)
          END DO

          fl = SUM(ff(jj(:,m)) * ww(:,l)) * rj(l,m)  &
               + SUM(4 * T1(jj(:,m)) * ww(:,l) - T2(jj(:,m))*ww(:,l))*rj(l,m) /(2* dt)

          DO ni = 1, n_w 
             u0(jj(ni,m)) = u0(jj(ni,m)) +  ww(ni,l)  * fl *ray
          ENDDO

       ENDDO
    ENDDO

  END SUBROUTINE qs_00_inst_3d


  SUBROUTINE qs_00_inst_init_3d(mesh,ff,T1,dt,u0)  !avec le rayon
    !=================================
    !second membre pour l'equation de la diffusion avec Euler retarde
    !  < w, f > +<w,T(n)/dt>     ===>   u0 pour les cos(1) et sin(2)

    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff   !analytique (noeud)
    REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: T1   !T(noeud) ,Tn
    REAL(KIND=8),                 INTENT(IN)  :: dt
    INTEGER ::  m, l , i , ni
    REAL(KIND=8) :: fl   


    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8)   :: ray

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    u0 = 0

    DO m = 1, me
       DO l = 1, l_G

          !--------On calcul le rayon du point gauss
          ray = 0
          DO ni = 1, n_w;  i = jj(ni,m)
             ray = ray + mesh%rr(1,i)*ww(ni,l)
          END DO

          fl = SUM(ff(jj(:,m)) * ww(:,l)) * rj(l,m) &
               +SUM(T1(jj(:,m)) * ww(:,l)) * rj(l,m) / dt


          DO ni = 1, n_w 
             u0(jj(ni,m)) = u0(jj(ni,m)) +  ww(ni,l) * ray * fl
          ENDDO

       ENDDO
    ENDDO

  END SUBROUTINE qs_00_inst_init_3d

  SUBROUTINE qs_00_inst_init_3d_ssr(mesh,ff,T1,dt,u0)  !sans le rayon
    !=================================
    !second membre pour l'equation de la diffusion avec Euler retarde
    !  < w, f > +<w,T(n)/dt>     ===>   u0 pour les cos(1) et sin(2)

    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff   !analytique (noeud)
    REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: T1   !T(noeud) ,Tn
    REAL(KIND=8),                 INTENT(IN)  :: dt
    INTEGER ::  m, l , i , ni
    REAL(KIND=8) :: fl   


    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8)   :: ray

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    u0 = 0

    DO m = 1, me
       DO l = 1, l_G

          !--------On calcul le rayon du point gauss
          ray = 0.d0
          DO ni = 1, n_w;  i = jj(ni,m)
             ray = ray + mesh%rr(1,i)*ww(ni,l)
          END DO

          fl = SUM(ff(jj(:,m)) * ww(:,l)) * rj(l,m) &
               +SUM(T1(jj(:,m)) * ww(:,l)) * rj(l,m) / dt  *ray


          DO ni = 1, n_w 
             u0(jj(ni,m)) = u0(jj(ni,m)) +  ww(ni,l)  * fl
          ENDDO

       ENDDO
    ENDDO

  END SUBROUTINE qs_00_inst_init_3d_ssr


  SUBROUTINE qs_00_adv_diff_3d(eps,mode,mesh,ff,T1,T2,gg,dt,u0)
    !=================================
    !second membre pour l'advection diffusion d'un scalaire avec un BDF2
    !  < w, f > + <w,(4T(n)-T(n-1))/2dt> + <w,eps* mV_theta /r (2T(n)-T(n-1)> 
    !                                         ===>   u0 pour les cos(1) et sin(2)
    !  le +ou- est donne par eps
    !eps = 1 si on traite les cos et -1 sinon


    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff   !analytique (noeud)
    REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: T1   !T(noeud) ,Tn
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: T2   !T(noeud) ,Tn-1
    REAL(KIND=8),                 INTENT(IN)  :: dt
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: gg   !champ de vitesse V_theta
    REAL(KIND=8),                 INTENT(IN)  :: eps
    INTEGER     ,                 INTENT(IN)  :: mode
    INTEGER ::  m, l , i , ni
    REAL(KIND=8) :: fs , ft , fexp   


    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8)   :: ray

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    u0 = 0

    DO m = 1, me
       DO l = 1, l_G

          !--------On calcul le rayon du point gauss
          ray = 0
          DO ni = 1, n_w;  i = jj(ni,m)
             ray = ray + mesh%rr(1,i)*ww(ni,l)
          END DO

          fs = SUM(ff(jj(:,m)) * ww(:,l)) * rj(l,m) 
          ft = SUM(4 * T1(jj(:,m)) * ww(:,l) - T2(jj(:,m))*ww(:,l))*rj(l,m) /(2* dt) 
          fexp = eps * mode * SUM(gg(3,jj(:,m)) * T1(jj(:,m)) * ww(:,l))*rj(l,m) / ray


          DO ni = 1, n_w 
             u0(jj(ni,m)) = u0(jj(ni,m)) +  ww(ni,l) * ray *( fs + ft + fexp)  
          ENDDO

       ENDDO
    ENDDO

  END SUBROUTINE qs_00_adv_diff_3d

  SUBROUTINE qs_00_adv_diff_3d_ssr(eps,mode,mesh,ff,T1,T2,gg,dt,u0) !sans le rayon dans l'integral
    !=================================
    !second membre pour l'advection diffusion d'un scalaire avec un BDF2
    !  (< w, f > + <w,(4T(n)-T(n-1))/2dt> + <w,esp* mV_theta /r (2T(n)-T(n-1)> )
    !                                         ===>   u0 pour les cos(1) et sin(2)
    !  le +ou- est donne par eps
    !eps = 1 si on traite les cos et -1 sinon


    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff   !analytique (noeud)
    REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: T1   !T(type,noeud) ,Tn
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: T2   !T(type,noeud) ,Tn-1
    REAL(KIND=8),                 INTENT(IN)  :: dt
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: gg   !champ de vitesse V_theta
    REAL(KIND=8),                 INTENT(IN)  :: eps
    INTEGER     ,                 INTENT(IN)  :: mode
    INTEGER ::  m, l , i , ni
    REAL(KIND=8) :: fs , ft , fexp   


    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8)   :: ray

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    u0 = 0

    DO m = 1, me
       DO l = 1, l_G

          !--------On calcul le rayon du point gauss
          ray = 0
          DO ni = 1, n_w;  i = jj(ni,m)
             ray = ray + mesh%rr(1,i)*ww(ni,l)
          END DO

          fs = SUM(ff(jj(:,m)) * ww(:,l)) * rj(l,m) 

          IF (eps.EQ.(-1.d0)) THEN !cas des cosinus 
             ft = SUM(4 * T1(1,jj(:,m)) * ww(:,l) - T2(1,jj(:,m))*ww(:,l))*rj(l,m) /(2* dt)
             fexp = eps * mode * SUM(gg(3,jj(:,m)) * (2*T1(2,jj(:,m))-T2(2,jj(:,m))) &
                  * ww(:,l))*rj(l,m)/ray
          ELSEIF  (eps.EQ.(1.d0)) THEN  !cas des sinus
             ft = SUM(4 * T1(2,jj(:,m)) * ww(:,l) - T2(2,jj(:,m))*ww(:,l))*rj(l,m) /(2* dt)
             fexp = eps * mode * SUM(gg(3,jj(:,m)) * (2*T1(1,jj(:,m))-T2(1,jj(:,m))) & 
                  * ww(:,l))*rj(l,m)/ray
          ELSE
             WRITE(*,*)  'probleme ds calcul second membre avec epsilon'
          ENDIF

          DO ni = 1, n_w 
             u0(jj(ni,m)) = u0(jj(ni,m)) +  ww(ni,l) * (fs + ray * (ft + fexp))  
          ENDDO

       ENDDO
    ENDDO

  END SUBROUTINE qs_00_adv_diff_3d_ssr



  SUBROUTINE qs_00_adv_diff_init_3d(eps,mode,mesh,ff,T1,dt,gg,u0)
    !=================================
    !second membre pour l'advection diffusion d'un scalaire avec un Euler retarde
    !  < w, f > +<w,T(n)/dt> + <w,eps* mV_theta /r T(n)>     
    !      ===>   u0 pour les cos(1) et sin(2)

    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff   !analytique (noeud)
    REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: T1   !T(noeud) ,Tn
    REAL(KIND=8),                 INTENT(IN)  :: dt 
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN):: gg   !champ de vitesse V_theta
    REAL(KIND=8),                 INTENT(IN)  :: eps
    INTEGER     ,                 INTENT(IN)  :: mode
    INTEGER ::  m, l , i , ni
    REAL(KIND=8) :: fs , ft , fexp   


    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8)   :: ray

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    u0 = 0

    DO m = 1, me
       DO l = 1, l_G

          !--------On calcul le rayon du point gauss
          ray = 0
          DO ni = 1, n_w;  i = jj(ni,m)
             ray = ray + mesh%rr(1,i)*ww(ni,l)
          END DO

          fs = SUM(ff(jj(:,m)) * ww(:,l)) * rj(l,m) 
          ft = SUM(T1(jj(:,m)) * ww(:,l)) * rj(l,m) / dt
          fexp = eps * mode * SUM(gg(3,jj(:,m)) * T1(jj(:,m)) * ww(:,l))*rj(l,m)/ray

          DO ni = 1, n_w 
             u0(jj(ni,m)) = u0(jj(ni,m)) +  ww(ni,l) * ray * (fs + ft + fexp)
          ENDDO

       ENDDO
    ENDDO

  END SUBROUTINE qs_00_adv_diff_init_3d

  SUBROUTINE qs_00_adv_diff_init_3d_ssr(eps,mode,mesh,ff,T1,dt,gg,u0)   !sans le rayon ds l'integrale
    !=================================
    !second membre pour l'advection diffusion d'un scalaire avec un Euler retarde
    !  < w, f > +<w,T(n)/dt> + <w,eps* mV_theta /r T(n)>     
    !      ===>   u0 pour les cos(1) et sin(2)

    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff   !analytique (noeud)
    REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: T1   !T(noeud) ,Tn
    REAL(KIND=8),                 INTENT(IN)  :: dt 
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg   !champ de vitesse V_theta
    REAL(KIND=8),                 INTENT(IN)  :: eps
    INTEGER     ,                 INTENT(IN)  :: mode
    INTEGER ::  m, l , i , ni
    REAL(KIND=8) :: fs , ft , fexp   


    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8)   :: ray

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    u0 = 0

    DO m = 1, me
       DO l = 1, l_G

          !--------On calcul le rayon du point gauss
          ray = 0
          DO ni = 1, n_w;  i = jj(ni,m)
             ray = ray + mesh%rr(1,i)*ww(ni,l)
          END DO

          IF (eps.EQ.(-1.d0)) THEN    !cas des cosinus
             fs = SUM(ff(jj(:,m)) * ww(:,l)) * rj(l,m) 
             ft = SUM(T1(1,jj(:,m)) * ww(:,l)) * rj(l,m) / dt
             fexp = eps * mode * SUM(gg(3,jj(:,m)) * T1(2,jj(:,m)) * ww(:,l))*rj(l,m)/ray 
          ELSEIF  (eps.EQ.(1.d0)) THEN    !cas des sinus
             fs = SUM(ff(jj(:,m)) * ww(:,l)) * rj(l,m) 
             ft = SUM(T1(2,jj(:,m)) * ww(:,l)) * rj(l,m) / dt
             fexp = eps * mode * SUM(gg(3,jj(:,m)) * T1(1,jj(:,m))* ww(:,l))*rj(l,m)/ray 
          ELSE
             WRITE(*,*) 'problem calc second membre init avec epsilon'
          ENDIF

          DO ni = 1, n_w 
             u0(jj(ni,m)) = u0(jj(ni,m)) +  ww(ni,l) * (fs + ray *(ft + fexp))
          ENDDO

       ENDDO
    ENDDO

  END SUBROUTINE qs_00_adv_diff_init_3d_ssr


  SUBROUTINE qs_00_lap(mesh,alpha,mode,ff,V2,V1,dt,u0) !avec  le rayon
    !=================================================
    !calcul du second membre en entier pour le laplacien scalaire 
    !et un BDF2
    !<w,f> + <w,(4V(n)-V(n-1))/2dt> ==> u0


    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: ff   !analytique (type,noeud)
    REAL(KIND=8), DIMENSION(:,:),   INTENT(OUT) :: u0   !u0(type,noeud)
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: V2   !V(type,noeud) ,Vn
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: V1   !V(type,noeud) ,V(n-1)
    REAL(KIND=8),                   INTENT(IN)  :: dt , alpha
    INTEGER                     ,   INTENT(IN)  :: mode
    INTEGER ::  m, l , i , ni , j
    REAL(KIND=8), DIMENSION(6) :: fs  , ft

    !fs pour le terme de forcage
    !ft pour le terme temporelle   

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8)   :: ray

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    u0 = 0

    DO m = 1, me
       DO l = 1, l_G

          !--------On calcul le rayon du point gauss
          ray = 0
          DO ni = 1, n_w;  i = jj(ni,m)
             ray = ray + mesh%rr(1,i)*ww(ni,l)
          END DO

          !--------calcul de la premiere composante : u0(1,:) <--> f(r,m,c)

          fs(1) = SUM(ff(1,jj(:,m)) * ww(:,l)) * rj(l,m)  
          ft(1) = 1/(2*dt) * SUM((4*V2(1,jj(:,m))-V1(1,jj(:,m))) * ww(:,l)) * rj(l,m)

          !--------calcul de la seconde composante : u0(2,:) <--> f(r,m,s)

          fs(2) = SUM(ff(2,jj(:,m)) * ww(:,l)) * rj(l,m)
          ft(2) = 1/(2*dt)* SUM((4*V2(2,jj(:,m))-V1(2,jj(:,m))) * ww(:,l)) * rj(l,m)

          !---------- calcul du second membre total

          DO j=5,6
             DO ni = 1, n_w 
                u0(j,jj(ni,m)) = u0(j,jj(ni,m)) +  ww(ni,l)  * ray * (fs(j) + ft(j)) 
             ENDDO
          ENDDO


       ENDDO
    ENDDO

  END SUBROUTINE qs_00_lap

  SUBROUTINE qs_00_lap_init(mesh,alpha,mode,ff,V,dt,u0) !avec  le rayon
    !=================================================
    !calcul du second membre en entier pour le laplacien scalaire 
    !et un Euler retarde
    !<w,f> + <w,(V(n))/dt> ==> u0


    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: ff   !analytique (type,noeud)
    REAL(KIND=8), DIMENSION(:,:),   INTENT(OUT) :: u0   !u0(type,noeud)
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: V   !V(type,noeud) ,Vn
    REAL(KIND=8),                   INTENT(IN)  :: dt , alpha
    INTEGER                     ,   INTENT(IN)  :: mode
    INTEGER ::  m, l , i , ni , j
    REAL(KIND=8), DIMENSION(6) :: fs  , ft

    !fs pour le terme de forcage
    !ft pour le terme temporelle   

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8)   :: ray

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    u0 = 0

    DO m = 1, me
       DO l = 1, l_G

          !--------On calcul le rayon du point gauss
          ray = 0
          DO ni = 1, n_w;  i = jj(ni,m)
             ray = ray + mesh%rr(1,i)*ww(ni,l)
          END DO

          !--------calcul de la premiere composante : u0(1,:) <--> f(r,m,c)

          fs(1) = SUM(ff(2,jj(:,m)) * ww(:,l)) * rj(l,m)  
          ft(1) = 1/(dt) * SUM(V(2,jj(:,m)) * ww(:,l)) * rj(l,m)

          !--------calcul de la seconde composante : u0(2,:) <--> f(r,m,s)

          fs(2) = SUM(ff(2,jj(:,m)) * ww(:,l)) * rj(l,m)
          ft(2) = 1/(dt) * SUM(V(2,jj(:,m)) * ww(:,l)) * rj(l,m)

          !---------- calcul du second membre total

          DO j=5,6
             DO ni = 1, n_w 
                u0(j,jj(ni,m)) = u0(j,jj(ni,m)) +  ww(ni,l)  * ray * (fs(j)  + ft(j)) 
             ENDDO
          ENDDO


       ENDDO
    ENDDO

  END SUBROUTINE qs_00_lap_init

  SUBROUTINE qs_00_rot(mesh, mode, V, Rot)

    USE Gauss_points
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:),   INTENT(OUT) :: Rot !Rot(type,noeud)
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: V   !V(type,noeud) ,Vn
    INTEGER                     ,   INTENT(IN)  :: mode
    INTEGER ::  m, l , i , ni , j
    REAL(KIND=8), DIMENSION(6) :: f
    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8)   :: ray, rayj
    REAL(KIND=8), DIMENSION(6,mesh%gauss%n_w) :: V_loc
    REAL(KIND=8), DIMENSION(2,mesh%gauss%n_w) :: dw_loc
    INTEGER,      DIMENSION(mesh%gauss%n_w) :: j_loc
    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    Rot = 0.d0

    DO m = 1, me
       j_loc = jj(:,m)
       V_loc = V(:,j_loc)
       DO l = 1, l_G
          dw_loc = dw(:,:,l,m)

          !--------On calcule le rayon du point gauss
          ray = 0
          DO ni = 1, n_w;  i = jj(ni,m)
             ray = ray + mesh%rr(1,i)*ww(ni,l)
          END DO
          rayj = ray*rj(l,m)
          !----------------------

          f(1) = ( mode/ray*SUM(V_loc(6,:)*ww(:,l))-SUM(V_loc(3,:)*dw_loc(2,:)))*rayj
          f(2) = (-mode/ray*SUM(V_loc(5,:)*ww(:,l))-SUM(V_loc(4,:)*dw_loc(2,:)))*rayj
          f(3) = (SUM(V_loc(1,:)*dw_loc(2,:))-SUM(V_loc(5,:)*dw_loc(1,:)))*rayj
          f(4) = (SUM(V_loc(2,:)*dw_loc(2,:))-SUM(V_loc(6,:)*dw_loc(1,:)))*rayj
          f(5) =   (1/ray*SUM(V_loc(3,:)*ww(:,l))+SUM(V_loc(3,:)*dw_loc(1,:))&
               -mode/ray*SUM(V_loc(2,:)*ww(:,l)))*rayj
          f(6) =   (1/ray*SUM(V_loc(4,:)*ww(:,l))+SUM(V_loc(4,:)*dw_loc(1,:))&
               +mode/ray*SUM(V_loc(1,:)*ww(:,l)))*rayj

          DO ni = 1, n_w 
             DO j= 1,6
                Rot(j,j_loc(ni)) = Rot(j,j_loc(ni)) +  ww(ni,l)*f(j)                 
             ENDDO
          ENDDO

       ENDDO
    ENDDO

  END SUBROUTINE qs_00_rot



  SUBROUTINE qs_00_div(mesh,mode,V,u0) 
    !===============================================
    !calucl de la divergence du champ V(6 composantes : cos,sin,cos,sin ...)
    !en spectral pour theta
    !<w,div(V)> ==> u0
    USE Gauss_points

    IMPLICIT NONE


    REAL(KIND=8), DIMENSION(:,:),   INTENT(OUT) :: u0   !u0(type,noeud)
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: V   !V(type,noeud) ,Vn
    INTEGER                     ,   INTENT(IN)  :: mode
    INTEGER ::  m, l , i , ni , j
    REAL(KIND=8), DIMENSION(2) :: fr  , fth , fz

    !fr pour le terme de forcage
    !fth pour le terme temporelle   
    !fz

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8)   :: ray

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    u0 = 0
    fr = 0.d0
    fth = 0.d0
    fz = 0.d0

    DO m = 1, me
       DO l = 1, l_G

          !--------On calcul le rayon du point gauss
          ray = 0
          DO ni = 1, n_w;  i = jj(ni,m)
             ray = ray + mesh%rr(1,i)*ww(ni,l)
          END DO
          !----------------------

          !fr(1) = -SUM(V(1,jj(:,m))*ww(:,l))* rj(l,m)
          !fth(1) = mode*SUM(V(4,jj(:,m))*ww(:,l))* rj(l,m)
          !fz(1) = -SUM(V(5,jj(:,m))*ww(:,l))* rj(l,m)

          !fr(2) = -SUM(V(2,jj(:,m))*ww(:,l))* rj(l,m)
          !fth(2) = - mode*SUM(V(3,jj(:,m))*ww(:,l))* rj(l,m)
          !fz(2) = -SUM(V(6,jj(:,m))*ww(:,l))* rj(l,m)

          fr(1) = -SUM(V(1,jj(:,m))*dw(1,:,l,m))* rj(l,m)
          fth(1) = mode*SUM(V(4,jj(:,m))*ww(:,l))* rj(l,m)/ray
          fz(1) = -SUM(V(5,jj(:,m))*dw(2,:,l,m))* rj(l,m)

          fr(2) = -SUM(V(2,jj(:,m))*dw(1,:,l,m))* rj(l,m)
          fth(2) = - mode*SUM(V(3,jj(:,m))*ww(:,l))* rj(l,m)/ray
          fz(2) = -SUM(V(6,jj(:,m))*dw(2,:,l,m))* rj(l,m)


          !---------- calcul du second membre total

          DO j=1,2
             DO ni = 1, n_w 
                !u0(j,jj(ni,m)) = u0(j,jj(ni,m)) +  ray*(ww(ni,l)*fth(j) + dw(1,ni,l,m)*fr(j) &
                !                 +dw(2,ni,l,m)*fz(j))
                u0(j,jj(ni,m)) = u0(j,jj(ni,m)) +  ray*ww(ni,l)*(fth(j) + fr(j) + fz(j))                 
             ENDDO
          ENDDO


       ENDDO
    ENDDO

  END SUBROUTINE qs_00_div



  !------------------------------------------------------------------------------
  SUBROUTINE qs_01_div_hybrid_old (uu_mesh, pp_mesh,mode, gg,  u0_c)
    !================================================

    !  < w_c, D.g >   ===>   u0_c

    USE Gauss_points

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: gg
    REAL(KIND=8), DIMENSION(:,:),   INTENT(OUT) :: u0_c
    INTEGER     ,                   INTENT(IN)  :: mode
    TYPE(mesh_type), TARGET                     :: uu_mesh, pp_mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj, jj_c
    INTEGER,      POINTER       :: me, nwc

    REAL(KIND=8), DIMENSION(pp_mesh%gauss%n_w, uu_mesh%gauss%l_G) ::  w_c
    INTEGER :: m, l, n, k,i,j
    REAL(KIND=8) :: dgl
    REAL(KIND=8), DIMENSION(6,uu_mesh%gauss%n_w) :: ggkn 
    REAL(KIND=8), DIMENSION(uu_mesh%gauss%k_d,uu_mesh%gauss%n_w) :: dwkn
    REAL(KIND=8), DIMENSION(2,uu_mesh%gauss%n_w) :: u0_cn
    REAL(KIND=8), DIMENSION(2) :: fr  , fth , fz
    REAL(KIND=8), DIMENSION(3,2) :: f   
    REAL(KIND=8)   :: ray
    INTEGER, DIMENSION(uu_mesh%gauss%n_w) :: j_loc   
    REAL(KIND=8), DIMENSION(2)   :: smb   
    REAL(KIND=8):: user_time
    !EXTERNAL :: user_time
    REAL(KIND=8) :: tps, dummy, tt,tps1     


    CALL gauss(uu_mesh)
    jj => uu_mesh%jj
    jj_c => pp_mesh%jj
    me => uu_mesh%me
    nwc => pp_mesh%gauss%n_w

    DO l = 1, l_G   
       w_c(1,l) = ww(1,l) + 0.5*(ww(n_w-1,l) + ww(n_w,l)) 
       w_c(2,l) = ww(2,l) + 0.5*(ww(n_w,l) + ww(4,l)) 
       w_c(3,l) = ww(3,l) + 0.5*(ww(4,l) + ww(n_w-1,l)) 
    END DO

    !write(*,*) 'me=',me
    !write(*,*) 'lg=',l_G
    !write(*,*) 'n_w=',n_w

    u0_c = 0.d0   
    tps = 0.d0

    DO m = 1, me
       j_loc(:) = jj(:,m)
       !u0_cn = 0.d0
       DO l = 1, l_G

          !--------On calcule le rayon du point gauss des P2
          ray = 0
          DO n = 1, n_w;  i = jj(n,m)
             ray = ray + uu_mesh%rr(1,i)*ww(n,l)
          END DO
          !----------------------

          !tt = user_time(dummy)
          !calcul de la divergence sur les cos
          f(1,1)  = (ray*SUM(gg(1,j_loc(:))*dw(1,:,l,m)) &
               + SUM(gg(1,j_loc(:))*ww(:,l))) 
          f(2,1) = mode*SUM(gg(4,j_loc(:))*ww(:,l)) 
          f(3,1)  = SUM(gg(5,j_loc(:))*dw(2,:,l,m)) * ray

          !calcul de la divergence sur les sin
          f(1,2)  = (ray*SUM(gg(2,j_loc(:))*dw(1,:,l,m)) &
               + SUM(gg(2,j_loc(:))*ww(:,l))) 
          f(2,2) = - mode*SUM(gg(3,j_loc(:))*ww(:,l))
          f(3,2)  = SUM(gg(6,j_loc(:))*dw(2,:,l,m)) * ray

          f = f *rj(l,m)



          !tps = tps +user_time(dummy) -tt

          DO j=1,2
             DO n=1, nwc
                u0_c(j,jj_c(n,m)) = u0_c(j,jj_c(n,m)) + w_c(n,l)*(f(1,j)+f(2,j)+f(3,j))
             ENDDO
          ENDDO

          !calcul de la divergence sur les cos
          !fr(1)  = (ray*SUM(gg(1,jj(:,m))*dw(1,:,l,m)) &
          !            + SUM(gg(1,jj(:,m))*ww(:,l))) 
          !fth(1) = mode*SUM(gg(4,jj(:,m))*ww(:,l)) 
          !fz(1)  = SUM(gg(5,jj(:,m))*dw(2,:,l,m)) * ray
          ! 
          ! !calcul de la divergence sur les sin
          ! fr(2)  = (ray*SUM(gg(2,jj(:,m))*dw(1,:,l,m)) &
          !             + SUM(gg(2,jj(:,m))*ww(:,l))) 
          ! fth(2) = - mode*SUM(gg(3,jj(:,m))*ww(:,l))
          ! fz(2)  = SUM(gg(6,jj(:,m))*dw(2,:,l,m)) * ray

          ! fr(:)  = fr(:)   * rj(l,m)  
          ! fth(:) = fth(:)  * rj(l,m)  
          ! fz(:)  = fz(:)   * rj(l,m) 
          !integration sur les P1
          !j=1 pour les cos et j=2 pour les sin

          !DO n = 1, nwc 
          !   DO j = 1, 2
          !      u0_cn(j,n) = u0_cn(j,n) + w_c(n,l)*(fr(j) + fz(j) + fth(j))
          !   ENDDO
          !END DO

       ENDDO

       !DO n = 1, nwc 
       !   DO j=1, 2
       !      u0_c(j,jj_c(n,m)) = u0_c(j,jj_c(n,m)) + u0_cn(j,n)
       !   ENDDO
       !END DO


    ENDDO
    !write(*,*) 'tps=',tps 

  END SUBROUTINE qs_01_div_hybrid_old


  !------------------------------------------------------------------------------
  SUBROUTINE qs_01_div_hybrid (uu_mesh, pp_mesh,mode, gg,  u0_c)
    !================================================

    !  < w_c, D.g >   ===>   u0_c

    USE Gauss_points

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: gg
    REAL(KIND=8), DIMENSION(:,:),   INTENT(OUT) :: u0_c
    INTEGER     ,                   INTENT(IN)  :: mode
    TYPE(mesh_type), TARGET                     :: uu_mesh, pp_mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj, jj_c
    INTEGER,      POINTER       :: me, nwc

    REAL(KIND=8), DIMENSION(pp_mesh%gauss%n_w, uu_mesh%gauss%l_G) ::  w_c
    INTEGER :: m, l, n, k,i,j
    REAL(KIND=8) :: dgl
    REAL(KIND=8), DIMENSION(6,uu_mesh%gauss%n_w) :: ggkn 
    REAL(KIND=8), DIMENSION(uu_mesh%gauss%k_d,uu_mesh%gauss%n_w) :: dwkn
    REAL(KIND=8), DIMENSION(2,uu_mesh%gauss%n_w) :: u0_cn
    REAL(KIND=8), DIMENSION(2) :: fr  , fth , fz
    REAL(KIND=8), DIMENSION(3,2) :: f   
    REAL(KIND=8)   :: ray
    INTEGER, DIMENSION(uu_mesh%gauss%n_w) :: j_loc   
    REAL(KIND=8), DIMENSION(2)   :: smb   
    REAL(KIND=8):: user_time
    !EXTERNAL :: user_time
    REAL(KIND=8) :: tps, dummy, tt,tps1     


    CALL gauss(uu_mesh)
    jj => uu_mesh%jj
    jj_c => pp_mesh%jj
    me => uu_mesh%me
    nwc => pp_mesh%gauss%n_w

    DO l = 1, l_G   
       w_c(1,l) = ww(1,l) + 0.5*(ww(n_w-1,l) + ww(n_w,l)) 
       w_c(2,l) = ww(2,l) + 0.5*(ww(n_w,l) + ww(4,l)) 
       w_c(3,l) = ww(3,l) + 0.5*(ww(4,l) + ww(n_w-1,l)) 
    END DO

    u0_c = 0.d0   
    tps = 0.d0

    DO m = 1, me
       j_loc(:) = jj(:,m)
       !u0_cn = 0.d0
       DO l = 1, l_G

          !--------On calcule le rayon du point gauss des P2
          ray = 0
          DO n = 1, n_w;  i = jj(n,m)
             ray = ray + uu_mesh%rr(1,i)*ww(n,l)
          END DO
          !----------------------
          tt = user_time(dummy)

          !calcul de la divergence sur les cos
          f(1,1)  = (ray*SUM(gg(1,j_loc(:))*dw(1,:,l,m)) &
               + SUM(gg(1,j_loc(:))*ww(:,l))) 
          f(2,1) = mode*SUM(gg(4,j_loc(:))*ww(:,l)) 
          f(3,1)  = SUM(gg(5,j_loc(:))*dw(2,:,l,m)) * ray

          !calcul de la divergence sur les sin
          f(1,2)  = (ray*SUM(gg(2,j_loc(:))*dw(1,:,l,m)) &
               + SUM(gg(2,j_loc(:))*ww(:,l))) 
          f(2,2) = - mode*SUM(gg(3,j_loc(:))*ww(:,l))
          f(3,2)  = SUM(gg(6,j_loc(:))*dw(2,:,l,m)) * ray

          f = f *rj(l,m)

          tps = tps +user_time(dummy) -tt

          DO j=1,2
             DO n=1, nwc
                u0_c(j,jj_c(n,m)) = u0_c(j,jj_c(n,m)) + w_c(n,l)*(f(1,j)+f(2,j)+f(3,j))
             ENDDO
          ENDDO

       ENDDO

    ENDDO


  END SUBROUTINE qs_01_div_hybrid

  SUBROUTINE qs_01_div_hybrid_2006 (uu_mesh, pp_mesh,mode, gg,  u0_c)
    !================================================

    !  < w_c, D.g >   ===>   u0_c

    USE Gauss_points

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: gg
    REAL(KIND=8), DIMENSION(:,:),   INTENT(OUT) :: u0_c
    INTEGER     ,                   INTENT(IN)  :: mode
    TYPE(mesh_type), TARGET                     :: uu_mesh, pp_mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj, jj_c
    INTEGER,      POINTER       :: me, nwc
    REAL(KIND=8) :: x
    REAL(KIND=8), DIMENSION(pp_mesh%gauss%n_w, uu_mesh%gauss%l_G) ::  w_c
    INTEGER :: m, l, n, k,i,j
    REAL(KIND=8), DIMENSION(3,2) :: f   
    REAL(KIND=8)   :: ray
    INTEGER, DIMENSION(uu_mesh%gauss%n_w) :: j_loc   


    CALL gauss(uu_mesh)
    jj => uu_mesh%jj
    jj_c => pp_mesh%jj
    me => uu_mesh%me
    nwc => pp_mesh%gauss%n_w



    DO l = 1, l_G 
       !cas P1 P2
       w_c(1,l) = ww(1,l) + 0.5*(ww(n_w-1,l) + ww(n_w,l)) 
       w_c(2,l) = ww(2,l) + 0.5*(ww(n_w,l) + ww(4,l)) 
       w_c(3,l) = ww(3,l) + 0.5*(ww(4,l) + ww(n_w-1,l)) 
       !cas P1 P1
       !w_c(1,l) = ww(1,l) 
       !w_c(2,l) = ww(2,l) 
       !w_c(3,l) = ww(3,l) 
    END DO

    u0_c = 0.d0   

    DO m = 1, me
       j_loc(:) = jj(:,m)
       DO l = 1, l_G

          !--------On calcule le rayon du point gauss des P2
          ray = 0
          DO n = 1, n_w;  i = jj(n,m)
             ray = ray + uu_mesh%rr(1,i)*ww(n,l)
          END DO
          !----------------------

          !calcul de la divergence sur les cos
          f(1,1) = (ray*SUM(gg(j_loc,1)*dw(1,:,l,m)) &
               + SUM(gg(j_loc,1)*ww(:,l))) 
          f(2,1) = mode*SUM(gg(j_loc,4)*ww(:,l)) 
          f(3,1) =      SUM(gg(j_loc,5)*dw(2,:,l,m)) * ray

          !calcul de la divergence sur les sin
          f(1,2)  = (ray*SUM(gg(j_loc,2)*dw(1,:,l,m)) &
               + SUM(gg(j_loc,2)*ww(:,l))) 
          f(2,2)  =-mode*SUM(gg(j_loc,3)*ww(:,l))
          f(3,2)  =      SUM(gg(j_loc,6)*dw(2,:,l,m)) * ray

          f = f *rj(l,m)

          DO j=1,2
             x = f(1,j)+f(2,j)+f(3,j)
             DO n=1, nwc
                u0_c(jj_c(n,m),j) = u0_c(jj_c(n,m),j) + w_c(n,l)*x
             ENDDO
          ENDDO

       ENDDO

    ENDDO

  END SUBROUTINE qs_01_div_hybrid_2006

  !-------------------------------------
  SUBROUTINE qs_01_grad(mesh,mode,pp,u0)
    !=====================================
    !calcul le second membre pour le calcul du gradient pour un mode
    USE Gauss_points

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)   :: pp
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT)  :: u0
    INTEGER                     , INTENT(IN)   :: mode   
    INTEGER ::  m, l , i , ni , j
    REAL(KIND=8),DIMENSION(6)                  :: fp


    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8)   :: ray

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    fp = 0.d0
    u0 = 0.d0

    DO m = 1, me
       DO l = 1, l_G
          !--------On calcule le rayon du point gauss
          ray = 0
          DO ni = 1, n_w;  i = jj(ni,m)
             ray = ray + mesh%rr(1,i)*ww(ni,l)
          END DO


          fp(1) = SUM(pp(1,jj(:,m))*dw(1,:,l,m)) * rj(l,m)
          fp(2) = SUM(pp(2,jj(:,m))*dw(1,:,l,m)) * rj(l,m)
          !JLG CN, 27 May 2009
          !Il y avait une inversion pp(1,...) et pp(2,...)
          fp(3) = mode/ray * SUM(pp(2,jj(:,m))*ww(:,l)) * rj(l,m)
          fp(4) = -mode/ray * SUM(pp(1,jj(:,m))*ww(:,l)) * rj(l,m)
          !JLG CN, 27 May 2009 
          fp(5) = SUM(pp(1,jj(:,m))*dw(2,:,l,m)) * rj(l,m)
          fp(6) = SUM(pp(2,jj(:,m))*dw(2,:,l,m)) * rj(l,m)

          DO ni = 1, n_w 
             DO j=1, 6
                u0(j,jj(ni,m)) = u0(j,jj(ni,m)) +  ww(ni,l) * ray * fp(j)
             ENDDO
          ENDDO

       ENDDO
    ENDDO

  END SUBROUTINE qs_01_grad

  SUBROUTINE qs_01_grad_gl(mesh,mod_max,pp,u0)
    !=====================================
    !calcul le second membre pour le calcul du gradient pour m=0 -> mod_max
    USE Gauss_points

    IMPLICIT NONE

    TYPE(mesh_type), TARGET                              :: mesh 
    INTEGER                     , INTENT(IN)             :: mod_max    
    REAL(KIND=8), DIMENSION(2,mesh%np,0:mod_max), INTENT(IN)   :: pp
    REAL(KIND=8), DIMENSION(6,mesh%np,0:mod_max), INTENT(OUT)  :: u0 
    INTEGER ::  m, l , i , ni , j, k
    REAL(KIND=8),DIMENSION(6)                            :: fp



    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8)   :: ray

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    fp = 0.d0
    u0 = 0.d0

    DO k=0, mod_max

       DO m = 1, me
          DO l = 1, l_G
             !--------On calcule le rayon du point gauss
             ray = 0
             DO ni = 1, n_w;  i = jj(ni,m)
                ray = ray + mesh%rr(1,i)*ww(ni,l)
             END DO

             !DO k=0, mod_max

             !je multiplie par le rayon ici, car ray peut etre tres petit
             fp(1) = SUM(pp(1,jj(:,m),k)*dw(1,:,l,m)) * rj(l,m)*ray
             fp(2) = SUM(pp(2,jj(:,m),k)*dw(1,:,l,m)) * rj(l,m)*ray

             fp(3) = k* SUM(pp(2,jj(:,m),k)*ww(:,l)) * rj(l,m)  !/ray * ray
             fp(4) = -k* SUM(pp(1,jj(:,m),k)*ww(:,l)) * rj(l,m) !/ray * ray

             fp(5) = SUM(pp(1,jj(:,m),k)*dw(2,:,l,m)) * rj(l,m)*ray
             fp(6) = SUM(pp(2,jj(:,m),k)*dw(2,:,l,m)) * rj(l,m)*ray

             DO ni = 1, n_w 
                DO j=1, 6
                   u0(j,jj(ni,m),k) = u0(j,jj(ni,m),k) +  ww(ni,l) * fp(j)
                ENDDO
             ENDDO

          ENDDO

       ENDDO

    ENDDO

  END SUBROUTINE qs_01_grad_gl
  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  !procedure pour le cas vectoriel
  !---------------------------------------------------------------------
  !---------------------------------------------------------------------


  SUBROUTINE qs_00_inst_vect_3d(mesh,alpha,mode,ff,V1,V2,dt,u0) !avec  le rayon
    !=================================
    !calcul le second membre en entier pour la diffusion d'un vecteur
    !avec un BDF2
    !  < w, f > +ou-alpha* <w , 2m/r^2*(2V(n)-V(n-1))> + (<w,(4V(n)-V(n-1))/2dt>
    !           ===>   u0 
    !cela depend des composantes et il y a en plus couplage entre les composantes
    !1 et 4, 2 et 3, 5 et 6


    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: ff   !analytique (type,noeud)
    REAL(KIND=8), DIMENSION(:,:),   INTENT(OUT) :: u0   !u0(type,noeud)
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: V1   !V(type,noeud) ,Vn-1
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: V2   !V(type,noeud) ,Vn
    REAL(KIND=8),                   INTENT(IN)  :: dt , alpha
    INTEGER                     ,   INTENT(IN)  :: mode
    INTEGER ::  m, l , i , ni , j
    REAL(KIND=8), DIMENSION(6) :: fs , fexp , ft

    !fs pour le terme de forcage
    !fexp pour le terme de vitesse explicite
    !ft pour le terme temporelle   


    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8)   :: ray

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    u0 = 0
    fs = 0.d0
    fexp = 0.d0
    ft = 0.d0

    DO m = 1, me
       DO l = 1, l_G

          !--------On calcul le rayon du point gauss
          ray = 0
          DO ni = 1, n_w;  i = jj(ni,m)
             ray = ray + mesh%rr(1,i)*ww(ni,l)
          END DO

          !--------calcul de la premiere composante : u0(1,:) <--> f(r,m,c)

          fs(1) = SUM(ff(1,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(1) = -(alpha*2*mode/ray**2)*SUM((2*V2(4,jj(:,m))-V1(4,jj(:,m)))* ww(:,l)) * rj(l,m) 
          ft(1) = 1/(2*dt) * SUM((4*V2(1,jj(:,m))-V1(1,jj(:,m))) * ww(:,l)) * rj(l,m)

          !--------calcul de la seconde composante : u0(2,:) <--> f(r,m,s)

          fs(2) = SUM(ff(2,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(2) = (alpha*2*mode/ray**2)*SUM((2*V2(3,jj(:,m))-V1(3,jj(:,m)))* ww(:,l)) * rj(l,m)  
          ft(2) = 1/(2*dt) * SUM((4*V2(2,jj(:,m))-V1(2,jj(:,m))) * ww(:,l)) * rj(l,m)

          !--------calcul de la troisieme composante : u0(3,:) <--> f(th,m,c)

          fs(3) = SUM(ff(3,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(3) = (alpha*2*mode/ray**2)*SUM((2*V2(2,jj(:,m))-V1(2,jj(:,m)))* ww(:,l)) * rj(l,m)  
          ft(3) = 1/(2*dt) * SUM((4*V2(3,jj(:,m))-V1(3,jj(:,m))) * ww(:,l)) * rj(l,m)

          !--------calcul de la quatrieme composante : u0(4,:) <--> f(th,m,s)

          fs(4) = SUM(ff(4,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(4) = -(alpha*2*mode/ray**2)*SUM((2*V2(1,jj(:,m))-V1(1,jj(:,m)))* ww(:,l)) * rj(l,m) 
          ft(4) = 1/(2*dt) * SUM((4*V2(4,jj(:,m))-V1(4,jj(:,m))) * ww(:,l)) * rj(l,m)

          !--------calcul de la cinquieme composante : u0(5,:) <--> f(z,m,c)

          fs(5) = SUM(ff(5,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(5) = 0.d0
          ft(5) = 1/(2*dt) * SUM((4*V2(5,jj(:,m))-V1(5,jj(:,m))) * ww(:,l)) * rj(l,m)

          !--------calcul de la sixieme composante : u0(6,:) <--> f(z,m,s)

          fs(6) = SUM(ff(6,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(6) = 0.d0
          ft(6) = 1/(2*dt) * SUM((4*V2(6,jj(:,m))-V1(6,jj(:,m))) * ww(:,l)) * rj(l,m)


          !---------- calcul du second membre total

          DO j=1,6
             DO ni = 1, n_w 
                u0(j,jj(ni,m)) = u0(j,jj(ni,m)) +  ww(ni,l)  * ray * (fs(j) + fexp(j) + ft(j)) 
             ENDDO
          ENDDO


       ENDDO
    ENDDO

  END SUBROUTINE qs_00_inst_vect_3d



  SUBROUTINE qs_00_inst_vect_3d_ssr(mesh,alpha,mode,ff,V1,V2,dt,u0) !sans  le rayon
    !pour le calculde l'integral du au forcage
    !=================================
    !calcul le second membre en entier pour la diffusion d'un vecteur
    !avec un BDF2
    !  < w, f > +ou-alpha* <w , 2m/r^2*(2V(n)-V(n-1))> + (<w,(4V(n)-V(n-1))/2dt>
    !           ===>   u0 
    !cela depend des composantes et il y a en plus couplage


    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: ff   !analytique (type,noeud)
    REAL(KIND=8), DIMENSION(:,:),   INTENT(OUT) :: u0   !u0(type,noeud)
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: V1   !T(type,noeud) ,Tn-1
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: V2   !T(type,noeud) ,Tn
    REAL(KIND=8),                   INTENT(IN)  :: dt , alpha
    INTEGER                     ,   INTENT(IN)  :: mode
    INTEGER ::  m, l , i , ni , j
    REAL(KIND=8), DIMENSION(6) :: fs , fexp , ft

    !fs pour le terme de forcage
    !fexp pour le terme de vitesse explicite
    !ft pour le terme temporelle   


    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8)   :: ray

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    u0 = 0

    DO m = 1, me
       DO l = 1, l_G

          !--------On calcul le rayon du point gauss
          ray = 0
          DO ni = 1, n_w;  i = jj(ni,m)
             ray = ray + mesh%rr(1,i)*ww(ni,l)
          END DO

          !--------calcul de la premiere composante : u0(1,:) <--> f(r,m,c)

          fs(1) = SUM(ff(1,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(1) = -(alpha*2*mode/ray**2)*SUM((2*V2(4,jj(:,m))-V1(4,jj(:,m)))* ww(:,l)) * rj(l,m) 
          ft(1) = 1/(2*dt) * SUM((4*V2(1,jj(:,m))-V1(1,jj(:,m))) * ww(:,l)) * rj(l,m)

          !--------calcul de la seconde composante : u0(2,:) <--> f(r,m,s)

          fs(2) = SUM(ff(2,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(2) = (alpha*2*mode/ray**2)*SUM((2*V2(3,jj(:,m))-V1(3,jj(:,m)))* ww(:,l)) * rj(l,m)  
          ft(2) = 1/(2*dt) * SUM((4*V2(2,jj(:,m))-V1(2,jj(:,m))) * ww(:,l)) * rj(l,m)

          !--------calcul de la troisieme composante : u0(3,:) <--> f(th,m,c)

          fs(3) = SUM(ff(3,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(3) = (alpha*2*mode/ray**2)*SUM((2*V2(2,jj(:,m))-V1(2,jj(:,m)))* ww(:,l)) * rj(l,m)  
          ft(3) = 1/(2*dt) * SUM((4*V2(3,jj(:,m))-V1(3,jj(:,m))) * ww(:,l)) * rj(l,m)

          !--------calcul de la quatrieme composante : u0(4,:) <--> f(th,m,s)

          fs(4) = SUM(ff(4,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(4) = -(alpha*2*mode/ray**2)*SUM((2*V2(1,jj(:,m))-V1(1,jj(:,m)))* ww(:,l)) * rj(l,m) 
          ft(4) = 1/(2*dt) * SUM((4*V2(4,jj(:,m))-V1(4,jj(:,m))) * ww(:,l)) * rj(l,m)

          !--------calcul de la cinquieme composante : u0(5,:) <--> f(z,m,c)

          fs(5) = SUM(ff(5,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(5) = 0.d0
          ft(5) = 1/(2*dt) * SUM((4*V2(5,jj(:,m))-V1(5,jj(:,m))) * ww(:,l)) * rj(l,m)

          !--------calcul de la sixieme composante : u0(6,:) <--> f(z,m,s)

          fs(6) = SUM(ff(6,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(6) = 0.d0
          ft(6) = 1/(2*dt) * SUM((4*V2(6,jj(:,m))-V1(6,jj(:,m))) * ww(:,l)) * rj(l,m)


          !---------- calcul du second membre total


          DO j=1,6
             DO ni = 1, n_w 
                u0(j,jj(ni,m)) = u0(j,jj(ni,m)) +  ww(ni,l)  *  (fs(j) + ray*(fexp(j) + ft(j))) 
             ENDDO
          ENDDO


       ENDDO
    ENDDO

  END SUBROUTINE qs_00_inst_vect_3d_ssr


  SUBROUTINE qs_00_inst_vect_3d_init(mesh,alpha,mode,ff,V2,dt,u0) !avec  le rayon
    !=================================
    !calcul le second membre en entier pour la diffusion d'un vecteur
    !avec un Euler retarde
    !calcul le second membre en entier avce les termes explicites
    !  < w, f > +ou-alpha* <w , 2m/r^2*T(n)> + (<w,V(n)/dt>
    !           ===>   u0 
    !cela depend des composantes et il y a en plus couplage


    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: ff   !analytique (type,noeud)
    REAL(KIND=8), DIMENSION(:,:),   INTENT(OUT) :: u0   !u0(type,noeud)
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: V2   !T(type,noeud) ,Tn
    REAL(KIND=8),                   INTENT(IN)  :: dt , alpha
    INTEGER                     ,   INTENT(IN)  :: mode
    INTEGER ::  m, l , i , ni , j
    REAL(KIND=8), DIMENSION(6) :: fs , fexp , ft

    !fs pour le terme de forcage
    !fexp pour le terme de vitesse explicite
    !ft pour le terme temporelle   


    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8)   :: ray

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    u0 = 0

    DO m = 1, me
       DO l = 1, l_G

          !--------On calcul le rayon du point gauss
          ray = 0
          DO ni = 1, n_w;  i = jj(ni,m)
             ray = ray + mesh%rr(1,i)*ww(ni,l)
          END DO

          !--------calcul de la premiere composante : u0(1,:) <--> f(r,m,c)

          fs(1) = SUM(ff(1,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(1) = (-alpha*2*mode/ray**2)*SUM(V2(4,jj(:,m))* ww(:,l)) * rj(l,m)  
          ft(1) = 1/(dt) * SUM(V2(1,jj(:,m)) * ww(:,l)) * rj(l,m)

          !--------calcul de la seconde composante : u0(2,:) <--> f(r,m,s)

          fs(2) = SUM(ff(2,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(2) = (alpha*2*mode/ray**2)*SUM(V2(3,jj(:,m))* ww(:,l)) * rj(l,m)  
          ft(2) = 1/(dt) * SUM(V2(2,jj(:,m)) * ww(:,l)) * rj(l,m)

          !--------calcul de la troisieme composante : u0(3,:) <--> f(th,m,c)

          fs(3) = SUM(ff(3,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(3) = (alpha*2*mode/ray**2)*SUM(V2(2,jj(:,m))* ww(:,l)) * rj(l,m)  
          ft(3) = 1/(dt) * SUM(V2(3,jj(:,m)) * ww(:,l)) * rj(l,m)

          !--------calcul de la quatrieme composante : u0(4,:) <--> f(th,m,s)

          fs(4) = SUM(ff(4,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(4) = (-alpha*2*mode/ray**2)*SUM(V2(1,jj(:,m))* ww(:,l)) * rj(l,m)  
          ft(4) = 1/(dt) * SUM(V2(4,jj(:,m)) * ww(:,l)) * rj(l,m)

          !--------calcul de la cinquieme composante : u0(5,:) <--> f(z,m,c)

          fs(5) = SUM(ff(5,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(5) = 0.d0
          ft(5) = 1/(dt) * SUM(V2(5,jj(:,m)) * ww(:,l)) * rj(l,m)

          !--------calcul de la sixieme composante : u0(6,:) <--> f(z,m,s)

          fs(6) = SUM(ff(6,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(6) = 0.d0
          ft(6) = 1/(dt) * SUM(V2(6,jj(:,m)) * ww(:,l)) * rj(l,m)

          !---------- calcul du second membre total

          DO j=1,6
             DO ni = 1, n_w 
                u0(j,jj(ni,m)) = u0(j,jj(ni,m)) +  ww(ni,l)  * ray * (fs(j) + fexp(j) + ft(j)) 
             ENDDO
          ENDDO


       ENDDO
    ENDDO

  END SUBROUTINE qs_00_inst_vect_3d_init


  SUBROUTINE qs_00_inst_vect_3d_init_ssr(mesh,alpha,mode,ff,V,dt,u0) !sans  le rayon
    !dans l'intgration du forcage
    !=================================
    !calcul le second membre en entier pour la diffusion d'un vecteur
    !avec un Euler retarde
    !calcul le second membre en entier avce les termes explicites
    !  < w, f > +ou-alpha* <w , 2m/r^2*T(n)> + (<w,V(n)/dt>
    !           ===>   u0 
    !cela depend des composantes et il y a en plus couplage


    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: ff   !analytique (type,noeud)
    REAL(KIND=8), DIMENSION(:,:),   INTENT(OUT) :: u0   !u0(type,noeud)
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: V   !T(type,noeud) ,Tn
    REAL(KIND=8),                   INTENT(IN)  :: dt , alpha
    INTEGER                     ,   INTENT(IN)  :: mode
    INTEGER ::  m, l , i , ni , j
    REAL(KIND=8), DIMENSION(6) :: fs , fexp , ft

    !fs pour le terme de forcage
    !fexp pour le terme de vitesse explicite
    !ft pour le terme temporelle   


    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8)   :: ray

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    u0 = 0
    fs = 0.d0
    ft = 0.d0
    fexp = 0.d0

    DO m = 1, me
       DO l = 1, l_G
          !--------On calcule le rayon du point gauss
          ray = 0
          DO ni = 1, n_w;  i = jj(ni,m)
             ray = ray + mesh%rr(1,i)*ww(ni,l)
          END DO

          !--------calcul de la premiere composante : u0(1,:) <--> f(r,m,c)

          fs(1) = SUM(ff(1,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(1) = -(alpha*2*mode/ray**2)*SUM(V(4,jj(:,m))* ww(:,l)) * rj(l,m)  
          ft(1) = 1/(dt) * SUM(V(1,jj(:,m)) * ww(:,l)) * rj(l,m)

          !--------calcul de la seconde composante : u0(2,:) <--> f(r,m,s)

          fs(2) = SUM(ff(2,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(2) = (alpha*2*mode/ray**2)*SUM(V(3,jj(:,m))* ww(:,l)) * rj(l,m)  
          ft(2) = 1/(dt) * SUM(V(2,jj(:,m)) * ww(:,l)) * rj(l,m)

          !--------calcul de la troisieme composante : u0(3,:) <--> f(th,m,c)

          fs(3) = SUM(ff(3,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(3) = (alpha*2*mode/ray**2)*SUM(V(2,jj(:,m))* ww(:,l)) * rj(l,m)  
          ft(3) = 1/(dt) * SUM(V(3,jj(:,m)) * ww(:,l)) * rj(l,m)

          !--------calcul de la quatrieme composante : u0(4,:) <--> f(th,m,s)

          fs(4) = SUM(ff(4,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(4) = -(alpha*2*mode/ray**2)*SUM(V(1,jj(:,m))* ww(:,l)) * rj(l,m)  
          ft(4) = 1/(dt) * SUM(V(4,jj(:,m)) * ww(:,l)) * rj(l,m)

          !--------calcul de la cinquieme composante : u0(5,:) <--> f(z,m,c)

          fs(5) = SUM(ff(5,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(5) = 0.d0
          ft(5) = 1/(dt) * SUM(V(5,jj(:,m)) * ww(:,l)) * rj(l,m)

          !--------calcul de la sixieme composante : u0(6,:) <--> f(z,m,s)

          fs(6) = SUM(ff(6,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(6) = 0.d0
          ft(6) = 1/(dt) * SUM(V(6,jj(:,m)) * ww(:,l)) * rj(l,m)        

          !---------- calcul du second membre total



          DO j=1,6
             DO ni = 1, n_w 
                u0(j,jj(ni,m)) = u0(j,jj(ni,m)) +  ww(ni,l)  *  (fs(j) + ray*(fexp(j) + ft(j))) 
             ENDDO
          ENDDO


       ENDDO
    ENDDO

  END SUBROUTINE qs_00_inst_vect_3d_init_ssr


  SUBROUTINE qs_00_adv_diff_vect_3d_init_ssr(mesh,alpha,mode,gg,ff,V,dt,u0) !sans  le rayon
    !dans l'intgration du forcage
    !=================================



    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: ff   !analytique (type,noeud)
    REAL(KIND=8), DIMENSION(:,:),   INTENT(OUT) :: u0   !u0(type,noeud)
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: V   !T(type,noeud) ,Tn
    REAL(KIND=8),                   INTENT(IN)  :: dt , alpha
    INTEGER                     ,   INTENT(IN)  :: mode
    INTEGER ::  m, l , i , ni , j
    REAL(KIND=8), DIMENSION(6) :: fs , fexp , ft , fv
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg   !champ de vitesse V_theta
    !fs pour le terme de forcage
    !fexp pour le terme de vitesse explicite
    !ft pour le terme temporelle   


    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8)   :: ray

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    u0 = 0
    fs = 0.d0
    ft = 0.d0
    fexp = 0.d0
    fv = 0.d0   
    DO m = 1, me
       DO l = 1, l_G
          !--------On calcule le rayon du point gauss
          ray = 0
          DO ni = 1, n_w;  i = jj(ni,m)
             ray = ray + mesh%rr(1,i)*ww(ni,l)
          END DO

          !--------calcul de la premiere composante : u0(1,:) <--> f(r,m,c)

          fs(1) = SUM(ff(1,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(1) = -(alpha*2*mode/ray**2)*SUM(V(4,jj(:,m))* ww(:,l)) * rj(l,m)  
          ft(1) = 1/(dt) * SUM(V(1,jj(:,m)) * ww(:,l)) * rj(l,m)
          fv(1) = -mode/ray*SUM(gg(3,jj(:,m))*V(2,jj(:,m))* ww(:,l))* rj(l,m)

          !--------calcul de la seconde composante : u0(2,:) <--> f(r,m,s)

          fs(2) = SUM(ff(2,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(2) = (alpha*2*mode/ray**2)*SUM(V(3,jj(:,m))* ww(:,l)) * rj(l,m)  
          ft(2) = 1/(dt) * SUM(V(2,jj(:,m)) * ww(:,l)) * rj(l,m)
          fv(2) = mode/ray*SUM(gg(3,jj(:,m))*V(1,jj(:,m))* ww(:,l))* rj(l,m)

          !--------calcul de la troisieme composante : u0(3,:) <--> f(th,m,c)

          fs(3) = SUM(ff(3,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(3) = (alpha*2*mode/ray**2)*SUM(V(2,jj(:,m))* ww(:,l)) * rj(l,m)  
          ft(3) = 1/(dt) * SUM(V(3,jj(:,m)) * ww(:,l)) * rj(l,m)
          fv(3) = -mode/ray*SUM(gg(3,jj(:,m))*V(4,jj(:,m))* ww(:,l))* rj(l,m)

          !--------calcul de la quatrieme composante : u0(4,:) <--> f(th,m,s)

          fs(4) = SUM(ff(4,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(4) = -(alpha*2*mode/ray**2)*SUM(V(1,jj(:,m))* ww(:,l)) * rj(l,m)  
          ft(4) = 1/(dt) * SUM(V(4,jj(:,m)) * ww(:,l)) * rj(l,m)
          fv(4) = mode/ray*SUM(gg(3,jj(:,m))*V(3,jj(:,m))* ww(:,l))* rj(l,m)

          !--------calcul de la cinquieme composante : u0(5,:) <--> f(z,m,c)

          fs(5) = SUM(ff(5,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(5) = 0.d0
          ft(5) = 1/(dt) * SUM(V(5,jj(:,m)) * ww(:,l)) * rj(l,m)
          fv(5) = -mode/ray*SUM(gg(3,jj(:,m))*V(6,jj(:,m))* ww(:,l))* rj(l,m)

          !--------calcul de la sixieme composante : u0(6,:) <--> f(z,m,s)

          fs(6) = SUM(ff(6,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(6) = 0.d0
          ft(6) = 1/(dt) * SUM(V(6,jj(:,m)) * ww(:,l)) * rj(l,m)        
          fv(6) = mode/ray*SUM(gg(3,jj(:,m))*V(5,jj(:,m))* ww(:,l))* rj(l,m)

          !---------- calcul du second membre total



          DO j=1,6
             DO ni = 1, n_w 
                u0(j,jj(ni,m)) = u0(j,jj(ni,m)) +  ww(ni,l)  *  (fs(j) + ray*(fexp(j) + ft(j) + fv(j))) 
             ENDDO
          ENDDO


       ENDDO
    ENDDO

  END SUBROUTINE qs_00_adv_diff_vect_3d_init_ssr


  SUBROUTINE qs_00_adv_diff_vect_3d_init(mesh,alpha,mode,gg,ff,V,dt,u0) !avec  le rayon
    !dans l'intgration du forcage
    !=================================



    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: ff   !analytique (type,noeud)
    REAL(KIND=8), DIMENSION(:,:),   INTENT(OUT) :: u0   !u0(type,noeud)
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: V   !T(type,noeud) ,Tn
    REAL(KIND=8),                   INTENT(IN)  :: dt , alpha
    INTEGER                     ,   INTENT(IN)  :: mode
    INTEGER ::  m, l , i , ni , j
    REAL(KIND=8), DIMENSION(6) :: fs , fexp , ft , fv
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg   !champ de vitesse V_theta
    !fs pour le terme de forcage
    !fexp pour le terme de vitesse explicite
    !ft pour le terme temporelle   


    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8)   :: ray

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    u0 = 0
    fs = 0.d0
    ft = 0.d0
    fexp = 0.d0
    fv = 0.d0   
    DO m = 1, me
       DO l = 1, l_G
          !--------On calcule le rayon du point gauss
          ray = 0
          DO ni = 1, n_w;  i = jj(ni,m)
             ray = ray + mesh%rr(1,i)*ww(ni,l)
          END DO

          !--------calcul de la premiere composante : u0(1,:) <--> f(r,m,c)

          fs(1) = SUM(ff(1,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(1) = -(alpha*2*mode/ray**2)*SUM(V(4,jj(:,m))* ww(:,l)) * rj(l,m)  
          ft(1) = 1/(dt) * SUM(V(1,jj(:,m)) * ww(:,l)) * rj(l,m)
          fv(1) = -mode/ray*SUM(gg(3,jj(:,m))*V(2,jj(:,m))* ww(:,l))* rj(l,m)

          !--------calcul de la seconde composante : u0(2,:) <--> f(r,m,s)

          fs(2) = SUM(ff(2,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(2) = (alpha*2*mode/ray**2)*SUM(V(3,jj(:,m))* ww(:,l)) * rj(l,m)  
          ft(2) = 1/(dt) * SUM(V(2,jj(:,m)) * ww(:,l)) * rj(l,m)
          fv(2) = mode/ray*SUM(gg(3,jj(:,m))*V(1,jj(:,m))* ww(:,l))* rj(l,m)

          !--------calcul de la troisieme composante : u0(3,:) <--> f(th,m,c)

          fs(3) = SUM(ff(3,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(3) = (alpha*2*mode/ray**2)*SUM(V(2,jj(:,m))* ww(:,l)) * rj(l,m)  
          ft(3) = 1/(dt) * SUM(V(3,jj(:,m)) * ww(:,l)) * rj(l,m)
          fv(3) = -mode/ray*SUM(gg(3,jj(:,m))*V(4,jj(:,m))* ww(:,l))* rj(l,m)

          !--------calcul de la quatrieme composante : u0(4,:) <--> f(th,m,s)

          fs(4) = SUM(ff(4,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(4) = -(alpha*2*mode/ray**2)*SUM(V(1,jj(:,m))* ww(:,l)) * rj(l,m)  
          ft(4) = 1/(dt) * SUM(V(4,jj(:,m)) * ww(:,l)) * rj(l,m)
          fv(4) = mode/ray*SUM(gg(3,jj(:,m))*V(3,jj(:,m))* ww(:,l))* rj(l,m)

          !--------calcul de la cinquieme composante : u0(5,:) <--> f(z,m,c)

          fs(5) = SUM(ff(5,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(5) = 0.d0
          ft(5) = 1/(dt) * SUM(V(5,jj(:,m)) * ww(:,l)) * rj(l,m)
          fv(5) = -mode/ray*SUM(gg(3,jj(:,m))*V(6,jj(:,m))* ww(:,l))* rj(l,m)

          !--------calcul de la sixieme composante : u0(6,:) <--> f(z,m,s)

          fs(6) = SUM(ff(6,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(6) = 0.d0
          ft(6) = 1/(dt) * SUM(V(6,jj(:,m)) * ww(:,l)) * rj(l,m)        
          fv(6) = mode/ray*SUM(gg(3,jj(:,m))*V(5,jj(:,m))* ww(:,l))* rj(l,m)

          !---------- calcul du second membre total



          DO j=1,6
             DO ni = 1, n_w 
                u0(j,jj(ni,m)) = u0(j,jj(ni,m)) +  ww(ni,l)  *  ray*(fs(j) + fexp(j) + ft(j) + fv(j)) 
             ENDDO
          ENDDO


       ENDDO
    ENDDO

  END SUBROUTINE qs_00_adv_diff_vect_3d_init

  SUBROUTINE qs_00_adv_diff_vect_3d_ssr(mesh,alpha,mode,gg,ff,V1,V2,dt,u0) !sans  le rayon
    !pour le calculde l'integral du au forcage
    !=================================


    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: ff   !analytique (type,noeud)
    REAL(KIND=8), DIMENSION(:,:),   INTENT(OUT) :: u0   !u0(type,noeud)
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: V1   !T(type,noeud) ,Tn-1
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: V2   !T(type,noeud) ,Tn
    REAL(KIND=8),                   INTENT(IN)  :: dt , alpha
    INTEGER                     ,   INTENT(IN)  :: mode
    INTEGER ::  m, l , i , ni , j
    REAL(KIND=8), DIMENSION(6) :: fs , fexp , ft , fv

    !fs pour le terme de forcage
    !fexp pour le terme de vitesse explicite
    !ft pour le terme temporelle   


    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8)   :: ray
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg   !champ de vitesse V_theta
    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    u0 = 0

    fs = 0.d0
    fexp = 0.d0
    ft = 0.d0
    fv = 0.d0

    DO m = 1, me
       DO l = 1, l_G

          !--------On calcul le rayon du point gauss
          ray = 0
          DO ni = 1, n_w;  i = jj(ni,m)
             ray = ray + mesh%rr(1,i)*ww(ni,l)
          END DO

          !--------calcul de la premiere composante : u0(1,:) <--> f(r,m,c)

          fs(1) = SUM(ff(1,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(1) = -(alpha*2*mode/ray**2)*SUM((2*V2(4,jj(:,m))-V1(4,jj(:,m)))* ww(:,l)) * rj(l,m) 
          ft(1) = 1/(2*dt) * SUM((4*V2(1,jj(:,m))-V1(1,jj(:,m))) * ww(:,l)) * rj(l,m)
          fv(1) = -mode/ray*SUM(gg(3,jj(:,m))*(2*V2(2,jj(:,m)) - V1(2,jj(:,m)))* ww(:,l))*rj(l,m)

          !--------calcul de la seconde composante : u0(2,:) <--> f(r,m,s)

          fs(2) = SUM(ff(2,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(2) = (alpha*2*mode/ray**2)*SUM((2*V2(3,jj(:,m))-V1(3,jj(:,m)))* ww(:,l))*rj(l,m)  
          ft(2) = 1/(2*dt) * SUM((4*V2(2,jj(:,m))-V1(2,jj(:,m))) * ww(:,l)) *rj(l,m)
          fv(2) = mode/ray*SUM(gg(3,jj(:,m))*(2*V2(1,jj(:,m)) - V1(1,jj(:,m)))* ww(:,l))*rj(l,m)

          !--------calcul de la troisieme composante : u0(3,:) <--> f(th,m,c)

          fs(3) = SUM(ff(3,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(3) = (alpha*2*mode/ray**2)*SUM((2*V2(2,jj(:,m))-V1(2,jj(:,m)))* ww(:,l)) * rj(l,m) 
          ft(3) = 1/(2*dt) * SUM((4*V2(3,jj(:,m))-V1(3,jj(:,m))) * ww(:,l)) * rj(l,m)
          fv(3) = -mode/ray*SUM(gg(3,jj(:,m))*(2*V2(4,jj(:,m)) - V1(4,jj(:,m)))* ww(:,l))*rj(l,m)

          !--------calcul de la quatrieme composante : u0(4,:) <--> f(th,m,s)

          fs(4) = SUM(ff(4,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(4) = -(alpha*2*mode/ray**2)*SUM((2*V2(1,jj(:,m))-V1(1,jj(:,m)))* ww(:,l)) * rj(l,m)         
          ft(4) = 1/(2*dt) * SUM((4*V2(4,jj(:,m))-V1(4,jj(:,m))) * ww(:,l)) * rj(l,m)
          fv(4) = mode/ray*SUM(gg(3,jj(:,m))*(2*V2(3,jj(:,m)) - V1(3,jj(:,m)))* ww(:,l))*rj(l,m)

          !--------calcul de la cinquieme composante : u0(5,:) <--> f(z,m,c)

          fs(5) = SUM(ff(5,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(5) = 0.d0
          ft(5) = 1/(2*dt) * SUM((4*V2(5,jj(:,m))-V1(5,jj(:,m))) * ww(:,l)) * rj(l,m)
          fv(5) = -mode/ray*SUM(gg(3,jj(:,m))*(2*V2(6,jj(:,m)) - V1(6,jj(:,m)))* ww(:,l))*rj(l,m)

          !--------calcul de la sixieme composante : u0(6,:) <--> f(z,m,s)

          fs(6) = SUM(ff(6,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(6) =0.d0
          ft(6) = 1/(2*dt) * SUM((4*V2(6,jj(:,m))-V1(6,jj(:,m))) * ww(:,l)) * rj(l,m)
          fv(6) = mode/ray*SUM(gg(3,jj(:,m))*(2*V2(5,jj(:,m)) - V1(5,jj(:,m)))* ww(:,l))*rj(l,m)

          !---------- calcul du second membre total


          DO j=1,6
             DO ni = 1, n_w 
                u0(j,jj(ni,m)) = u0(j,jj(ni,m)) +  ww(ni,l)  *  (fs(j) + ray*(fexp(j) + ft(j) + fv(j))) 
             ENDDO
          ENDDO


       ENDDO
    ENDDO

  END SUBROUTINE qs_00_adv_diff_vect_3d_ssr

  SUBROUTINE qs_00_adv_diff_vect_3d(mesh,alpha,mode,gg,ff,V1,V2,dt,u0) !sans  le rayon
    !pour le calculde l'integral du au forcage
    !=================================


    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: ff   !analytique (type,noeud)
    REAL(KIND=8), DIMENSION(:,:),   INTENT(OUT) :: u0   !u0(type,noeud)
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: V1   !T(type,noeud) ,Tn-1
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: V2   !T(type,noeud) ,Tn
    REAL(KIND=8),                   INTENT(IN)  :: dt , alpha
    INTEGER                     ,   INTENT(IN)  :: mode
    INTEGER ::  m, l , i , ni , j
    REAL(KIND=8), DIMENSION(6) :: fs , fexp , ft , fv

    !fs pour le terme de forcage
    !fexp pour le terme de vitesse explicite
    !ft pour le terme temporelle   


    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8)   :: ray
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg   !champ de vitesse V_theta
    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    u0 = 0

    fs = 0.d0
    fexp = 0.d0
    ft = 0.d0
    fv = 0.d0

    DO m = 1, me
       DO l = 1, l_G

          !--------On calcul le rayon du point gauss
          ray = 0
          DO ni = 1, n_w;  i = jj(ni,m)
             ray = ray + mesh%rr(1,i)*ww(ni,l)
          END DO

          !--------calcul de la premiere composante : u0(1,:) <--> f(r,m,c)

          fs(1) = SUM(ff(1,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(1) = -(alpha*2*mode/ray**2)*SUM((2*V2(4,jj(:,m))-V1(4,jj(:,m)))* ww(:,l)) * rj(l,m) 
          ft(1) = 1/(2*dt) * SUM((4*V2(1,jj(:,m))-V1(1,jj(:,m))) * ww(:,l)) * rj(l,m)
          fv(1) = -mode/ray*SUM(gg(3,jj(:,m))*(2*V2(2,jj(:,m)) - V1(2,jj(:,m)))* ww(:,l))*rj(l,m)

          !--------calcul de la seconde composante : u0(2,:) <--> f(r,m,s)

          fs(2) = SUM(ff(2,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(2) = (alpha*2*mode/ray**2)*SUM((2*V2(3,jj(:,m))-V1(3,jj(:,m)))* ww(:,l))*rj(l,m)  
          ft(2) = 1/(2*dt) * SUM((4*V2(2,jj(:,m))-V1(2,jj(:,m))) * ww(:,l)) *rj(l,m)
          fv(2) = mode/ray*SUM(gg(3,jj(:,m))*(2*V2(1,jj(:,m)) - V1(1,jj(:,m)))* ww(:,l))*rj(l,m)

          !--------calcul de la troisieme composante : u0(3,:) <--> f(th,m,c)

          fs(3) = SUM(ff(3,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(3) = (alpha*2*mode/ray**2)*SUM((2*V2(2,jj(:,m))-V1(2,jj(:,m)))* ww(:,l)) * rj(l,m) 
          ft(3) = 1/(2*dt) * SUM((4*V2(3,jj(:,m))-V1(3,jj(:,m))) * ww(:,l)) * rj(l,m)
          fv(3) = -mode/ray*SUM(gg(3,jj(:,m))*(2*V2(4,jj(:,m)) - V1(4,jj(:,m)))* ww(:,l))*rj(l,m)

          !--------calcul de la quatrieme composante : u0(4,:) <--> f(th,m,s)

          fs(4) = SUM(ff(4,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(4) = -(alpha*2*mode/ray**2)*SUM((2*V2(1,jj(:,m))-V1(1,jj(:,m)))* ww(:,l)) * rj(l,m)         
          ft(4) = 1/(2*dt) * SUM((4*V2(4,jj(:,m))-V1(4,jj(:,m))) * ww(:,l)) * rj(l,m)
          fv(4) = mode/ray*SUM(gg(3,jj(:,m))*(2*V2(3,jj(:,m)) - V1(3,jj(:,m)))* ww(:,l))*rj(l,m)

          !--------calcul de la cinquieme composante : u0(5,:) <--> f(z,m,c)

          fs(5) = SUM(ff(5,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(5) = 0.d0
          ft(5) = 1/(2*dt) * SUM((4*V2(5,jj(:,m))-V1(5,jj(:,m))) * ww(:,l)) * rj(l,m)
          fv(5) = -mode/ray*SUM(gg(3,jj(:,m))*(2*V2(6,jj(:,m)) - V1(6,jj(:,m)))* ww(:,l))*rj(l,m)

          !--------calcul de la sixieme composante : u0(6,:) <--> f(z,m,s)

          fs(6) = SUM(ff(6,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(6) =0.d0
          ft(6) = 1/(2*dt) * SUM((4*V2(6,jj(:,m))-V1(6,jj(:,m))) * ww(:,l)) * rj(l,m)
          fv(6) = mode/ray*SUM(gg(3,jj(:,m))*(2*V2(5,jj(:,m)) - V1(5,jj(:,m)))* ww(:,l))*rj(l,m)

          !---------- calcul du second membre total


          DO j=1,6
             DO ni = 1, n_w 
                u0(j,jj(ni,m)) = u0(j,jj(ni,m)) +  ww(ni,l)  * ray* (fs(j) + fexp(j) + ft(j) + fv(j)) 
             ENDDO
          ENDDO


       ENDDO
    ENDDO

  END SUBROUTINE qs_00_adv_diff_vect_3d


  SUBROUTINE qs_00_stokes_3d_init(mesh,alpha,mode,ff,P,V,dt,u0) !sans  le rayon
    !dans l'intgration du forcage
    !=================================
    !calcul le second membre en entier pour le probleme de Stokes avec un Euler retarde
    !  < w, f > +ou-alpha* <w , 2m/r^2*T(n)> + (<w,(V(n))/dt> + <w,grad(P)>
    !           ===>   u0 
    !cela depend des composantes et il y a en plus couplage


    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: ff   !analytique (type,noeud)
    REAL(KIND=8), DIMENSION(:,:),   INTENT(OUT) :: u0   !u0(type,noeud)
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: V   !T(type,noeud) ,Tn
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: P
    REAL(KIND=8),                   INTENT(IN)  :: dt , alpha
    INTEGER                     ,   INTENT(IN)  :: mode
    INTEGER ::  m, l , i , ni , j
    REAL(KIND=8), DIMENSION(6) :: fs , fexp , ft , fp

    !fs pour le terme de forcage
    !fexp pour le terme de vitesse explicite
    !ft pour le terme temporelle   
    !fp gradient de pression

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8)   :: ray

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    u0 = 0
    fs = 0.d0
    ft = 0.d0
    fexp = 0.d0
    fp = 0.d0

    DO m = 1, me
       DO l = 1, l_G
          !--------On calcule le rayon du point gauss
          ray = 0
          DO ni = 1, n_w;  i = jj(ni,m)
             ray = ray + mesh%rr(1,i)*ww(ni,l)
          END DO

          !--------calcul de la premiere composante : u0(1,:) <--> f(r,m,c)

          fs(1) = SUM(ff(1,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(1) = -(alpha*2*mode/ray**2)*SUM(V(4,jj(:,m))* ww(:,l)) * rj(l,m)  
          ft(1) = 1/(dt) * SUM(V(1,jj(:,m)) * ww(:,l)) * rj(l,m)
          !fp(1) = SUM(P(1,jj(:,m))*ww(:,l)) * rj(l,m)
          fp(1) = SUM(P(1,jj(:,m))*dw(1,:,l,m)) * rj(l,m)
          !--------calcul de la seconde composante : u0(2,:) <--> f(r,m,s)

          fs(2) = SUM(ff(2,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(2) = (alpha*2*mode/ray**2)*SUM(V(3,jj(:,m))* ww(:,l)) * rj(l,m)  
          ft(2) = 1/(dt) * SUM(V(2,jj(:,m)) * ww(:,l)) * rj(l,m)
          !fp(2) = SUM(P(2,jj(:,m))*ww(:,l)) * rj(l,m)
          fp(2) = SUM(P(2,jj(:,m))*dw(1,:,l,m)) * rj(l,m)
          !--------calcul de la troisieme composante : u0(3,:) <--> f(th,m,c)

          fs(3) = SUM(ff(3,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(3) = (alpha*2*mode/ray**2)*SUM(V(2,jj(:,m))* ww(:,l)) * rj(l,m)  
          ft(3) = 1/(dt) * SUM(V(3,jj(:,m)) * ww(:,l)) * rj(l,m)
          fp(3) = SUM(P(2,jj(:,m))*ww(:,l)) * rj(l,m)/ray*mode

          !--------calcul de la quatrieme composante : u0(4,:) <--> f(th,m,s)

          fs(4) = SUM(ff(4,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(4) = -(alpha*2*mode/ray**2)*SUM(V(1,jj(:,m))* ww(:,l)) * rj(l,m)  
          ft(4) = 1/(dt) * SUM(V(4,jj(:,m)) * ww(:,l)) * rj(l,m)
          fp(4) = -SUM(P(1,jj(:,m))*ww(:,l)) * rj(l,m)/ray *mode

          !--------calcul de la cinquieme composante : u0(5,:) <--> f(z,m,c)

          fs(5) = SUM(ff(5,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(5) = 0.d0
          ft(5) = 1/(dt) * SUM(V(5,jj(:,m)) * ww(:,l)) * rj(l,m)
          !fp(5) = SUM(P(1,jj(:,m))*ww(:,l)) * rj(l,m)
          fp(5) = SUM(P(1,jj(:,m))*dw(2,:,l,m)) * rj(l,m)
          !--------calcul de la sixieme composante : u0(6,:) <--> f(z,m,s)

          fs(6) = SUM(ff(6,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(6) = 0.d0
          ft(6) = 1/(dt) * SUM(V(6,jj(:,m)) * ww(:,l)) * rj(l,m)        
          !fp(6) = SUM(P(2,jj(:,m))*ww(:,l)) * rj(l,m) 
          fp(6) = SUM(P(2,jj(:,m))*dw(2,:,l,m)) * rj(l,m)
          !---------- calcul du second membre total



          DO j=1,6
             DO ni = 1, n_w 
                u0(j,jj(ni,m)) = u0(j,jj(ni,m)) +  ww(ni,l)  * ray* (fs(j) + fexp(j) + ft(j) + fp(j)) 

                !IF (j == 3.or.j == 4) THEN
                !   u0(j,jj(ni,m)) = u0(j,jj(ni,m)) +  ww(ni,l) * ray *fp(j)
                !ELSEIF (j == 1.or. j == 2) THEN
                !   u0(j,jj(ni,m)) = u0(j,jj(ni,m)) +  ww(ni,l) * ray * fp(j)
                !ELSEIF (j == 5.or. j == 6) THEN
                !   u0(j,jj(ni,m)) = u0(j,jj(ni,m)) +  ww(ni,l) * ray * fp(j)
                !ELSE
                !ELSEIF (j == 1.or. j == 2) THEN
                !   u0(j,jj(ni,m)) = u0(j,jj(ni,m)) +  dw(1,ni,l,m) * ray * fp(j)
                !ELSEIF (j == 5.or. j == 6) THEN
                !   u0(j,jj(ni,m)) = u0(j,jj(ni,m)) +  dw(2,ni,l,m) * ray * fp(j)
                !ELSE
                !   WRITE(*,*) 'probleme ds calcul second membre avec j,subrout stokes'
                !ENDIF
             ENDDO
          ENDDO


       ENDDO
    ENDDO

  END SUBROUTINE qs_00_stokes_3d_init


  SUBROUTINE qs_00_stokes_3d(mesh,alpha,mode,ff,V1,V2,P,dt,u0) !avec  le rayon
    !=================================



    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: ff   !analytique (type,noeud)
    REAL(KIND=8), DIMENSION(:,:),   INTENT(OUT) :: u0   !u0(type,noeud)
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: V1   !V(type,noeud) ,Vn-1
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: V2   !V(type,noeud) ,Vn
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: P
    REAL(KIND=8),                   INTENT(IN)  :: dt , alpha
    INTEGER                     ,   INTENT(IN)  :: mode
    INTEGER ::  m, l , i , ni , j
    REAL(KIND=8), DIMENSION(6) :: fs , fexp , ft , fp

    !fs pour le terme de forcage
    !fexp pour le terme de vitesse explicite
    !ft pour le terme temporelle   
    !fp pour le gradient de pression

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8)   :: ray

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    u0   = 0.d0
    fs   = 0.d0
    fexp = 0.d0
    ft   = 0.d0
    fp   = 0.d0
    DO m = 1, me
       DO l = 1, l_G

          !--------On calcul le rayon du point gauss
          ray = 0
          DO ni = 1, n_w;  i = jj(ni,m)
             ray = ray + mesh%rr(1,i)*ww(ni,l)
          END DO

          !--------calcul de la premiere composante : u0(1,:) <--> f(r,m,c)

          fs(1) = SUM(ff(1,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(1) = (-alpha*2*mode/ray**2)*SUM((2*V2(4,jj(:,m))-V1(4,jj(:,m)))* ww(:,l)) * rj(l,m)  
          ft(1) = 1/(2*dt) * SUM((4*V2(1,jj(:,m))-V1(1,jj(:,m))) * ww(:,l)) * rj(l,m)
          fp(1) = -SUM(P(1,jj(:,m))*dw(1,:,l,m)) * rj(l,m)

          !--------calcul de la seconde composante : u0(2,:) <--> f(r,m,s)

          fs(2) = SUM(ff(2,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(2) = (alpha*2*mode/ray**2)*SUM((2*V2(3,jj(:,m))-V1(3,jj(:,m)))* ww(:,l)) * rj(l,m)  
          ft(2) = 1/(2*dt) * SUM((4*V2(2,jj(:,m))-V1(2,jj(:,m))) * ww(:,l)) * rj(l,m)
          fp(2) = -SUM(P(2,jj(:,m))*dw(1,:,l,m)) * rj(l,m)

          !--------calcul de la troisieme composante : u0(3,:) <--> f(th,m,c)

          fs(3) = SUM(ff(3,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(3) = (alpha*2*mode/ray**2)*SUM((2*V2(2,jj(:,m))-V1(2,jj(:,m)))* ww(:,l)) * rj(l,m)  
          ft(3) = 1/(2*dt) * SUM((4*V2(3,jj(:,m))-V1(3,jj(:,m))) * ww(:,l)) * rj(l,m)
          fp(3) = -SUM(P(2,jj(:,m))*ww(:,l)) * rj(l,m)/ray*mode

          !--------calcul de la quatrieme composante : u0(4,:) <--> f(th,m,s)

          fs(4) = SUM(ff(4,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(4) = (-alpha*2*mode/ray**2)*SUM((2*V2(1,jj(:,m))-V1(1,jj(:,m)))* ww(:,l)) * rj(l,m)  
          ft(4) = 1/(2*dt) * SUM((4*V2(4,jj(:,m))-V1(4,jj(:,m))) * ww(:,l)) * rj(l,m)
          fp(4) = SUM(P(1,jj(:,m))*ww(:,l)) * rj(l,m)/ray *mode

          !--------calcul de la cinquieme composante : u0(5,:) <--> f(z,m,c)

          fs(5) = SUM(ff(5,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(5) = 0.d0
          ft(5) = 1/(2*dt) * SUM((4*V2(5,jj(:,m))-V1(5,jj(:,m))) * ww(:,l)) * rj(l,m)
          fp(5) = -SUM(P(1,jj(:,m))*dw(2,:,l,m)) * rj(l,m)

          !--------calcul de la sixieme composante : u0(6,:) <--> f(z,m,s)

          fs(6) = SUM(ff(6,jj(:,m)) * ww(:,l)) * rj(l,m)
          fexp(6) = 0.d0
          ft(6) = 1/(2*dt) * SUM((4*V2(6,jj(:,m))-V1(6,jj(:,m))) * ww(:,l)) * rj(l,m)
          fp(6) = -SUM(P(2,jj(:,m))*dw(2,:,l,m)) * rj(l,m)

          !---------- calcul du second membre total

          DO j=1,6
             DO ni = 1, n_w 
                u0(j,jj(ni,m)) = u0(j,jj(ni,m)) +  ww(ni,l)  * ray * (fs(j) + fexp(j) + ft(j) + fp(j)) 
             ENDDO
          ENDDO


       ENDDO
    ENDDO

  END SUBROUTINE qs_00_stokes_3d

  SUBROUTINE qs_00_stokes_3d_new_ssr(mesh,alpha,mode,ff,V1,V2,P,dt,u0)
    !=================================
    !sans le terme de couplage

    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: ff   !analytique (type,noeud)
    REAL(KIND=8), DIMENSION(:,:),   INTENT(OUT) :: u0   !u0(type,noeud)
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: V1   !V(type,noeud) ,Vn-1
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: V2   !V(type,noeud) ,Vn
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: P
    REAL(KIND=8),                   INTENT(IN)  :: dt , alpha
    INTEGER                     ,   INTENT(IN)  :: mode
    INTEGER ::  m, l , i , ni , j
    REAL(KIND=8), DIMENSION(6) :: fs , ft , fp

    !fs pour le terme de forcage
    !fexp pour le terme de vitesse explicite
    !ft pour le terme temporelle
    !fp pour le gradient de pression

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8)   :: ray

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    u0   = 0.d0
    fs   = 0.d0
    ft   = 0.d0
    fp   = 0.d0

    DO m = 1, me
       DO l = 1, l_G

          !--------On calcul le rayon du point gauss
          ray = 0
          DO ni = 1, n_w;  i = jj(ni,m)
             ray = ray + mesh%rr(1,i)*ww(ni,l)
          END DO

          !--------calcul de la premiere composante : u0(1,:) <--> f(r,m,c)

          fs(1) = SUM(ff(1,jj(:,m)) * ww(:,l)) * rj(l,m)
          ft(1) = 1/(2*dt) * SUM((4*V2(1,jj(:,m))-V1(1,jj(:,m))) * ww(:,l)) * rj(l,m)
          fp(1) = -SUM(P(1,jj(:,m))*dw(1,:,l,m)) * rj(l,m)

          !--------calcul de la seconde composante : u0(2,:) <--> f(r,m,s)

          fs(2) = SUM(ff(2,jj(:,m)) * ww(:,l)) * rj(l,m)
          ft(2) = 1/(2*dt) * SUM((4*V2(2,jj(:,m))-V1(2,jj(:,m))) * ww(:,l)) * rj(l,m)
          fp(2) = -SUM(P(2,jj(:,m))*dw(1,:,l,m)) * rj(l,m)

          !--------calcul de la troisieme composante : u0(3,:) <--> f(th,m,c)

          fs(3) = SUM(ff(3,jj(:,m)) * ww(:,l)) * rj(l,m)
          ft(3) = 1/(2*dt) * SUM((4*V2(3,jj(:,m))-V1(3,jj(:,m))) * ww(:,l)) * rj(l,m)
          fp(3) = -SUM(P(2,jj(:,m))*ww(:,l)) * rj(l,m)/ray*mode

          !--------calcul de la quatrieme composante : u0(4,:) <--> f(th,m,s)

          fs(4) = SUM(ff(4,jj(:,m)) * ww(:,l)) * rj(l,m)
          ft(4) = 1/(2*dt) * SUM((4*V2(4,jj(:,m))-V1(4,jj(:,m))) * ww(:,l)) * rj(l,m)
          fp(4) = SUM(P(1,jj(:,m))*ww(:,l)) * rj(l,m)/ray *mode

          !--------calcul de la cinquieme composante : u0(5,:) <--> f(z,m,c)

          fs(5) = SUM(ff(5,jj(:,m)) * ww(:,l)) * rj(l,m)
          ft(5) = 1/(2*dt) * SUM((4*V2(5,jj(:,m))-V1(5,jj(:,m))) * ww(:,l)) * rj(l,m)
          fp(5) = -SUM(P(1,jj(:,m))*dw(2,:,l,m)) * rj(l,m)

          !--------calcul de la sixieme composante : u0(6,:) <--> f(z,m,s)

          fs(6) = SUM(ff(6,jj(:,m)) * ww(:,l)) * rj(l,m)
          ft(6) = 1/(2*dt) * SUM((4*V2(6,jj(:,m))-V1(6,jj(:,m))) * ww(:,l)) * rj(l,m)
          fp(6) = -SUM(P(2,jj(:,m))*dw(2,:,l,m)) * rj(l,m)

          !---------- calcul du second membre total

          DO j=1,6
             DO ni = 1, n_w
                u0(j,jj(ni,m)) = u0(j,jj(ni,m)) +  ww(ni,l)   * (fs(j)  + ray*(ft(j) + fp(j)))
             ENDDO
          ENDDO


       ENDDO
    ENDDO

  END SUBROUTINE qs_00_stokes_3d_new_ssr


  SUBROUTINE qs_00_stokes_3d_new(mesh,alpha,mode,ff,V1,V2,P,dt,u0)
    !=================================
    !sans le terme de couplage


    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: ff   !analytique (type,noeud)
    REAL(KIND=8), DIMENSION(:,:),   INTENT(OUT) :: u0   !u0(type,noeud)
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: V1   !V(type,noeud) ,Vn-1
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: V2   !V(type,noeud) ,Vn
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: P
    REAL(KIND=8),                   INTENT(IN)  :: dt , alpha
    INTEGER                     ,   INTENT(IN)  :: mode
    INTEGER ::  m, l , i , ni , j
    REAL(KIND=8), DIMENSION(6) :: fs , ft , fp


    !fs pour le terme de forcage
    !fexp pour le terme de vitesse explicite
    !ft pour le terme temporelle   
    !fp pour le gradient de pression

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8)   :: ray

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    u0   = 0.d0
    fs   = 0.d0
    ft   = 0.d0
    fp   = 0.d0

    DO m = 1, me
       DO l = 1, l_G

          !--------On calcul le rayon du point gauss
          ray = 0
          DO ni = 1, n_w;  i = jj(ni,m)
             ray = ray + mesh%rr(1,i)*ww(ni,l)
          END DO

          !--------calcul de la premiere composante : u0(1,:) <--> f(r,m,c)

          fs(1) = SUM(ff(1,jj(:,m)) * ww(:,l)) * rj(l,m)
          ft(1) = 1/(2*dt) * SUM((4*V2(1,jj(:,m))-V1(1,jj(:,m))) * ww(:,l)) * rj(l,m)
          fp(1) = -SUM(P(1,jj(:,m))*dw(1,:,l,m)) * rj(l,m)

          !--------calcul de la seconde composante : u0(2,:) <--> f(r,m,s)

          fs(2) = SUM(ff(2,jj(:,m)) * ww(:,l)) * rj(l,m)
          ft(2) = 1/(2*dt) * SUM((4*V2(2,jj(:,m))-V1(2,jj(:,m))) * ww(:,l)) * rj(l,m)
          fp(2) = -SUM(P(2,jj(:,m))*dw(1,:,l,m)) * rj(l,m)

          !--------calcul de la troisieme composante : u0(3,:) <--> f(th,m,c)

          !       fs(3) = SUM(ff(3,jj(:,m)) * ww(:,l)) * rj(l,m)/ray
          fs(3) = SUM(ff(3,jj(:,m)) * ww(:,l)) * rj(l,m)
          ft(3) = 1/(2*dt) * SUM((4*V2(3,jj(:,m))-V1(3,jj(:,m))) * ww(:,l)) * rj(l,m)
          fp(3) = -SUM(P(2,jj(:,m))*ww(:,l)) * rj(l,m)/ray*mode

          !--------calcul de la quatrieme composante : u0(4,:) <--> f(th,m,s)

          !        fs(4) = SUM(ff(4,jj(:,m)) * ww(:,l)) * rj(l,m)/ray
          fs(4) = SUM(ff(4,jj(:,m)) * ww(:,l)) * rj(l,m)
          ft(4) = 1/(2*dt) * SUM((4*V2(4,jj(:,m))-V1(4,jj(:,m))) * ww(:,l)) * rj(l,m)
          fp(4) = SUM(P(1,jj(:,m))*ww(:,l)) * rj(l,m)/ray *mode

          !--------calcul de la cinquieme composante : u0(5,:) <--> f(z,m,c)

          fs(5) = SUM(ff(5,jj(:,m)) * ww(:,l)) * rj(l,m)
          ft(5) = 1/(2*dt) * SUM((4*V2(5,jj(:,m))-V1(5,jj(:,m))) * ww(:,l)) * rj(l,m)
          fp(5) = -SUM(P(1,jj(:,m))*dw(2,:,l,m)) * rj(l,m)

          !--------calcul de la sixieme composante : u0(6,:) <--> f(z,m,s)

          fs(6) = SUM(ff(6,jj(:,m)) * ww(:,l)) * rj(l,m)
          ft(6) = 1/(2*dt) * SUM((4*V2(6,jj(:,m))-V1(6,jj(:,m))) * ww(:,l)) * rj(l,m)
          fp(6) = -SUM(P(2,jj(:,m))*dw(2,:,l,m)) * rj(l,m)

          !---------- calcul du second membre total

          DO j=1,6
             DO ni = 1, n_w 
                u0(j,jj(ni,m)) = u0(j,jj(ni,m)) +  ww(ni,l)  * ray * (fs(j)  + ft(j) + fp(j)) 
             ENDDO
          ENDDO


       ENDDO
    ENDDO

  END SUBROUTINE qs_00_stokes_3d_new


  SUBROUTINE qs_00_NS_inline(mesh,mode,m_max,ff,V1,V2,P,dt,u0,meth,tps_sft)
    !=================================
    !Calcul du second membre pour le probleme de Navier_Stokes
    !avec terme non lineaire et source (ff)
    !Le terme non lineaire est calcule sur les points de Gauss 
    !directement dans cette procedure. Le calcul peut s'effectuer avec
    !Matmul ou avec des multiplications classiques


    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:),               INTENT(IN)   :: ff   !analytique (type,noeud)
    REAL(KIND=8), DIMENSION(:,:),               INTENT(OUT)  :: u0   !u0(type,noeud)
    TYPE(mesh_type), TARGET                                  :: mesh
    INTEGER                     ,               INTENT(IN)   :: m_max   
    REAL(KIND=8), DIMENSION(6,mesh%np,0:m_max), INTENT(IN)   :: V1   !V(type,noeud) ,Vn-1
    REAL(KIND=8), DIMENSION(6,mesh%np,0:m_max), INTENT(IN)   :: V2   !V(type,noeud) ,Vn   
    REAL(KIND=8), DIMENSION(:,:),               INTENT(IN)   :: P
    REAL(KIND=8),                               INTENT(IN)   :: dt 
    INTEGER                     ,               INTENT(IN)   :: mode 
    CHARACTER(len=64),                          INTENT(IN)   :: meth 
    REAL(KIND=8), OPTIONAL,                     INTENT(INOUT):: tps_sft   

    INTEGER                                                  :: m, l , i , ni , j
    REAL(KIND=8), DIMENSION(6)                               :: fs , ft , fp, fnl, smb
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:,:),SAVE       :: W, RotV, ProdI
    INTEGER                                                  :: N
    INTEGER, ALLOCATABLE,DIMENSION(:),            SAVE       :: j_loc
    REAL(KIND=8), ALLOCATABLE,DIMENSION(:,:),     SAVE       :: dw_loc 
    REAL(KIND=8), DIMENSION(6,mesh%np)                       :: V1m, V2m, Vt   
    REAL(KIND=8), DIMENSION(6,mesh%np,0:m_max)               :: Vs   
    LOGICAL,                                      SAVE       :: once = .TRUE.
    LOGICAL,                                      SAVE       :: once_sft = .TRUE.
    REAL(KIND=8):: user_time
    !EXTERNAL :: user_time
    REAL(KIND=8) :: tps, dummy, tt,tps1                     

    !fs pour le terme de forcage
    !ft pour le terme temporelle   
    !fp pour le gradient de pression
    !fnl pour les termes non lineaires

    !W1, W2 : vitesses dans le format accepte par fft_prod ; 
    !         W1(1,:)=Vr, W1(2,:)=Vth, W1(3,:)=Vz; 1<->tps n-1, 2<->tps n 
    !RotV1, RotV2 : pareil pour le rotationnel
    !Prod1, Prod2 : pareil pour le produit
    !Prod : produits intermediaire pour le calcul du produit vect

    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8)   :: ray


    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:), SAVE                 :: aps, asp        
    REAL(KIND=8)                                                    :: PI
    INTEGER                                                         :: li, ki   
    REAL(KIND=8),              DIMENSION(0:2*m_max)                 :: vr_p, vth_p, vz_p, rotr_p, rotth_p, rotz_p   
    REAL(KIND=8),              DIMENSION(3,0:2*m_max)               :: prod   

    !--------------------calcul du produit vectoriel----------------------------
    !{rot(V) x V}r  = {rotV}z*Vth - {rotV}th*Vz = prodr
    !{rot(V) x V}th = {rotV}r*Vz  - {rotV}z*Vr  = prodth
    !{rot(V) x V}z  = {rotV}th*Vr - {rotV}r*Vth = prodz
    !soit 9 sft
    !----------------------------------------------------------------------------

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me    

    IF (once) THEN

       once =.FALSE.
       !CALCUL DES MATRICES UTILES AUX SFT
       PI = ACOS(-1.d0)
       ALLOCATE(asp(0:2*m_max, 0:2*m_max))
       ALLOCATE(aps(0:2*m_max, 0:2*m_max))    

       !calcul de asp : spectral > physique

       DO li= 0, 2*m_max
          DO ki= 0, 2*m_max
             IF(ki < m_max+1) THEN
                asp(li, ki) = COS(li*2.d0*PI/(2*m_max+1)*ki)
             ELSE
                asp(li, ki) = SIN(li*2.d0*PI/(2*m_max+1)*(ki-m_max))
             ENDIF
          ENDDO
       ENDDO

       !calcul de aps : physique -> spectrale 

       DO ki= 0, 2*m_max
          DO li= 0, 2*m_max
             IF(ki < m_max+1) THEN
                IF (ki == 0) THEN
                   aps(ki,li) = 0.5d0
                ELSE
                   aps(ki, li) = COS(li*2.d0*PI/(2*m_max+1)*ki) 
                ENDIF
             ELSE
                aps(ki, li) = SIN(li*2.d0*PI/(2*m_max+1)*(ki-m_max))
             ENDIF
             aps(ki, li) = aps(ki, li) * 2.d0/(FLOAT(2*m_max)+1.d0)
          ENDDO
       ENDDO
       WRITE(*,*) 'calcul matrices sft done'

       ALLOCATE(W(3,0:2*m_max,l_G,me))
       ALLOCATE(RotV(3,0:2*m_max,l_G,me))
       ALLOCATE(ProdI(3,0:2*m_max,l_G,me))
       ALLOCATE(j_loc(mesh%gauss%n_w))   
       ALLOCATE(dw_loc(mesh%gauss%k_d,mesh%gauss%n_w))

    ENDIF


    IF (mode == 0) THEN  !pour le calcul des sft
       Vs = 2*V2-V1
    ENDIF

    V1m = V1(:,:,mode)
    V2m = V2(:,:,mode)      

    u0   = 0.d0

    tps = 0  
    tps1 = 0 
    DO m = 1, me
       j_loc = jj(:,m)
       DO l = 1, l_G

          !tt = user_time(dummy)
          dw_loc = dw(:,:,l,m)

          !--------On calcul le rayon du point gauss
          ray = SUM(mesh%rr(1,j_loc)*ww(:,l))

          DO j=1, 6
             fs(j) = SUM(ff(j,j_loc) * ww(:,l))
             ft(j) = 1/(2*dt)*SUM((4*V2m(j,j_loc)-V1m(j,j_loc)) * ww(:,l))
          ENDDO

          fp(1) = -SUM(P(1,j_loc)*dw_loc(1,:))
          fp(2) = -SUM(P(2,j_loc)*dw_loc(1,:))
          fp(3) = -SUM(P(2,j_loc)*ww(:,l))/ray*mode
          fp(4) = SUM(P(1,j_loc)*ww(:,l))/ray *mode 
          fp(5) = -SUM(P(1,j_loc)*dw_loc(2,:))
          fp(6) = -SUM(P(2,j_loc)*dw_loc(2,:)) 


          IF (mode == 0) THEN

             DO j = 0, m_max 

                !tt = user_time(dummy)
                !-----------------rotationnel sur les points de Gauss-------------------

                !coeff sur les cosinus 
                RotV(1,j,l,m) = j/ray*SUM(Vs(6,j_loc,j)*ww(:,l)) &
                     -SUM(Vs(3,j_loc,j)*dw_loc(2,:))
                RotV(2,j,l,m) = SUM(Vs(1,j_loc,j)*dw_loc(2,:)) &
                     -SUM(Vs(5,j_loc,j)*dw_loc(1,:))
                RotV(3,j,l,m) = 1/ray*SUM(Vs(3,j_loc,j)*ww(:,l))+SUM(Vs(3,j_loc,j) &
                     *dw_loc(1,:))-j/ray*SUM(Vs(2,j_loc,j)*ww(:,l))
                !coeff sur les sinus 
                IF (j /= 0) THEN         
                   RotV(1,j+m_max,l,m) = -j/ray*SUM(Vs(5,j_loc,j)*ww(:,l)) &
                        -SUM(Vs(4,j_loc,j)*dw_loc(2,:))
                   RotV(2,j+m_max,l,m) = SUM(Vs(2,j_loc,j)*dw_loc(2,:)) &
                        -SUM(Vs(6,j_loc,j)*dw_loc(1,:))
                   RotV(3,j+m_max,l,m)  = 1/ray*SUM(Vs(4,j_loc,j)*ww(:,l))+SUM(Vs(4,j_loc,j) &
                        *dw_loc(1,:))+j/ray*SUM(Vs(1,j_loc,j)*ww(:,l))
                ENDIF
                !-----------------vitesse sur les points de Gauss---------------------------
                W(1,j,l,m) = SUM(Vs(1,j_loc,j)*ww(:,l))
                W(2,j,l,m) = SUM(Vs(3,j_loc,j)*ww(:,l))
                W(3,j,l,m) = SUM(Vs(5,j_loc,j)*ww(:,l))
                IF (j /= 0) THEN
                   W(1,j+m_max,l,m) = SUM(Vs(2,j_loc,j)*ww(:,l))
                   W(2,j+m_max,l,m) = SUM(Vs(4,j_loc,j)*ww(:,l))
                   W(3,j+m_max,l,m) = SUM(Vs(6,j_loc,j)*ww(:,l))
                ENDIF

                !tps1 = tps1+user_time(dummy) - tt
             ENDDO

             tt = user_time(dummy)

             IF (meth == 'sum') THEN
                DO li=0, 2*m_max
                   vr_p(li)    = SUM(asp(li,:)*W(1,:,l,m))
                   vth_p(li)   = SUM(asp(li,:)*W(2,:,l,m))
                   vz_p(li)    = SUM(asp(li,:)*W(3,:,l,m))
                   rotr_p(li)  = SUM(asp(li,:)*RotV(1,:,l,m))
                   rotth_p(li) = SUM(asp(li,:)*RotV(2,:,l,m))
                   rotz_p(li)  = SUM(asp(li,:)*RotV(3,:,l,m))
                ENDDO
             ELSEIF (meth == 'matmul') THEN   
                vr_p    = MATMUL(asp,W(1,:,l,m) )
                vth_p   = MATMUL(asp,W(2,:,l,m))
                vz_p    = MATMUL(asp,W(3,:,l,m))

                rotr_p  = MATMUL(asp,RotV(1,:,l,m)) 
                rotth_p = MATMUL(asp,RotV(2,:,l,m))
                rotz_p  = MATMUL(asp,RotV(3,:,l,m))
             ENDIF

             prod(1,:) = rotz_p(:)  * vth_p(:)  -  rotth_p(:) * vz_p(:)
             prod(2,:) = rotr_p(:)  * vz_p(:)   -  rotz_p(:)  * vr_p(:)
             prod(3,:) = rotth_p(:) * vr_p(:)   -  rotr_p(:)  * vth_p(:)

             IF (meth == 'sum') THEN
                DO li=0, 2*m_max
                   ProdI(1,li,l,m) = SUM(aps(li,:)*prod(1,:))
                   ProdI(2,li,l,m) = SUM(aps(li,:)*prod(2,:))
                   ProdI(3,li,l,m) = SUM(aps(li,:)*prod(3,:))
                ENDDO
             ELSEIF (meth == 'matmul') THEN   
                ProdI(1,:,l,m)  = MATMUL(aps, prod(1,:))
                ProdI(2,:,l,m)  = MATMUL(aps, prod(2,:))
                ProdI(3,:,l,m)  = MATMUL(aps, prod(3,:))   
             ENDIF

             IF (PRESENT(tps_sft)) THEN             
                tps_sft = tps_sft + user_time(dummy) - tt 
             ENDIF

          ENDIF

          !-----------------assemblage du second membre du terme NL--------------

          fnl(1) = ProdI(1,mode,l,m)
          fnl(3) = ProdI(2,mode,l,m)
          fnl(5) = ProdI(3,mode,l,m)
          IF (mode /= 0) THEN
             fnl(2) = ProdI(1,mode+m_max,l,m)
             fnl(4) = ProdI(2,mode+m_max,l,m)
             fnl(6) = ProdI(3,mode+m_max,l,m)
          ENDIF

          smb = (fnl+ft+fp+fs)*ray*rj(l,m)      

          !---------- calcul du second membre total

          DO j=1,6
             DO ni = 1, n_w 
                u0(j,jj(ni,m)) = u0(j,jj(ni,m)) +  ww(ni,l)*smb(j)
             ENDDO
          ENDDO

       ENDDO
    ENDDO

    !write(*,*) 'temps d une boucle ss terme nl=',tps    
    !write(*,*) 'temps du calcul de V et RotV aux points de Gauss=',tps1 

  END SUBROUTINE qs_00_NS_inline


  SUBROUTINE qs_stokes(mesh,mode,ff,V1m,V2m,P,dt,u0)
    !=================================
    !Calcul du second membre pour le probleme de Stokes
    !sans terme non lineaire et source (ff)

    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:),               INTENT(IN)   :: ff   !analytique (type,noeud)
    REAL(KIND=8), DIMENSION(:,:),               INTENT(OUT)  :: u0   !u0(type,noeud)
    TYPE(mesh_type), TARGET                                  :: mesh  
    REAL(KIND=8), DIMENSION(6,mesh%np),         INTENT(IN)   :: V1m, V2m   !V(type,noeud)
    REAL(KIND=8), DIMENSION(:,:),               INTENT(IN)   :: P
    REAL(KIND=8),                               INTENT(IN)   :: dt 
    INTEGER                     ,               INTENT(IN)   :: mode  

    INTEGER                                                  ::  m, l , i , ni , j
    REAL(KIND=8), DIMENSION(6)                               :: fs , ft , fp, fnl, smb
    INTEGER                                                  :: N
    INTEGER, DIMENSION(mesh%gauss%n_w)                       :: j_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w)   :: dw_loc     
    LOGICAL,                                      SAVE       :: once = .TRUE.
    LOGICAL,                                      SAVE       :: once_sft = .TRUE.

    !fs pour le terme de forcage
    !ft pour le terme temporelle   
    !fp pour le gradient de pression

    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8)   :: ray

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    u0   = 0.d0

    DO m = 1, me
       j_loc = jj(:,m)
       DO l = 1, l_G
          dw_loc = dw(:,:,l,m)

          !--------On calcul le rayon du point gauss
          ray = SUM(mesh%rr(1,j_loc)*ww(:,l))

          !-calcul des seconds membres pr les termes de forçage, temporels et de gradP 
          !--------calcul de la premiere composante : u0(1,:) <--> f(r,m,c)

          fs(1) = SUM(ff(1,j_loc) * ww(:,l))
          ft(1) = 1/(2*dt) * SUM((4*V2m(1,j_loc)-V1m(1,j_loc)) * ww(:,l)) 
          fp(1) = -SUM(P(1,j_loc)*dw_loc(1,:))

          !--------calcul de la seconde composante : u0(2,:) <--> f(r,m,s)

          fs(2) = SUM(ff(2,j_loc) * ww(:,l)) !* rj(l,m)
          ft(2) = 1/(2*dt) * SUM((4*V2m(2,j_loc)-V1m(2,j_loc)) * ww(:,l)) !* rj(l,m)
          fp(2) = -SUM(P(2,j_loc)*dw_loc(1,:)) !* rj(l,m)

          !--------calcul de la troisieme composante : u0(3,:) <--> f(th,m,c)

          fs(3) = SUM(ff(3,j_loc) * ww(:,l)) !* rj(l,m)
          ft(3) = 1/(2*dt) * SUM((4*V2m(3,j_loc)-V1m(3,j_loc)) * ww(:,l)) !* rj(l,m)
          fp(3) = -SUM(P(2,j_loc)*ww(:,l))/ray*mode !* rj(l,m)

          !--------calcul de la quatrieme composante : u0(4,:) <--> f(th,m,s)

          fs(4) = SUM(ff(4,j_loc) * ww(:,l)) !* rj(l,m)
          ft(4) = 1/(2*dt) * SUM((4*V2m(4,j_loc)-V1m(4,j_loc)) * ww(:,l)) !* rj(l,m)
          fp(4) = SUM(P(1,j_loc)*ww(:,l))/ray *mode !* rj(l,m)

          !--------calcul de la cinquieme composante : u0(5,:) <--> f(z,m,c)

          fs(5) = SUM(ff(5,j_loc) * ww(:,l)) !* rj(l,m)
          ft(5) = 1/(2*dt) * SUM((4*V2m(5,j_loc)-V1m(5,j_loc)) * ww(:,l)) !* rj(l,m)
          fp(5) = -SUM(P(1,j_loc)*dw_loc(2,:)) !* rj(l,m)

          !--------calcul de la sixieme composante : u0(6,:) <--> f(z,m,s)

          fs(6) = SUM(ff(6,j_loc) * ww(:,l)) !* rj(l,m)
          ft(6) = 1/(2*dt) * SUM((4*V2m(6,j_loc)-V1m(6,j_loc)) * ww(:,l)) !* rj(l,m)
          fp(6) = -SUM(P(2,j_loc)*dw_loc(2,:)) !* rj(l,m)

          !-------calcul du second membre pour le terme nonlineaire------------------------
          !-------------------------------------------------------------------------------

          smb =  (ft+fp+fs)*ray*rj(l,m)      
          !write(95,*) m,l,ft(3)
          !---------- calcul du second membre total
          !smb = 1.d0
          DO j=1,6
             DO ni = 1, n_w 
                u0(j,jj(ni,m)) = u0(j,jj(ni,m)) +  ww(ni,l)*smb(j)
             ENDDO
          ENDDO

       ENDDO
    ENDDO


  END SUBROUTINE qs_stokes

  SUBROUTINE qs_navier_stokes_2006(mesh,mode,ff,V1m,P,dt,u0,rotv_v)
    !=================================
    !Calcul du second membre pour le probleme de Stokes
    !sans terme non lineaire et source (ff)

    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:),               INTENT(IN) :: ff, rotv_v   !analytique (type,noeud)
    REAL(KIND=8), DIMENSION(:,:),               INTENT(OUT):: u0   !u0(type,noeud)
    TYPE(mesh_type), TARGET                                :: mesh  
    REAL(KIND=8), DIMENSION(:,:),               INTENT(IN) :: V1m   !V(type,noeud) 
    REAL(KIND=8), DIMENSION(:,:),               INTENT(IN) :: P
    REAL(KIND=8),                               INTENT(IN) :: dt 
    INTEGER,                                    INTENT(IN) :: mode  
    INTEGER                                                :: m, l , i , ni , j, index
    REAL(KIND=8), DIMENSION(6)                             :: fs , ft , fp, smb
    INTEGER, DIMENSION(mesh%gauss%n_w)                     :: j_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w) :: dw_loc     

    !fs pour le terme de forcage
    !ft pour le terme temporelle   
    !fp pour le gradient de pression

    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8)   :: ray

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    u0   = 0.d0


    index = 0
    DO m = 1, me
       j_loc = jj(:,m)
       DO l = 1, l_G
          index  = index +1
          dw_loc = dw(:,:,l,m)

          !--------On calcul le rayon du point gauss
          ray = SUM(mesh%rr(1,j_loc)*ww(:,l))

          !-calcul des seconds membres pr les termes de forçage, temporels et de gradP 
          !--------calcul de la premiere composante : u0(1,:) <--> f(r,m,c)

          fs(1) = SUM(ff(j_loc,1) * ww(:,l))
          ft(1) = SUM(V1m(j_loc,1) * ww(:,l)) 
          fp(1) = -SUM(P(j_loc,1)*dw_loc(1,:))

          !--------calcul de la seconde composante : u0(2,:) <--> f(r,m,s)

          fs(2) = SUM(ff(j_loc,2) * ww(:,l)) !* rj(l,m)
          ft(2) = SUM(V1m(j_loc,2) * ww(:,l)) !* rj(l,m)
          fp(2) = -SUM(P(j_loc,2)*dw_loc(1,:)) !* rj(l,m)

          !--------calcul de la troisieme composante : u0(3,:) <--> f(th,m,c)

          fs(3) = SUM(ff(j_loc,3) * ww(:,l)) !* rj(l,m)
          ft(3) = SUM(V1m(j_loc,3) * ww(:,l)) !* rj(l,m)
          fp(3) = -SUM(P(j_loc,2)*ww(:,l))/ray*mode !* rj(l,m)

          !--------calcul de la quatrieme composante : u0(4,:) <--> f(th,m,s)

          fs(4) = SUM(ff(j_loc,4) * ww(:,l)) !* rj(l,m)
          ft(4) = SUM(V1m(j_loc,4) * ww(:,l)) !* rj(l,m)
          fp(4) = SUM(P(j_loc,1)*ww(:,l))/ray *mode !* rj(l,m)

          !--------calcul de la cinquieme composante : u0(5,:) <--> f(z,m,c)

          fs(5) = SUM(ff(j_loc,5) * ww(:,l)) !* rj(l,m)
          ft(5) = SUM(V1m(j_loc,5) * ww(:,l)) !* rj(l,m)
          fp(5) = -SUM(P(j_loc,1)*dw_loc(2,:)) !* rj(l,m)

          !--------calcul de la sixieme composante : u0(6,:) <--> f(z,m,s)

          fs(6) = SUM(ff(j_loc,6) * ww(:,l)) !* rj(l,m)
          ft(6) = SUM(V1m(j_loc,6) * ww(:,l)) !* rj(l,m)
          fp(6) = -SUM(P(j_loc,2)*dw_loc(2,:)) !* rj(l,m)

          !-------calcul du second membre pour le terme nonlineaire------------------------
          !-------------------------------------------------------------------------------

          smb =  (ft+fp+fs-rotv_v(:,index))*ray*rj(l,m)      

          !---------- calcul du second membre total
          DO j=1,6
             DO ni = 1, n_w 
                u0(jj(ni,m),j) = u0(jj(ni,m),j) +  ww(ni,l)*smb(j)
             ENDDO
          ENDDO

       ENDDO
    ENDDO


  END SUBROUTINE qs_navier_stokes_2006

  SUBROUTINE qs_ns_2006(mesh,mode,ff,V1m,P,dt,u0,rotv_v)
    !=================================
    !Calcul du second membre pour le probleme de Stokes
    !sans terme non lineaire et source (ff)

    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:),               INTENT(IN) :: ff, rotv_v   !analytique (type,noeud)
    REAL(KIND=8), DIMENSION(:,:),               INTENT(OUT):: u0   !u0(type,noeud)
    TYPE(mesh_type), TARGET                                :: mesh  
    REAL(KIND=8), DIMENSION(:,:),               INTENT(IN) :: V1m   !V(type,noeud) 
    REAL(KIND=8), DIMENSION(:,:),               INTENT(IN) :: P
    REAL(KIND=8),                               INTENT(IN) :: dt 
    INTEGER,                                    INTENT(IN) :: mode  
    INTEGER                                                :: m, l , i , ni , j, index
    REAL(KIND=8), DIMENSION(6)                             :: fs , ft , fp, smb
    INTEGER, DIMENSION(mesh%gauss%n_w)                     :: j_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w) :: dw_loc     

    !fs pour le terme de forcage
    !ft pour le terme temporelle   
    !fp pour le gradient de pression

    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8)   :: ray

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    u0   = 0.d0


    index = 0
    DO m = 1, me
       j_loc = jj(:,m)
       DO l = 1, l_G
          index  = index +1
          dw_loc = dw(:,:,l,m)

          !--------On calcul le rayon du point gauss
          ray = SUM(mesh%rr(1,j_loc)*ww(:,l))

          !-calcul des seconds membres pr les termes de forçage, temporels et de gradP 
          !--------calcul de la premiere composante : u0(1,:) <--> f(r,m,c)

          fs(1) = SUM(ff(j_loc,1) * ww(:,l))
          ft(1) = SUM(V1m(j_loc,1) * ww(:,l)) 
          fp(1) = -SUM(P(j_loc,1)*dw_loc(1,:))

          !--------calcul de la seconde composante : u0(2,:) <--> f(r,m,s)

          fs(2) = SUM(ff(j_loc,2) * ww(:,l)) !* rj(l,m)
          ft(2) = SUM(V1m(j_loc,2) * ww(:,l)) !* rj(l,m)
          fp(2) = -SUM(P(j_loc,2)*dw_loc(1,:)) !* rj(l,m)

          !--------calcul de la troisieme composante : u0(3,:) <--> f(th,m,c)

          fs(3) = SUM(ff(j_loc,3) * ww(:,l)) !* rj(l,m)
          ft(3) = SUM(V1m(j_loc,3) * ww(:,l)) !* rj(l,m)
          fp(3) = -SUM(P(j_loc,2)*ww(:,l))/ray*mode !* rj(l,m)

          !--------calcul de la quatrieme composante : u0(4,:) <--> f(th,m,s)

          fs(4) = SUM(ff(j_loc,4) * ww(:,l)) !* rj(l,m)
          ft(4) = SUM(V1m(j_loc,4) * ww(:,l)) !* rj(l,m)
          fp(4) = SUM(P(j_loc,1)*ww(:,l))/ray *mode !* rj(l,m)

          !--------calcul de la cinquieme composante : u0(5,:) <--> f(z,m,c)

          fs(5) = SUM(ff(j_loc,5) * ww(:,l)) !* rj(l,m)
          ft(5) = SUM(V1m(j_loc,5) * ww(:,l)) !* rj(l,m)
          fp(5) = -SUM(P(j_loc,1)*dw_loc(2,:)) !* rj(l,m)

          !--------calcul de la sixieme composante : u0(6,:) <--> f(z,m,s)

          fs(6) = SUM(ff(j_loc,6) * ww(:,l)) !* rj(l,m)
          ft(6) = SUM(V1m(j_loc,6) * ww(:,l)) !* rj(l,m)
          fp(6) = -SUM(P(j_loc,2)*dw_loc(2,:)) !* rj(l,m)

          !-------calcul du second membre pour le terme nonlineaire------------------------
          !-------------------------------------------------------------------------------

          smb =  (ft+fp+fs-rotv_v(index,:))*ray*rj(l,m)      

          !---------- calcul du second membre total
          DO j=1,6
             DO ni = 1, n_w 
                u0(jj(ni,m),j) = u0(jj(ni,m),j) +  ww(ni,l)*smb(j)
             ENDDO
          ENDDO

       ENDDO
    ENDDO


  END SUBROUTINE qs_ns_2006

  SUBROUTINE qs_ns_stab_new(mesh,mode,ff,vel_tot,V1m,vit,P,dudt,phalf,nlhalf,dt,u0,rotv_v)
    !=================================
    !Calcul du second membre pour le probleme de Stokes
    !sans terme non lineaire et source (ff)

    USE Gauss_points
    USE chaine_caractere
    USE boundary !To have access to test_de_convergence
    USE sub_plot
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:),               INTENT(IN) :: ff, rotv_v, dudt, nlhalf 
    REAL(KIND=8), DIMENSION(:),                 INTENT(IN) :: vel_tot      !(noeud)
    REAL(KIND=8), DIMENSION(:,:),               INTENT(OUT):: u0           !u0(noeud, type)
    TYPE(mesh_type), TARGET                                :: mesh  
    REAL(KIND=8), DIMENSION(:,:),               INTENT(IN) :: V1m, vit     !V(noeud, type) 
    REAL(KIND=8), DIMENSION(:,:),               INTENT(IN) :: P, phalf
    REAL(KIND=8),                               INTENT(IN) :: dt 
    INTEGER,                                    INTENT(IN) :: mode  
    INTEGER                                                :: m, l , i , ni , j, index, index2, type, k
    REAL(KIND=8), DIMENSION(6)                             :: fs , ft , fp, smb, fv, mult
    REAL(KIND=8), DIMENSION(2,6,mesh%gauss%l_G)            :: grad
    REAL(KIND=8), DIMENSION(6,mesh%gauss%l_G)              :: vitloc
    INTEGER, DIMENSION(mesh%gauss%n_w)                     :: j_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w) :: dw_loc  
    REAL(KIND=8), DIMENSION(mesh%me) :: visc_plot
    !fs pour le terme de forcage
    !ft pour le terme temporelle   
    !fp pour le gradient de pression

    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8), DIMENSION(6)   :: visc1, visc2
    REAL(KIND=8)   :: ray, div1, div2, hloc, cfl, vloc, normal_vit
    REAL(KIND=8), SAVE ::  coeff1, coeff2, coeff_ed_st, R_eff, visc_eff, surf, nu_loc, coeff_visc_ordre_un
    LOGICAL, SAVE :: once = .true.
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE :: h

    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%l_Gs, 2) :: dwni_loc
    INTEGER, DIMENSION(mesh%gauss%n_w,2) :: jji_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,2) :: u0loci, uloci
    REAL(KIND=8), DIMENSION(6,mesh%me) :: viscosity
    REAL(KIND=8), DIMENSION(6) :: norm_vit
    REAL(KIND=8) :: dul, h2, type_fe, ed_st, coeff3, coeff4
    INTEGER :: ms, ls, cotei
    LOGICAL, SAVE :: LES

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    IF (once) THEN
       once =.FALSE.
       IF (test_de_convergence) THEN
          coeff1=0.d0
          coeff2=0.d0
          coeff3=0.d0
          coeff4=0.d0
       ELSE
          OPEN(UNIT = 21, FILE = 'data', FORM = 'formatted', STATUS = 'unknown')
          CALL read_until(21, 'data_stab_LES_NS')
          READ(21,*) LES
          IF (LES) THEN
             READ(21,*) coeff1 ! Used to help work with large cfl 
             READ(21,*) coeff2 ! Ratio of LES viscosity    (infty=first-order, 0=zero LES)
             READ(21,*) coeff3 ! Ratio of edge stabilization (0=no edge stabilization)
             READ(21,*) coeff4 ! Ratio of first-order viscosity
          ELSE
             coeff1=0.d0
             coeff2=0.d0
             coeff3=0.d0
             coeff4=0.d0
          END IF
          CLOSE(21)
       END IF

       surf = 0.d0
       DO m = 1, me
          DO l = 1, l_G
             ray = SUM(mesh%rr(1,mesh%jj(:,m))*mesh%gauss%ww(:,l))
             surf = surf + rj(l,m)*ray    
          END DO
       END DO

       IF (mesh%gauss%n_w==3) THEN
          type_fe = 1
       ELSE 
          type_fe = 2 
       END IF
       coeff_ed_st         = coeff3*0.02d0/type_fe
       coeff_visc_ordre_un = coeff4*0.25d0
       ALLOCATE(h(mesh%mi))
       DO ms = 1, mesh%mi
          h(ms)  = SUM(mesh%gauss%rji(:,ms))/type_fe
       END DO
    END IF

    !ATTENTION: JLG Jan 25 2010
    !ATTENTION: coeff1 is assumed to be of the order of the convective velocity
    !that simplifies the semi-implicit treatment of the LES viscosity
    normal_vit = MAXVAL(vel_tot)
    DO type = 1, 6
       norm_vit(type) = SUM(ABS(vit(:,type)))/mesh%np + 1.d-14
    END DO
    !ATTENTION: JLG Jan 25 2010
    u0   = 0.d0
    R_eff = 0.d0
    cfl = 0
    index = 0
    index2 = 0
    IF (LES) THEN
       DO m = 1, me
          j_loc = jj(:,m)
          hloc = SQRT(SUM(rj(:,m)))
          vloc = maxval(vel_tot(j_loc))
          cfl = max(vloc*dt/hloc,cfl)
          visc1 = 0
          visc2 = 0
          DO l = 1, l_G
             index2  = index2 +1
             dw_loc = dw(:,:,l,m)

             !--------On calcule le rayon du point gauss
             ray = SUM(mesh%rr(1,j_loc)*ww(:,l))

             !--------calcul de la premiere composante 
             ft(1) = SUM(dudt(j_loc,1) * ww(:,l))
             fp(1) = SUM(phalf(j_loc,1)*dw_loc(1,:))
             !--------calcul de la seconde composante 
             ft(2) = SUM(dudt(j_loc,2) * ww(:,l)) 
             fp(2) = SUM(phalf(j_loc,2)*dw_loc(1,:)) 
             !--------calcul de la troisieme composante 
             ft(3) = SUM(dudt(j_loc,3) * ww(:,l)) 
             fp(3) = SUM(phalf(j_loc,2)*ww(:,l))*mode/ray 
             !--------calcul de la quatrieme composante 
             ft(4) = SUM(dudt(j_loc,4) * ww(:,l)) 
             fp(4) = -SUM(phalf(j_loc,1)*ww(:,l))*mode/ray 
             !--------calcul de la cinquieme composante 
             ft(5) = SUM(dudt(j_loc,5) * ww(:,l)) 
             fp(5) = SUM(phalf(j_loc,1)*dw_loc(2,:)) 
             !--------calcul de la sixieme composante 
             ft(6) = SUM(dudt(j_loc,6) * ww(:,l)) 
             fp(6) = SUM(phalf(j_loc,2)*dw_loc(2,:)) 
             !-------calcul du second membre pour le terme nonlineaire------------------------

             visc1= MAX(visc1,ABS(ft+fp+nlhalf(index2,:)))

             !--------Calcul du gradient de la vitesse
             DO type = 1, 6
                DO k = 1 ,2
                   grad(k,type,l) = SUM(vit(j_loc,type)*dw_loc(k,:))
                END DO
                vitloc(type,l) = SUM(vit(j_loc,type)*ww(:,l))
             END DO

             !--------Calcul de la divergence
             div1 = abs(grad(1,1,l) + vitloc(1,l)/ray + grad(2,5,l) + vitloc(4,l)*mode/ray)
             div2 = abs(grad(1,2,l) + vitloc(2,l)/ray + grad(2,6,l) - vitloc(3,l)*mode/ray)
             visc2(1) = max(visc2(1),div1)
             visc2(2) = max(visc2(1),div2)
          END DO
          visc2(4) = visc2(1); visc2(5) = visc2(1) 
          visc2(3) = visc2(2); visc2(6) = visc2(2)

          nu_loc = 0.d0
          DO type = 1, 6
             !visc1(type) = MAX(visc1(type)/normal_vit,visc2(type))
             !visc1(type) = MAX(visc1(type),visc2(type)*normal_vit)/MAXVAL(ABS(vitloc(type,:)))
             visc1(type) = MAX(visc1(type),2*visc2(type)*normal_vit)/norm_vit(type)
             visc_eff = hloc*MIN(coeff_visc_ordre_un*normal_vit,coeff2*hloc*visc1(type))
             nu_loc = nu_loc + visc_eff
             !======Semi-implicit version==========
             viscosity(type,m) = coeff1*hloc - visc_eff
             !======Semi-implicit version==========
          END DO
          R_eff = R_eff + (nu_loc + dt**2*hloc)*hloc**2*ray 
          visc_plot(m) = (nu_loc/6)/(coeff_visc_ordre_un*hloc*normal_vit)
       END DO
    ELSE
       viscosity = 0.d0
       DO m = 1, me
          j_loc = jj(:,m)
          hloc = SQRT(SUM(rj(:,m)))
          vloc = maxval(vel_tot(j_loc))
          cfl = max(vloc*dt/hloc,cfl)
       END DO
    END IF
    !DO type = 1, 6
    !   CALL average(mesh,viscosity(type,:))
    !END DO

    DO m = 1, me
       mult(:)= viscosity(:,m)
       j_loc = jj(:,m)
       DO l = 1, l_G
          index  = index +1
          dw_loc = dw(:,:,l,m)
          DO type = 1, 6
             DO k = 1 ,2
                grad(k,type,l) = SUM(vit(j_loc,type)*dw_loc(k,:))
             END DO
             vitloc(type,l) = SUM(vit(j_loc,type)*ww(:,l))
          END DO

          !--------On calcule le rayon du point gauss
          ray = SUM(mesh%rr(1,j_loc)*ww(:,l))

          !-calcul des seconds membres pr les termes de forçage, temporels et de gradP 
          !--------calcul de la premiere composante : u0(1,:) <--> f(r,m,c)
          fs(1) = SUM(ff(j_loc,1) * ww(:,l))
          ft(1) = SUM(V1m(j_loc,1) * ww(:,l)) 
          fp(1) = -SUM(P(j_loc,1)*dw_loc(1,:))
          fv(1) = ((mode*vitloc(1,l)+vitloc(4,l))*mode +mode*vitloc(4,l)+vitloc(1,l))/ray**2
          !--------calcul de la seconde composante : u0(2,:) <--> f(r,m,s)
          fs(2) = SUM(ff(j_loc,2) * ww(:,l)) 
          ft(2) = SUM(V1m(j_loc,2) * ww(:,l)) 
          fp(2) = -SUM(P(j_loc,2)*dw_loc(1,:)) 
          fv(2) = ((mode*vitloc(2,l)-vitloc(3,l))*mode -mode*vitloc(3,l)+vitloc(2,l))/ray**2
          !--------calcul de la troisieme composante : u0(3,:) <--> f(th,m,c)
          fs(3) = SUM(ff(j_loc,3) * ww(:,l)) 
          ft(3) = SUM(V1m(j_loc,3) * ww(:,l)) 
          fp(3) = -SUM(P(j_loc,2)*ww(:,l))/ray*mode 
          fv(3) = (-mode*vitloc(2,l)+vitloc(3,l) +(mode*vitloc(3,l)-vitloc(2,l))*mode)/ray**2
          !--------calcul de la quatrieme composante : u0(4,:) <--> f(th,m,s)
          fs(4) = SUM(ff(j_loc,4) * ww(:,l)) 
          ft(4) = SUM(V1m(j_loc,4) * ww(:,l)) 
          fp(4) = SUM(P(j_loc,1)*ww(:,l))/ray *mode 
          fv(4) =  (mode*vitloc(1,l)+vitloc(4,l) +(mode*vitloc(4,l)+vitloc(1,l))*mode)/ray**2
          !--------calcul de la cinquieme composante : u0(5,:) <--> f(z,m,c)
          fs(5) = SUM(ff(j_loc,5) * ww(:,l)) 
          ft(5) = SUM(V1m(j_loc,5) * ww(:,l)) 
          fp(5) = -SUM(P(j_loc,1)*dw_loc(2,:)) 
          fv(5) =  vitloc(5,l)*(mode/ray)**2
          !--------calcul de la sixieme composante : u0(6,:) <--> f(z,m,s)
          fs(6) = SUM(ff(j_loc,6) * ww(:,l)) 
          ft(6) = SUM(V1m(j_loc,6) * ww(:,l)) 
          fp(6) = -SUM(P(j_loc,2)*dw_loc(2,:)) 
          fv(6) = vitloc(6,l)*(mode/ray)**2
          !-------calcul du second membre pour le terme nonlineaire------------------------
          !-------------------------------------------------------------------------------

          fv = mult*fv

          smb =  (ft+fp+fs+fv-rotv_v(index,:))*ray*rj(l,m)
          DO type = 1, 6
             grad(:,type,l) =  mult(type)*grad(:,type,l)*ray*rj(l,m)
          END DO
          !---------- calcul du second membre total
          DO j=1,6
             DO ni = 1, n_w 
                u0(j_loc(ni),j) = u0(j_loc(ni),j) +  ww(ni,l)*smb(j) + SUM(dw_loc(:,ni)*grad(:,j,l))
             ENDDO
          ENDDO

       ENDDO
    ENDDO


    WRITE(*,*) ' CFL = ', cfl, ' LES ', LES
    IF (LES) THEN
       WRITE(*,*) ' R_eff for mode', mode, normal_vit*6*surf/R_eff
    END IF
    !IF (mode==0) THEN 
    !CALL plot_const_p1_label(mesh%jj, mesh%rr, visc_plot, 'tttt_0.plt')
    !END IF

    IF (LES) THEN
       ed_st = normal_vit*coeff_ed_st
       DO ms = 1, mesh%mi
          dwni_loc = mesh%gauss%dwni(:,:,:,ms)
          jji_loc = mesh%jji(:,:,ms)
          h2 = -ed_st*h(ms)**2

          j_loc(1:n_ws) = mesh%jjsi(1:n_ws,ms)
          DO type = 1, 6
             DO cotei = 1, 2
                uloci(:,cotei) = vit(jji_loc(:,cotei),type)
             END DO

             u0loci = 0.d0
             DO ls = 1, mesh%gauss%l_Gs
                !--------On calcule le rayon du point gauss au milieu de la face
                ray = mesh%rr(1,j_loc(3))

                !--------Calcul du saut de la derivee normale de la vitesse
                dul = SUM(dwni_loc(:,ls,:)*uloci)*mesh%gauss%rji(ls,ms)*h2*ray
                DO cotei = 1, 2
                   DO ni = 1, mesh%gauss%n_w
                      u0loci(ni, cotei) =  u0loci(ni, cotei) + dwni_loc(ni,ls,cotei)*dul
                   END DO
                END DO
             END DO

             DO cotei = 1, 2
                DO ni = 1, mesh%gauss%n_w
                   u0(jji_loc(ni,cotei),type) = u0(jji_loc(ni,cotei),type) + u0loci(ni,cotei) 
                END DO
             END DO
          END DO

       END DO
    END IF

  END SUBROUTINE qs_ns_stab_new

  SUBROUTINE qs_ns_stab_2010(mesh, pp_mesh, mode,ff,vel_tot,V1m,vit,P,dudt,phalf,nlhalf,dt,u0,rotv_v)
    !=================================
    !Calcul du second membre pour le probleme de Stokes
    !sans terme non lineaire et source (ff)

    USE Gauss_points
    USE chaine_caractere
    USE boundary !To have access to test_de_convergence
    USE sub_plot
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:),               INTENT(IN) :: ff, rotv_v, dudt, nlhalf 
    REAL(KIND=8), DIMENSION(:),                 INTENT(IN) :: vel_tot      !(noeud)
    REAL(KIND=8), DIMENSION(:,:),               INTENT(OUT):: u0           !u0(noeud, type)
    TYPE(mesh_type), TARGET                                :: mesh, pp_mesh  
    REAL(KIND=8), DIMENSION(:,:),               INTENT(IN) :: V1m, vit     !V(noeud, type) 
    REAL(KIND=8), DIMENSION(:,:),               INTENT(IN) :: P, phalf
    REAL(KIND=8),                               INTENT(IN) :: dt 
    INTEGER,                                    INTENT(IN) :: mode  
    INTEGER                                                :: m, l , i , ni , j, index, index2, type, k
    REAL(KIND=8), DIMENSION(6)                             :: fs , ft , fp, smb, fv, mult
    REAL(KIND=8), DIMENSION(2,6,mesh%gauss%l_G)            :: grad
    REAL(KIND=8), DIMENSION(6,mesh%gauss%l_G)              :: vitloc
    INTEGER, DIMENSION(mesh%gauss%n_w)                     :: j_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w) :: dw_loc  
    REAL(KIND=8), DIMENSION(mesh%me) :: visc_plot
    !fs pour le terme de forcage
    !ft pour le terme temporelle   
    !fp pour le gradient de pression

    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8), DIMENSION(6)   :: visc1, visc2
    REAL(KIND=8)   :: ray, div1, div2, hloc, cfl, vloc, normal_vit
    REAL(KIND=8), SAVE ::  coeff1, coeff2, coeff_ed_st, R_eff, visc_eff, surf, nu_loc
    LOGICAL, SAVE :: once = .true.
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE :: h

    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%l_Gs, 2) :: dwni_loc
    INTEGER, DIMENSION(mesh%gauss%n_w,2) :: jji_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,2) :: u0loci, uloci
    REAL(KIND=8), DIMENSION(6,mesh%me)    :: viscosity
    REAL(KIND=8), DIMENSION(6,pp_mesh%np) :: viscosity_clement
    REAL(KIND=8), DIMENSION(mesh%np,6)    :: visc_p2
    REAL(KIND=8), DIMENSION(6)            :: norm_vit
    REAL(KIND=8) :: dul, h2, type_fe, ed_st
    INTEGER :: ms, ls, cotei

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    IF (once) THEN
       once =.FALSE.
       IF (test_de_convergence) THEN
          coeff1=0.d0
          coeff2=0.d0
       ELSE
          OPEN(UNIT = 21, FILE = 'data', FORM = 'formatted', STATUS = 'unknown')
          CALL read_until(21, 'data_stab_LES_NS')
          READ(21,*) coeff1, coeff2
          CLOSE(21)
       END IF
       surf = SUM(rj(:,:))

       IF (mesh%gauss%n_w==3) THEN
          type_fe = 1
       ELSE 
          type_fe = 2 
       END IF
       coeff_ed_st = 0.05d0/type_fe

       ALLOCATE(h(mesh%mi))
       DO ms = 1, mesh%mi
          h(ms)  = SUM(mesh%gauss%rji(:,ms))/type_fe
       END DO
    END IF

    !ATTENTION: JLG Jan 25 2010
    !ATTENTION: coeff1 is assumed to be of the order of the convective velocity
    !that simplifies the semi-implicit treatment of the LES viscosity
    normal_vit = MAXVAL(vel_tot)
    DO type = 1, 6
       norm_vit(type) = SUM(ABS(vit(:,type)))/mesh%np + 1.d-14
    END DO
    !ATTENTION: JLG Jan 25 2010
    u0   = 0.d0
    R_eff = 0.d0
    cfl = 0
    index = 0
    index2 = 0
    DO m = 1, me
       j_loc = jj(:,m)
       hloc = SQRT(SUM(rj(:,m)))
       !vloc = MIN(maxval(vel_tot(j_loc)),1.d0)
       vloc = maxval(vel_tot(j_loc))
       cfl = max(vloc*dt/hloc,cfl)
       visc1 = 0
       visc2 = 0
       DO l = 1, l_G
          index2  = index2 +1
          dw_loc = dw(:,:,l,m)

          !--------On calcule le rayon du point gauss
          ray = SUM(mesh%rr(1,j_loc)*ww(:,l))

          !--------calcul de la premiere composante 
          ft(1) = SUM(dudt(j_loc,1) * ww(:,l))
          fp(1) = SUM(phalf(j_loc,1)*dw_loc(1,:))
          !--------calcul de la seconde composante 
          ft(2) = SUM(dudt(j_loc,2) * ww(:,l)) 
          fp(2) = SUM(phalf(j_loc,2)*dw_loc(1,:)) 
          !--------calcul de la troisieme composante 
          ft(3) = SUM(dudt(j_loc,3) * ww(:,l)) 
          fp(3) = SUM(phalf(j_loc,2)*ww(:,l))*mode/ray 
          !--------calcul de la quatrieme composante 
          ft(4) = SUM(dudt(j_loc,4) * ww(:,l)) 
          fp(4) = -SUM(phalf(j_loc,1)*ww(:,l))*mode/ray 
          !--------calcul de la cinquieme composante 
          ft(5) = SUM(dudt(j_loc,5) * ww(:,l)) 
          fp(5) = SUM(phalf(j_loc,1)*dw_loc(2,:)) 
          !--------calcul de la sixieme composante 
          ft(6) = SUM(dudt(j_loc,6) * ww(:,l)) 
          fp(6) = SUM(phalf(j_loc,2)*dw_loc(2,:)) 
          !-------calcul du second membre pour le terme nonlineaire------------------------

          visc1= MAX(visc1,ABS(ft+fp+nlhalf(index2,:)))

          !--------Calcul du gradient de la vitesse
          DO type = 1, 6
             DO k = 1 ,2
                grad(k,type,l) = SUM(vit(j_loc,type)*dw_loc(k,:))
             END DO
             vitloc(type,l) = SUM(vit(j_loc,type)*ww(:,l))
          END DO

          !--------Calcul de la divergence
          div1 = abs(grad(1,1,l) + vitloc(1,l)/ray + grad(2,5,l) + vitloc(4,l)*mode/ray)
          div2 = abs(grad(1,2,l) + vitloc(2,l)/ray + grad(2,6,l) - vitloc(3,l)*mode/ray)
          visc2(1) = div1; visc2(4) = div1; visc2(5) = div1
          visc2(2) = div2; visc2(3) = div2; visc2(6) = div2
       END DO

       nu_loc = 0.d0
       DO type = 1, 6
          !visc1(type) = MAX(visc1(type)/normal_vit,visc2(type))
          visc1(type) = MAX(visc1(type),visc2(type)*normal_vit)/norm_vit(type)
          !visc1(type) = visc1(type)/norm_vit(type)
          visc_eff = hloc*MIN(0.25d0*normal_vit,coeff2*hloc*visc1(type))
          nu_loc = nu_loc + visc_eff
          !======Semi-implicit version==========
          viscosity(type,m) = coeff1*hloc - visc_eff
          !visc1(type) = coeff1*hloc - visc_eff
          !======Semi-implicit version==========
          !======Totally explicit version=======
          !visc1(type) = -visc_eff 
          !======Totally explicit version=======
       END DO
       R_eff = R_eff + (nu_loc + dt**2*hloc)*hloc**2
       visc_plot(m) = (nu_loc/6)/(0.25*hloc*normal_vit)
    END DO

    DO type = 1, 6
       !CALL average(mesh,viscosity(type,:))
       CALL clement_c(pp_mesh,viscosity(type,:),viscosity_clement(type,:))
       CALL inject_clement(pp_mesh%jj, mesh%jj, viscosity_clement(type,:), visc_p2(:,type))
    END DO
!CALL plot_const_p1_label(mesh%jj, mesh%rr, viscosity(1,:), 'tttt_0.plt')
!CALL plot_scalar_field(pp_mesh%jj, pp_mesh%rr, viscosity_clement(1,:), 'v1_p1.plt')
!CALL plot_scalar_field(mesh%jj, mesh%rr, visc_p2(:,1)/(0.25*hloc*normal_vit), 'v1.plt')
!CALL plot_scalar_field(mesh%jj, mesh%rr, visc_p2(:,3)/(0.25*hloc*normal_vit), 'v3.plt')
!CALL plot_scalar_field(mesh%jj, mesh%rr, visc_p2(:,5)/(0.25*hloc*normal_vit), 'v5.plt')
!stop
    DO m = 1, me
       !mult(:)= viscosity(:,m)
       j_loc = jj(:,m)
       DO l = 1, l_G
          
          index  = index +1
          dw_loc = dw(:,:,l,m)
          DO type = 1, 6
             mult(type)= SUM(visc_p2(j_loc(:),type)*ww(:,l))
             DO k = 1 ,2
                grad(k,type,l) = SUM(vit(j_loc,type)*dw_loc(k,:))
             END DO
             vitloc(type,l) = SUM(vit(j_loc,type)*ww(:,l))
          END DO

          !--------On calcule le rayon du point gauss
          ray = SUM(mesh%rr(1,j_loc)*ww(:,l))

          !-calcul des seconds membres pr les termes de forçage, temporels et de gradP 
          !--------calcul de la premiere composante : u0(1,:) <--> f(r,m,c)
          fs(1) = SUM(ff(j_loc,1) * ww(:,l))
          ft(1) = SUM(V1m(j_loc,1) * ww(:,l)) 
          fp(1) = -SUM(P(j_loc,1)*dw_loc(1,:))
          fv(1) = ((mode*vitloc(1,l)+vitloc(4,l))*mode +mode*vitloc(4,l)+vitloc(1,l))/ray**2
          !--------calcul de la seconde composante : u0(2,:) <--> f(r,m,s)
          fs(2) = SUM(ff(j_loc,2) * ww(:,l)) 
          ft(2) = SUM(V1m(j_loc,2) * ww(:,l)) 
          fp(2) = -SUM(P(j_loc,2)*dw_loc(1,:)) 
          fv(2) = ((mode*vitloc(2,l)-vitloc(3,l))*mode -mode*vitloc(3,l)+vitloc(2,l))/ray**2
          !--------calcul de la troisieme composante : u0(3,:) <--> f(th,m,c)
          fs(3) = SUM(ff(j_loc,3) * ww(:,l)) 
          ft(3) = SUM(V1m(j_loc,3) * ww(:,l)) 
          fp(3) = -SUM(P(j_loc,2)*ww(:,l))/ray*mode 
          fv(3) = (-mode*vitloc(2,l)+vitloc(3,l) +(mode*vitloc(3,l)-vitloc(2,l))*mode)/ray**2
          !--------calcul de la quatrieme composante : u0(4,:) <--> f(th,m,s)
          fs(4) = SUM(ff(j_loc,4) * ww(:,l)) 
          ft(4) = SUM(V1m(j_loc,4) * ww(:,l)) 
          fp(4) = SUM(P(j_loc,1)*ww(:,l))/ray *mode 
          fv(4) =  (mode*vitloc(1,l)+vitloc(4,l) +(mode*vitloc(4,l)+vitloc(1,l))*mode)/ray**2
          !--------calcul de la cinquieme composante : u0(5,:) <--> f(z,m,c)
          fs(5) = SUM(ff(j_loc,5) * ww(:,l)) 
          ft(5) = SUM(V1m(j_loc,5) * ww(:,l)) 
          fp(5) = -SUM(P(j_loc,1)*dw_loc(2,:)) 
          fv(5) =  vitloc(5,l)*(mode/ray)**2
          !--------calcul de la sixieme composante : u0(6,:) <--> f(z,m,s)
          fs(6) = SUM(ff(j_loc,6) * ww(:,l)) 
          ft(6) = SUM(V1m(j_loc,6) * ww(:,l)) 
          fp(6) = -SUM(P(j_loc,2)*dw_loc(2,:)) 
          fv(6) = vitloc(6,l)*(mode/ray)**2
          !-------calcul du second membre pour le terme nonlineaire------------------------
          !-------------------------------------------------------------------------------

          !mult = visc1
          fv = mult*fv

          smb =  (ft+fp+fs+fv-rotv_v(index,:))*ray*rj(l,m)
          DO type = 1, 6
             grad(:,type,l) =  mult(type)*grad(:,type,l)*ray*rj(l,m)
          END DO
          !---------- calcul du second membre total
          DO j=1,6
             DO ni = 1, n_w 
                u0(j_loc(ni),j) = u0(j_loc(ni),j) +  ww(ni,l)*smb(j) + SUM(dw_loc(:,ni)*grad(:,j,l))
             ENDDO
          ENDDO

       ENDDO
    ENDDO

    WRITE(*,*) ' CFL = ', cfl, ' R_eff for mode', mode, 6*surf/R_eff
    !IF (mode==0) THEN 
    !   CALL plot_const_p1_label(mesh%jj, mesh%rr, visc_plot, 'tttt_0.plt')
    !END IF

    ed_st = normal_vit*coeff_ed_st
    DO ms = 1, mesh%mi
       dwni_loc = mesh%gauss%dwni(:,:,:,ms)
       jji_loc = mesh%jji(:,:,ms)
       h2 = -ed_st*h(ms)**2

       j_loc(1:n_ws) = mesh%jjsi(1:n_ws,ms)
       DO type = 1, 6
          DO cotei = 1, 2
             uloci(:,cotei) = vit(jji_loc(:,cotei),type)
          END DO

          u0loci = 0.d0
          DO ls = 1, mesh%gauss%l_Gs
             !--------On calcule le rayon du point gauss courant
             !ray = SUM(mesh%rr(1,j_loc(1:n_ws))*wws(:,ls))
             ray = mesh%rr(1,j_loc(3))

             !--------Calcul du saut de la derivee normale de la vitesse
             dul = SUM(dwni_loc(:,ls,:)*uloci)*mesh%gauss%rji(ls,ms)*h2*ray
             DO cotei = 1, 2
                DO ni = 1, mesh%gauss%n_w
                   u0loci(ni, cotei) =  u0loci(ni, cotei) + dwni_loc(ni,ls,cotei)*dul
                END DO
             END DO
          END DO

          DO cotei = 1, 2
             DO ni = 1, mesh%gauss%n_w
                u0(jji_loc(ni,cotei),type) = u0(jji_loc(ni,cotei),type) + u0loci(ni,cotei) 
             END DO
          END DO
       END DO

    END DO

  END SUBROUTINE qs_ns_stab_2010


  SUBROUTINE qs_ns_stab_2008(mesh,mode,ff,vel_tot,V1m,vit,P,dt,u0,rotv_v)
    !=================================
    !Calcul du second membre pour le probleme de Stokes
    !sans terme non lineaire et source (ff)

    USE Gauss_points
    USE chaine_caractere
    USE boundary !To have access to test_de_convergence
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:),               INTENT(IN) :: ff, rotv_v   !analytique (type,noeud)
    REAL(KIND=8), DIMENSION(:),                 INTENT(IN) :: vel_tot      !(noeud)
    REAL(KIND=8), DIMENSION(:,:),               INTENT(OUT):: u0           !u0(noeud, type)
    TYPE(mesh_type), TARGET                                :: mesh  
    REAL(KIND=8), DIMENSION(:,:),               INTENT(IN) :: V1m, vit     !V(noeud, type) 
    REAL(KIND=8), DIMENSION(:,:),               INTENT(IN) :: P
    REAL(KIND=8),                               INTENT(IN) :: dt 
    INTEGER,                                    INTENT(IN) :: mode  
    INTEGER                                                :: m, l , i , ni , j, index, type, k
    REAL(KIND=8), DIMENSION(6)                             :: fs , ft , fp, smb, fv, mult
    REAL(KIND=8), DIMENSION(2,6,mesh%gauss%l_G)            :: grad
    REAL(KIND=8), DIMENSION(6,mesh%gauss%l_G)              :: vitloc
    INTEGER, DIMENSION(mesh%gauss%n_w)                     :: j_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w) :: dw_loc     

    !fs pour le terme de forcage
    !ft pour le terme temporelle   
    !fp pour le gradient de pression

    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8)   :: ray, div1, div2, vit1=0.005d0, vit2=.075d0, visc1, visc2, hloc, cfl, vloc
    REAL(KIND=8), SAVE ::  coeff1, coeff2, R_eff, visc_eff
    LOGICAL, SAVE :: once = .true.

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    IF (once) THEN
       once =.FALSE.
       IF (test_de_convergence) THEN
          coeff1=0.d0
          coeff2=0.d0
       ELSE
          OPEN(UNIT = 21, FILE = 'data', FORM = 'formatted', STATUS = 'unknown')
          CALL read_until(21, 'data_stab_LES_NS')
          READ(21,*) coeff1, coeff2
          CLOSE(21)
       END IF
    END IF

    u0   = 0.d0
    R_eff = 0.d0
    cfl = 0
    index = 0
    DO m = 1, me
       j_loc = jj(:,m)
       hloc = SQRT(SUM(rj(:,m)))
       vloc = MIN(maxval(vel_tot(j_loc)),1.d0)
       !vloc = maxval(vel_tot(j_loc))
       cfl = max(vloc*dt/hloc,cfl)
       visc1 = 0
       visc2 = 0
       DO l = 1, l_G
          dw_loc = dw(:,:,l,m)

          !--------On calcule le rayon du point gauss
          ray = SUM(mesh%rr(1,j_loc)*ww(:,l))

          !--------Calcul du gradient de la vitesse
          DO type = 1, 6
             DO k = 1 ,2
                grad(k,type,l) = SUM(vit(j_loc,type)*dw_loc(k,:))
             END DO
             vitloc(type,l) = SUM(vit(j_loc,type)*ww(:,l))
          END DO

          !--------Calcul de la divergence
          div1 = grad(1,1,l) + vitloc(1,l)/ray + grad(2,5,l) + vitloc(4,l)*mode/ray
          div2 = grad(1,2,l) + vitloc(2,l)/ray + grad(2,6,l) - vitloc(3,l)*mode/ray
          visc1 = MAX(visc1,abs(div1))
          visc2 = MAX(visc2,abs(div2))
       END DO
       visc1 = MAX(visc1,visc2)
       !visc1 = 0.d0
       !visc_eff = hloc*MIN(coeff1,coeff2*hloc*visc1) - coeff1*hloc
       visc_eff = hloc*MIN(coeff1,coeff2*visc1)
       !======Semi-implicit version==========
       !visc1 = visc_eff - coeff1*hloc
       !======Semi-implicit version==========
       !======Totally explicit version=======
       visc1 = visc_eff 
       !======Totally explicit version=======
       R_eff = R_eff + visc_eff*SUM(rj(:,m))
       visc2 = visc1

       DO l = 1, l_G
          index  = index +1
          dw_loc = dw(:,:,l,m)

          !--------On calcule le rayon du point gauss
          ray = SUM(mesh%rr(1,j_loc)*ww(:,l))

          !-calcul des seconds membres pr les termes de forçage, temporels et de gradP 
          !--------calcul de la premiere composante : u0(1,:) <--> f(r,m,c)

          fs(1) = SUM(ff(j_loc,1) * ww(:,l))
          ft(1) = SUM(V1m(j_loc,1) * ww(:,l)) 
          fp(1) = -SUM(P(j_loc,1)*dw_loc(1,:))
          fv(1) = ((mode*vitloc(1,l)+vitloc(4,l))*mode +mode*vitloc(4,l)+vitloc(1,l))/ray**2
          !--------calcul de la seconde composante : u0(2,:) <--> f(r,m,s)

          fs(2) = SUM(ff(j_loc,2) * ww(:,l)) !* rj(l,m)
          ft(2) = SUM(V1m(j_loc,2) * ww(:,l)) !* rj(l,m)
          fp(2) = -SUM(P(j_loc,2)*dw_loc(1,:)) !* rj(l,m)
          fv(2) = ((mode*vitloc(2,l)-vitloc(3,l))*mode -mode*vitloc(3,l)+vitloc(2,l))/ray**2
          !--------calcul de la troisieme composante : u0(3,:) <--> f(th,m,c)

          fs(3) = SUM(ff(j_loc,3) * ww(:,l)) !* rj(l,m)
          ft(3) = SUM(V1m(j_loc,3) * ww(:,l)) !* rj(l,m)
          fp(3) = -SUM(P(j_loc,2)*ww(:,l))/ray*mode !* rj(l,m)
          fv(3) = (-mode*vitloc(2,l)+vitloc(3,l) +(mode*vitloc(3,l)-vitloc(2,l))*mode)/ray**2
          !--------calcul de la quatrieme composante : u0(4,:) <--> f(th,m,s)

          fs(4) = SUM(ff(j_loc,4) * ww(:,l)) !* rj(l,m)
          ft(4) = SUM(V1m(j_loc,4) * ww(:,l)) !* rj(l,m)
          fp(4) = SUM(P(j_loc,1)*ww(:,l))/ray *mode !* rj(l,m)
          fv(4) =  (mode*vitloc(1,l)+vitloc(4,l) +(mode*vitloc(4,l)+vitloc(1,l))*mode)/ray**2
          !--------calcul de la cinquieme composante : u0(5,:) <--> f(z,m,c)

          fs(5) = SUM(ff(j_loc,5) * ww(:,l)) !* rj(l,m)
          ft(5) = SUM(V1m(j_loc,5) * ww(:,l)) !* rj(l,m)
          fp(5) = -SUM(P(j_loc,1)*dw_loc(2,:)) !* rj(l,m)
          fv(5) =  vitloc(5,l)*(mode/ray)**2
          !--------calcul de la sixieme composante : u0(6,:) <--> f(z,m,s)

          fs(6) = SUM(ff(j_loc,6) * ww(:,l)) !* rj(l,m)
          ft(6) = SUM(V1m(j_loc,6) * ww(:,l)) !* rj(l,m)
          fp(6) = -SUM(P(j_loc,2)*dw_loc(2,:)) !* rj(l,m)
          fv(6) = vitloc(6,l)*(mode/ray)**2
          !-------calcul du second membre pour le terme nonlineaire------------------------
          !-------------------------------------------------------------------------------

          mult(1) = visc1; mult(4) = visc1; mult(5) = visc1
          mult(2) = visc2; mult(3) = visc2; mult(6) = visc2
          fv = mult*fv

          smb =  (ft+fp+fs+fv-rotv_v(index,:))*ray*rj(l,m)
          DO type = 1, 6
             grad(:,type,l) =  mult(type)*grad(:,type,l)*ray*rj(l,m)
          END DO
          !---------- calcul du second membre total
          DO j=1,6
             DO ni = 1, n_w 
                u0(jj(ni,m),j) = u0(jj(ni,m),j) +  ww(ni,l)*smb(j) + SUM(dw_loc(:,ni)*grad(:,j,l))
             ENDDO
          ENDDO

       ENDDO
    ENDDO

    WRITE(*,*) ' CFL = ', cfl, ' R_eff ', 1/R_eff
    !IF (cfl > 10.d0) THEN
    !   WRITE(*,*) ' CFL too large '
    !   STOP
    !END IF


  END SUBROUTINE qs_ns_stab_2008

  SUBROUTINE qs_stokes_2006(mesh,mode,ff,V1m,P,dt,u0)
    !=================================
    !Calcul du second membre pour le probleme de Stokes
    !sans terme non lineaire et source (ff)

    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:),               INTENT(IN)   :: ff  !analytique (type,noeud)
    REAL(KIND=8), DIMENSION(:,:),               INTENT(OUT)  :: u0   !u0(type,noeud)
    TYPE(mesh_type), TARGET                                  :: mesh  
    REAL(KIND=8), DIMENSION(:,:),               INTENT(IN)   :: V1m   !V(type,noeud) 
    REAL(KIND=8), DIMENSION(:,:),               INTENT(IN)   :: P
    REAL(KIND=8),                               INTENT(IN)   :: dt 
    INTEGER                     ,               INTENT(IN)   :: mode  

    INTEGER                                                  ::  m, l , i , ni , j
    REAL(KIND=8), DIMENSION(6)                               :: fs , ft , fp, fnl, smb
    INTEGER                                                  :: N
    INTEGER, DIMENSION(mesh%gauss%n_w)                       :: j_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w)   :: dw_loc     
    LOGICAL,                                      SAVE       :: once = .TRUE.
    LOGICAL,                                      SAVE       :: once_sft = .TRUE.
    REAL(KIND=8):: user_time
    !EXTERNAL :: user_time
    REAL(KIND=8) :: tps, dummy, tt                     

    !fs pour le terme de forcage
    !ft pour le terme temporelle   
    !fp pour le gradient de pression

    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8)   :: ray

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    tps = 0.d0
    tt = 0.d0

    u0   = 0.d0
    fs   = 0.d0
    ft   = 0.d0
    fp   = 0.d0
    fnl  = 0.d0

    DO m = 1, me
       j_loc = jj(:,m)
       DO l = 1, l_G
          dw_loc = dw(:,:,l,m)

          !--------On calcul le rayon du point gauss
          ray = SUM(mesh%rr(1,j_loc)*ww(:,l))

          !-calcul des seconds membres pr les termes de forçage, temporels et de gradP 
          !--------calcul de la premiere composante : u0(1,:) <--> f(r,m,c)

          fs(1) = SUM(ff(j_loc,1) * ww(:,l))
          ft(1) = 1/(2*dt) * SUM(V1m(j_loc,1) * ww(:,l)) 
          fp(1) = -SUM(P(j_loc,1)*dw_loc(1,:))

          !--------calcul de la seconde composante : u0(2,:) <--> f(r,m,s)

          fs(2) = SUM(ff(j_loc,2) * ww(:,l)) !* rj(l,m)
          ft(2) = 1/(2*dt) * SUM(V1m(j_loc,2) * ww(:,l)) !* rj(l,m)
          fp(2) = -SUM(P(j_loc,2)*dw_loc(1,:)) !* rj(l,m)

          !--------calcul de la troisieme composante : u0(3,:) <--> f(th,m,c)

          fs(3) = SUM(ff(j_loc,3) * ww(:,l)) !* rj(l,m)
          ft(3) = 1/(2*dt) * SUM(V1m(j_loc,3) * ww(:,l)) !* rj(l,m)
          fp(3) = -SUM(P(j_loc,2)*ww(:,l))/ray*mode !* rj(l,m)

          !--------calcul de la quatrieme composante : u0(4,:) <--> f(th,m,s)

          fs(4) = SUM(ff(j_loc,4) * ww(:,l)) !* rj(l,m)
          ft(4) = 1/(2*dt) * SUM(V1m(j_loc,4) * ww(:,l)) !* rj(l,m)
          fp(4) = SUM(P(j_loc,1)*ww(:,l))/ray *mode !* rj(l,m)

          !--------calcul de la cinquieme composante : u0(5,:) <--> f(z,m,c)

          fs(5) = SUM(ff(j_loc,5) * ww(:,l)) !* rj(l,m)
          ft(5) = 1/(2*dt) * SUM(V1m(j_loc,5) * ww(:,l)) !* rj(l,m)
          fp(5) = -SUM(P(j_loc,1)*dw_loc(2,:)) !* rj(l,m)

          !--------calcul de la sixieme composante : u0(6,:) <--> f(z,m,s)

          fs(6) = SUM(ff(j_loc,6) * ww(:,l)) !* rj(l,m)
          ft(6) = 1/(2*dt) * SUM(V1m(j_loc,6) * ww(:,l)) !* rj(l,m)
          fp(6) = -SUM(P(j_loc,2)*dw_loc(2,:)) !* rj(l,m)

          !-------calcul du second membre pour le terme nonlineaire------------------------
          !-------------------------------------------------------------------------------

          smb =  (ft+fp+fs)*ray*rj(l,m)      
          !write(95,*) m,l,ft(3)
          !---------- calcul du second membre total
          !smb = 1.d0
          DO j=1,6
             DO ni = 1, n_w 
                u0(jj(ni,m),j) = u0(jj(ni,m),j) +  ww(ni,l)*smb(j)
             ENDDO
          ENDDO

       ENDDO
    ENDDO


  END SUBROUTINE qs_stokes_2006



  SUBROUTINE qs_01_rot(mesh,mode,V,u0)
    !=================================
    !sans le terme de couplage


    USE Gauss_points

    TYPE(mesh_type), TARGET                                  :: mesh
    INTEGER,                                    INTENT(IN)   :: mode
    REAL(KIND=8), DIMENSION(mesh%np,6),         INTENT(IN)   :: V  
    REAL(KIND=8), DIMENSION(:,:),               INTENT(OUT)  :: u0

    INTEGER                                                  ::  m, l , i , ni , j
    REAL(KIND=8), DIMENSION(6)                               :: frot
    REAL(KIND=8):: user_time
    !EXTERNAL :: user_time
    REAL(KIND=8) :: tps, dummy, tt  
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8)   :: ray

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    u0   = 0.d0
    frot   = 0.d0      

    DO m = 1, me

       DO l = 1, l_G

          !--------On calcul le rayon du point gauss
          ray = 0
          DO ni = 1, n_w;  i = jj(ni,m)
             ray = ray + mesh%rr(1,i)*ww(ni,l)
          ENDDO
          !-----------------rotationnel ------------------------------
          !coeff sur les cosinus 
          !ancien format 
          !frot(1)  = mode/ray*SUM(V(6,jj(:,m))*ww(:,l)) &
          !            -SUM(V(3,jj(:,m))*dw(2,:,l,m))
          !frot(3)  = SUM(V(1,jj(:,m))*dw(2,:,l,m)) &
          !           -SUM(V(5,jj(:,m))*dw(1,:,l,m))
          !frot(5)  = 1/ray*SUM(V(3,jj(:,m))*ww(:,l))+SUM(V(3,jj(:,m)) &
          !           *dw(1,:,l,m))-mode/ray*SUM(V(2,jj(:,m))*ww(:,l))
          !                       !coeff sur les sinus 
          !IF (mode /= 0) THEN         
          !   frot(2) = -mode/ray*SUM(V(5,jj(:,m))*ww(:,l)) &
          !            -SUM(V(4,jj(:,m))*dw(2,:,l,m))
          !   frot(4) = SUM(V(2,jj(:,m))*dw(2,:,l,m)) &
          !             -SUM(V(6,jj(:,m))*dw(1,:,l,m))
          !   frot(6) = 1/ray*SUM(V(4,jj(:,m))*ww(:,l))+SUM(V(4,jj(:,m)) &
          !             *dw(1,:,l,m))+mode/ray*SUM(V(1,jj(:,m))*ww(:,l))
          !ENDIF


          frot(1)  = mode/ray*SUM(V(jj(:,m),6)*ww(:,l)) &
               -SUM(V(jj(:,m),3)*dw(2,:,l,m))
          frot(3)  = SUM(V(jj(:,m),1)*dw(2,:,l,m)) &
               -SUM(V(jj(:,m),5)*dw(1,:,l,m))
          frot(5)  = 1/ray*SUM(V(jj(:,m),3)*ww(:,l))+SUM(V(jj(:,m),3) &
               *dw(1,:,l,m))-mode/ray*SUM(V(jj(:,m),2)*ww(:,l))
          !coeff sur les sinus 
          IF (mode /= 0) THEN         
             frot(2) = -mode/ray*SUM(V(jj(:,m),5)*ww(:,l)) &
                  -SUM(V(jj(:,m),4)*dw(2,:,l,m))
             frot(4) = SUM(V(jj(:,m),2)*dw(2,:,l,m)) &
                  -SUM(V(jj(:,m),6)*dw(1,:,l,m))
             frot(6) = 1/ray*SUM(V(jj(:,m),4)*ww(:,l))+SUM(V(jj(:,m),4) &
                  *dw(1,:,l,m))+mode/ray*SUM(V(jj(:,m),1)*ww(:,l))
          ENDIF
          frot = frot * rj(l,m)

          !---------- calcul du second membre total

          DO j=1,6
             DO ni = 1, n_w 
                u0(jj(ni,m),j) = u0(jj(ni,m),j) +  ww(ni,l)  * ray * frot(j)  
             ENDDO
          ENDDO

       ENDDO
    ENDDO

  END SUBROUTINE qs_01_rot


  SUBROUTINE average(mesh,visc)
    USE Gauss_points
    USE sub_plot
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:),   INTENT(OUT)   :: visc
    TYPE(mesh_type)                             :: mesh
    REAL(KIND=8), DIMENSION(mesh%me)            :: vint
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE  :: surf, sloc
    LOGICAL, SAVE :: once=.TRUE.
    INTEGER :: m, n, mp, it

    IF (once) THEN
       once =.FALSE.
       ALLOCATE(surf(mesh%me),sloc(mesh%me))
       DO m = 1, mesh%me
          surf(m) = SUM(mesh%gauss%rj(:,m))
       END DO
       DO m = 1, mesh%me
          sloc(m) = 0.d0
          DO n = 1, 3
             IF (mesh%neigh(n,m)==0) CYCLE
             mp = mesh%neigh(n,m)
             sloc(m) = sloc(m) + surf(m)+surf(mp)
          END DO
       END DO

    END IF
    vint = 0.d0
    DO m = 1, mesh%me
       DO n = 1, 3
          IF (mesh%neigh(n,m)==0) CYCLE
          mp = mesh%neigh(n,m)
          vint(m) = vint(m)+ visc(m)*surf(m) + visc(mp)*surf(mp)
       END DO
       vint(m) = vint(m)/sloc(m)
    END DO
    visc = vint

  END SUBROUTINE average

  SUBROUTINE clement_c(mesh,phi,phi_s)
    !USE GAUSS_POINTS
    USE def_type_mesh
    IMPLICIT NONE
    TYPE(mesh_type)                             :: mesh
    REAL(KIND=8), DIMENSION(:)                  :: phi, phi_s

    INTEGER :: m, l, k, i
    REAL(KIND=8), DIMENSION(mesh%np)                :: a0, b1, b2
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE   :: s, a11, a22, a12, oneoverdet
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, SAVE :: rg
    INTEGER,      DIMENSION(mesh%gauss%n_w)         :: j_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w)         :: phi_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d)         :: r_loc
    REAL(KIND=8) :: x1l, x2l, phil, a1, a2, det
    LOGICAL :: const
    LOGICAL, SAVE :: once=.TRUE.

    IF (mesh%gauss%n_w/=3) THEN
       WRITE(*,*) ' BUG CLEMENT: Only P1 programmed yet'
       phi_s= phi
       RETURN
    END IF
    CONST=.FALSE.
    IF (SIZE(phi)==mesh%me) const=.TRUE.

    IF (once) THEN
       once = .FALSE.
       ALLOCATE(rg(mesh%gauss%k_d,mesh%np), s(mesh%np))
       s = 0.d0
       rg = 0.d0
       DO m = 1, mesh%me
          j_loc = mesh%jj(:,m)
          s(j_loc(:)) = s(j_loc(:)) + SUM(mesh%gauss%rj(:,m))
          DO k = 1, mesh%gauss%k_d
             r_loc(k) = 0.d0
             DO l = 1, mesh%gauss%l_G
                r_loc(k) =  r_loc(k) + SUM(mesh%rr(k,j_loc(:))*mesh%gauss%ww(:,l))*mesh%gauss%rj(l,m)
             END DO
             rg(k,j_loc(:)) = rg(k,j_loc(:)) + r_loc(k) 
          END DO
       END DO
       DO i = 1, mesh%np
          rg(:,i) = rg(:,i)/s(i)
       END DO
       ALLOCATE(a11(mesh%np), a22(mesh%np), a12(mesh%np), oneoverdet(mesh%np))
       a11 = 0.d0
       a12 = 0.d0
       a22 = 0.d0
       DO m = 1, mesh%me
          j_loc = mesh%jj(:,m)
          DO l = 1, mesh%gauss%l_G
             x1l = SUM(mesh%rr(1,j_loc(:))*mesh%gauss%ww(:,l))
             x2l = SUM(mesh%rr(2,j_loc(:))*mesh%gauss%ww(:,l))
             a11(j_loc(:)) = a11(j_loc(:)) + (x1l-rg(1,j_loc(:)))**2*mesh%gauss%rj(l,m)
             a22(j_loc(:)) = a22(j_loc(:)) + (x2l-rg(2,j_loc(:)))**2*mesh%gauss%rj(l,m)
             a12(j_loc(:)) = a12(j_loc(:)) + (x1l-rg(1,j_loc(:)))*(x2l-rg(2,j_loc(:)))*mesh%gauss%rj(l,m)
          END DO
       END DO
       oneoverdet = 1.d0/(a11*a22-a12*a12)
    END IF


    a0 = 0.d0
    b1 = 0.d0
    b2 = 0.d0
    DO m = 1, mesh%me
       j_loc = mesh%jj(:,m)
       !IF (const) THEN
          phi_loc = phi(m)
       !ELSE
       !   phi_loc = phi(j_loc)
       !END IF
       DO l = 1, mesh%gauss%l_G
          x1l = SUM(mesh%rr(1,j_loc(:))*mesh%gauss%ww(:,l))
          x2l = SUM(mesh%rr(2,j_loc(:))*mesh%gauss%ww(:,l))
          phil = SUM(phi_loc*mesh%gauss%ww(:,l))
          b1(j_loc(:)) =  b1(j_loc(:)) + (x1l-rg(1,j_loc(:)))*phil*mesh%gauss%rj(l,m)
          b2(j_loc(:)) =  b2(j_loc(:)) + (x2l-rg(2,j_loc(:)))*phil*mesh%gauss%rj(l,m)
          a0(j_loc) = a0(j_loc) + phil*mesh%gauss%rj(l,m)
       END DO
    END DO

    DO i = 1, mesh%np
       a0(i)=a0(i)/s(i)
       !det = a11(i)*a22(i)-a12(i)*a12(i)
       a1 =(b1(i)*a22(i)-b2(i)*a12(i))*oneoverdet(i)
       a2 =(b2(i)*a11(i)-b1(i)*a12(i))*oneoverdet(i)
       phi_s(i) = a0(i)+a1*(mesh%rr(1,i)-rg(1,i))+a2*(mesh%rr(2,i)-rg(2,i))
    END DO

  END SUBROUTINE clement_c


  SUBROUTINE inject_clement(jj_c, jj_f, pp_c, pp_f)

    IMPLICIT NONE

    INTEGER,      DIMENSION(:,:), INTENT(IN)  :: jj_c, jj_f
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: pp_c
    REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: pp_f
    REAL(KIND=8) :: half = 0.5
    INTEGER:: m

    DO m = 1, SIZE(jj_f,2)
       pp_f(jj_f(1:3,m)) =  pp_c(jj_c(:,m))
       pp_f(jj_f(4,m)) = (pp_c(jj_c(2,m)) + pp_c(jj_c(3,m)))*half
       pp_f(jj_f(5,m)) = (pp_c(jj_c(3,m)) + pp_c(jj_c(1,m)))*half
       pp_f(jj_f(6,m)) = (pp_c(jj_c(1,m)) + pp_c(jj_c(2,m)))*half 
    END DO

  END SUBROUTINE inject_clement
END MODULE fem_s_axi
