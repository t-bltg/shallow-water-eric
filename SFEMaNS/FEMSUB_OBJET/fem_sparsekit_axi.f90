!
!
!Authors: Jean-Luc Guermond, Raphael Laguerre, Luigi Quartapelle, Copyrights 1996, 2000, 2004
!
MODULE fem_s_axi_M

CONTAINS



  SUBROUTINE qs_11_M (mesh, alpha, ia, ja,  a0)      ! laplacien en (r,z)
    !=================================================

    !  alpha << (Dw), (D_) >>   ===>   A0    incremental accumulation

    USE Gauss_points

    IMPLICIT NONE

    REAL(KIND=8),                 INTENT(IN)    :: alpha
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia 
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ja 
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0

    INTEGER :: m, l, ni, nj, i, j, p
    REAL(KIND=8) :: al, x

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8)    :: ray
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w) :: aij
    REAL(KIND=8), DIMENSION(:,:) , POINTER   :: rr


    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me
    rr => mesh%rr

    DO m = 1, me

       aij = 0

       DO l = 1, l_G
          !--------On calcule le rayon du point gauss
          ray = 0
          DO ni = 1, n_w;  i = jj(ni,m)
             ray = ray + rr(1,i)*ww(ni,l)
          END DO

          al = alpha * rj(l,m) * ray

          DO nj = 1, n_w;  j = jj(nj, m)
             DO ni = 1, n_w;  i = jj(ni, m)

                aij(ni,nj) = aij(ni,nj) +  al * SUM(dw(:,nj,l,m) &
                     * dw(:,ni,l,m))

             ENDDO
          ENDDO

       ENDDO

       DO ni = 1, n_w; i = jj(ni, m)
          DO nj = 1, n_w;  j = jj(nj, m)
             DO p = ia(i),  ia(i+1) - 1
                IF (ja(p) == j) THEN;  a0(p) = a0(p) + aij(ni,nj);  EXIT
                ENDIF
             ENDDO
          END DO
       END DO


    ENDDO

  END SUBROUTINE qs_11_M

  SUBROUTINE qs_00_M(mesh, alpha, ia, ja, a0)      ! masse
    !=================================================

    !  alpha < r w, _ >   ===>   A0    incremental accumulation

    USE Gauss_points

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:,:)  , POINTER     :: rr
    REAL(KIND=8),                 INTENT(IN)    :: alpha
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0

    INTEGER :: mm, m, l, ni, nj, i, j, p
    REAL(KIND=8) :: al, x, ray

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w) :: aij

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me
    rr => mesh%rr

    DO m = 1, me

       aij = 0

       DO l = 1, l_G

          !--------On calcule le rayon du point gauss
          ray = 0
          DO ni = 1, n_w;  i = jj(ni,m)
             ray = ray + rr(1,i)*ww(ni,l)
          END DO


          al = alpha * rj(l,m) * ray

          DO nj = 1, n_w;  j = jj(nj, m)
             DO ni = 1, n_w;  i = jj(ni, m)

                aij(ni,nj) = aij(ni,nj) +  ww(nj,l) * al * ww(ni,l)

             ENDDO
          ENDDO


       ENDDO

       DO ni = 1, n_w; i = jj(ni, m)
          DO nj = 1, n_w;  j = jj(nj, m)
             DO p = ia(i),  ia(i+1) - 1
                IF (ja(p) == j) THEN;  a0(p) = a0(p) + aij(ni,nj);  EXIT
                ENDIF
             ENDDO
          END DO
       END DO

    ENDDO

  END SUBROUTINE qs_00_M


  SUBROUTINE qs_mass_6comp(mesh, mass, ia, ja, a0)
    !========================================================================
    !mass pour les 6 composantes

    USE gauss_points

    IMPLICIT NONE

    REAL(KIND=8),                 INTENT(IN)    :: mass
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8), DIMENSION(:,:) , POINTER      :: rr

    INTEGER :: k, l, m, ni, nj, nk, i, j, p, np, ib, jb, kb, ki, kj, kk
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d) :: gl
    REAL(KIND=8) :: xij, masslm, viscolm, viscorlm
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w,mesh%gauss%l_G) :: wwprod
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w) :: aij   
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w) :: y
    REAL(KIND=8)    :: ray, eps1, eps2

    CALL gauss(mesh)

    jj => mesh%jj
    me => mesh%me
    rr => mesh%rr
    np = SIZE(mesh%rr,2)


    DO l = 1, l_G 
       DO ni = 1, n_w 
          DO nj = 1, n_w 
             wwprod(ni,nj,l) = ww(ni,l)*ww(nj,l)
          END DO
       END DO
    END DO

    DO m = 1, me

       aij = 0.d0

       DO l = 1, l_G


          !--------On calcule le rayon du point gauss
          ray = 0
          DO ni = 1, n_w;  i = jj(ni,m)
             ray = ray + rr(1,i)*ww(ni,l)
          END DO

          DO nj = 1, n_w; j = jj(nj, m)
             DO ni = 1, n_w;  i = jj(ni, m)

                !blocs diagonaux
                aij(ni,nj) =  aij(ni,nj) + mass*ray*wwprod(ni,nj,l)*rj(l,m)

             ENDDO
          ENDDO


       ENDDO

       DO ki= 1, 6  
          DO ni = 1, n_w 
             i = jj(ni, m)
             ib = i + (ki-1)*np 
             DO kj = 1, 6
                DO nj = 1, n_w
                   j = jj(nj, m)
                   jb = j + (kj-1)*np 
                   DO p = ia(ib),  ia(ib+1) - 1
                      IF (ja(p) == jb) THEN

                         IF (ki==kj) THEN
                            a0(p) = a0(p) + aij(ni,nj)
                         END IF

                         EXIT
                      ENDIF
                   END DO
                END DO
             END DO
          END DO
       END DO

    ENDDO



  END SUBROUTINE qs_mass_6comp


  SUBROUTINE qs_diff_mass_adv_M(mesh, visco, alpha, gg, ia, ja, a0)   !advection diffusion
    !==========================================================

    !  << visco (Dw), D_) >>
    !  +  alpha * < w, _ >
    !  +  < w, (gg.D)_ >
    !  ===>  a0

    USE gauss_points

    IMPLICIT NONE

    REAL(KIND=8),                 INTENT(IN)    :: visco, alpha
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)    :: gg
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8), DIMENSION(:,:) , POINTER      :: rr

    INTEGER :: k, l, m, ni, nj, i, j, p
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d) :: gl
    REAL(KIND=8) :: xij, masslm, viscolm
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w,mesh%gauss%l_G) :: wwprod
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w) :: aij
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w) :: y
    REAL(KIND=8)    :: ray

    CALL gauss(mesh)

    jj => mesh%jj
    me => mesh%me
    rr => mesh%rr

    DO l = 1, l_G 
       DO ni = 1, n_w 
          DO nj = 1, n_w 
             wwprod(ni,nj,l) = ww(ni,l)*ww(nj,l)
          END DO
       END DO
    END DO

    DO m = 1, me

       aij = 0
       DO l = 1, l_G


          !--------On calcule le rayon du point gauss
          ray = 0
          DO ni = 1, n_w;  i = jj(ni,m)
             ray = ray + rr(1,i)*ww(ni,l)
          END DO


          gl = 0.
          DO k = 1, k_d
             DO ni =1 ,n_w
                gl(k) = gl(k) + gg(k, jj(ni,m)) * ww(ni,l)
             END DO
          ENDDO
          viscolm = visco*rj(l,m)
          masslm  = alpha*rj(l,m)

          y = 0.
          DO ni = 1, n_w; 
             DO k = 1, k_d
                y(ni) =  y(ni) + gl(k) * dw(k,ni,l,m)
             END DO
          END DO
          y = rj(l,m)*y

          DO ni = 1, n_w; i = jj(ni, m)
             DO nj = 1, n_w;  j = jj(nj, m)

                xij = 0.
                DO k = 1, k_d
                   xij =  xij + dw(k,nj,l,m) * dw(k,ni,l,m)
                END DO

                aij(ni,nj) =  aij(ni,nj) + ray * (viscolm*xij + masslm*wwprod(ni,nj,l) + &
                     y(nj)*ww(ni,l))

             ENDDO
          ENDDO

       ENDDO

       DO ni = 1, n_w; i = jj(ni, m)
          DO nj = 1, n_w;  j = jj(nj, m)
             DO p = ia(i),  ia(i+1) - 1
                IF (ja(p) == j) THEN;  a0(p) = a0(p) + aij(ni,nj);  EXIT
                ENDIF
             ENDDO
          END DO
       END DO

    ENDDO

  END SUBROUTINE  qs_diff_mass_adv_M




  SUBROUTINE qs_11_M_theta(mesh, alpha, ia, ja, a0) ! masse pour le cas spectral
    !=================================================
    ! 
    !  alpha < r w, _ > / r^2   ===>   A0    incremental accumulation

    USE Gauss_points

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:,:)  , POINTER     :: rr
    REAL(KIND=8),                 INTENT(IN)    :: alpha
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0

    INTEGER :: mm, m, l, ni, nj, i, j, p
    REAL(KIND=8) :: al, x, ray

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w) :: aij

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me
    rr => mesh%rr

    DO m = 1, me

       aij = 0

       DO l = 1, l_G

          !--------On calcule le rayon du point gauss
          ray = 0
          DO ni = 1, n_w;  i = jj(ni,m)
             ray = ray + rr(1,i)*ww(ni,l)
          END DO


          al = alpha * rj(l,m)

          DO nj = 1, n_w;  j = jj(nj, m)
             DO ni = 1, n_w;  i = jj(ni, m)

                aij(ni,nj) = aij(ni,nj) +  ww(nj,l) * al * ww(ni,l) / ray

             ENDDO
          ENDDO


       ENDDO

       DO ni = 1, n_w; i = jj(ni, m)
          DO nj = 1, n_w;  j = jj(nj, m)
             DO p = ia(i),  ia(i+1) - 1
                IF (ja(p) == j) THEN;  a0(p) = a0(p) + aij(ni,nj);  EXIT
                ENDIF
             ENDDO
          END DO
       END DO

    ENDDO

  END SUBROUTINE qs_11_M_theta


  SUBROUTINE qs_inst_mass_diff(mesh, visco,beta, ia, ja, a0)
    !=========================================================
    !On calcule l'operateur (beta - visco*lap)(.)

    USE gauss_points

    IMPLICIT NONE

    REAL(KIND=8),                 INTENT(IN)    :: visco, beta
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8), DIMENSION(:,:) , POINTER      :: rr

    INTEGER :: k, l, m, ni, nj, i, j, p
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d) :: gl
    REAL(KIND=8) :: xij, masslm, viscolm
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w,mesh%gauss%l_G) :: wwprod
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w) :: aij
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w) :: y
    REAL(KIND=8)    :: ray

    CALL gauss(mesh)

    jj => mesh%jj
    me => mesh%me
    rr => mesh%rr

    DO l = 1, l_G 
       DO ni = 1, n_w 
          DO nj = 1, n_w 
             wwprod(ni,nj,l) = ww(ni,l)*ww(nj,l)
          END DO
       END DO
    END DO

    DO m = 1, me

       aij = 0
       DO l = 1, l_G


          !--------On calcule le rayon du point gauss
          ray = 0
          DO ni = 1, n_w;  i = jj(ni,m)
             ray = ray + rr(1,i)*ww(ni,l)
          END DO

          viscolm = visco*rj(l,m)
          masslm  = beta*rj(l,m)

          DO nj = 1, n_w; j = jj(nj, m)
             DO ni = 1, n_w;  i = jj(ni, m)

                xij = 0.
                DO k = 1, k_d
                   xij =  xij + dw(k,nj,l,m) * dw(k,ni,l,m)
                END DO

                aij(ni,nj) =  aij(ni,nj) + ray * (viscolm*xij + masslm*wwprod(ni,nj,l))

             ENDDO
          ENDDO

       ENDDO

       DO ni = 1, n_w; i = jj(ni, m)
          DO nj = 1, n_w;  j = jj(nj, m)
             DO p = ia(i),  ia(i+1) - 1
                IF (ja(p) == j) THEN;  a0(p) = a0(p) + aij(ni,nj);  EXIT
                ENDIF
             ENDDO
          END DO
       END DO

    ENDDO

  END SUBROUTINE qs_inst_mass_diff



  SUBROUTINE qs_total_3d(aa0,aa1,aam,mode)

    !regroupement de la matrice total pour le mode m dans *m , a partir des
    !matrices constantes *0 et *1(=qs_00_M_theta)
    !aa0 + m^2 aa1 ==> aam

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:),   INTENT(IN)    :: aa0
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)    :: aa1
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: aam
    INTEGER ::   j
    INTEGER                   ,   INTENT(IN)    :: mode

    DO j=1,SIZE(aa0(:))
       aam(j) = aa0(j) + (mode**2)*aa1(j)
    ENDDO

  END SUBROUTINE qs_total_3d

  SUBROUTINE qs_adv_diff_axi_3d(mesh, visco,beta,gg, ia, ja, a0)



    USE gauss_points

    IMPLICIT NONE

    REAL(KIND=8),                 INTENT(IN)    :: visco, beta
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8), DIMENSION(:,:) , POINTER      :: rr
    REAL(KIND=8), DIMENSION(:,:) , INTENT(IN)   :: gg

    INTEGER :: k, l, m, ni, nj, i, j, p
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d) :: gl
    REAL(KIND=8) :: xij, masslm, viscolm
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w,mesh%gauss%l_G) :: wwprod
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w) :: aij
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w) :: y
    REAL(KIND=8)    :: ray

    CALL gauss(mesh)

    jj => mesh%jj
    me => mesh%me
    rr => mesh%rr

    DO l = 1, l_G 
       DO ni = 1, n_w 
          DO nj = 1, n_w 
             wwprod(ni,nj,l) = ww(ni,l)*ww(nj,l)
          END DO
       END DO
    END DO

    DO m = 1, me

       aij = 0
       DO l = 1, l_G


          !--------On calcule le rayon du point gauss
          ray = 0
          DO ni = 1, n_w;  i = jj(ni,m)
             ray = ray + rr(1,i)*ww(ni,l)
          END DO

          gl = 0.
          DO k = 1, k_d
             DO ni =1 ,n_w
                gl(k) = gl(k) + gg(k, jj(ni,m)) * ww(ni,l)
             END DO
          ENDDO

          y = 0.
          DO ni = 1, n_w; 
             DO k = 1, k_d
                y(ni) =  y(ni) + gl(k) * dw(k,ni,l,m)
             END DO
          END DO
          y = rj(l,m)*y

          viscolm = visco*rj(l,m)
          masslm  = beta*rj(l,m)

          DO ni = 1, n_w; i = jj(ni, m)
             DO nj = 1, n_w;  j = jj(nj, m)

                xij = 0.
                DO k = 1, k_d
                   xij =  xij + dw(k,nj,l,m) * dw(k,ni,l,m)
                END DO

                aij(ni,nj) =  aij(ni,nj) + ray * (viscolm*xij + masslm*wwprod(ni,nj,l) &
                     + y(nj)*ww(ni,l))

             ENDDO
          ENDDO

       ENDDO

       DO ni = 1, n_w; i = jj(ni, m)
          DO nj = 1, n_w;  j = jj(nj, m)
             DO p = ia(i),  ia(i+1) - 1
                IF (ja(p) == j) THEN;  a0(p) = a0(p) + aij(ni,nj);  EXIT
                ENDIF
             ENDDO
          END DO
       END DO

    ENDDO

  END SUBROUTINE qs_adv_diff_axi_3d


  !------------------------------------------------------------------------
  !subroutine pour le cas vectoriel
  !----------------------------------------------------------------------

  SUBROUTINE qs_inst_mass_diff_vect(mesh, visco,beta, ia, ja, a0)
    !=========================================================
    ! on calcule l'operateur (beta - visco*(lap-1/r^2))(.)
    !lap correspond ici a la premiere partie du laplacien, qui ne 
    !varie pas avec le mode. L'autre partie est calculee avec la procedure
    !qs_11_mtheta

    USE gauss_points

    IMPLICIT NONE

    REAL(KIND=8),                 INTENT(IN)    :: visco, beta
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8), DIMENSION(:,:) , POINTER      :: rr

    INTEGER :: k, l, m, ni, nj, i, j, p
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d) :: gl
    REAL(KIND=8) :: xij, masslm, viscolm
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w,mesh%gauss%l_G) :: wwprod
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w) :: aij
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w) :: y
    REAL(KIND=8)    :: ray

    CALL gauss(mesh)

    jj => mesh%jj
    me => mesh%me
    rr => mesh%rr

    DO l = 1, l_G 
       DO ni = 1, n_w 
          DO nj = 1, n_w 
             wwprod(ni,nj,l) = ww(ni,l)*ww(nj,l)
          END DO
       END DO
    END DO

    DO m = 1, me

       aij = 0
       DO l = 1, l_G


          !--------On calcule le rayon du point gauss
          ray = 0
          DO ni = 1, n_w;  i = jj(ni,m)
             ray = ray + rr(1,i)*ww(ni,l)
          END DO

          viscolm = visco*rj(l,m)
          masslm  = beta*rj(l,m)

          DO nj = 1, n_w; j = jj(nj, m)
             DO ni = 1, n_w;  i = jj(ni, m)

                xij = 0.
                DO k = 1, k_d
                   xij =  xij + dw(k,nj,l,m) * dw(k,ni,l,m)
                END DO

                aij(ni,nj) =  aij(ni,nj) + ray * (viscolm*(xij + wwprod(ni,nj,l)/(ray**2)) + masslm*wwprod(ni,nj,l))

             ENDDO
          ENDDO

       ENDDO

       DO ni = 1, n_w; i = jj(ni, m)
          DO nj = 1, n_w;  j = jj(nj, m)
             DO p = ia(i),  ia(i+1) - 1
                IF (ja(p) == j) THEN;  a0(p) = a0(p) + aij(ni,nj);  EXIT
                ENDIF
             ENDDO
          END DO
       END DO

    ENDDO

  END SUBROUTINE qs_inst_mass_diff_vect

  SUBROUTINE qs_mass_diff_vect(type_op, mesh, visco, mass, stab, mode, ia, ja, a0)
    !========================================================================
    ! laplacien vectoriel qui renvoie une matrice 2np * 2np, pour deux composantes
    ! (V1,V4), (V2,V3) et une structure scaleire np*np pour (V5,V6)
    ! calcul de l'operateur mass-visco(lap-eps1*1/r^2) et eps2*2m*visco/r^2
    ! eps1 et eps2 sont des parametres qui definissent le type d'operateur
    ! type de l'operateur : 1 pour les composantes 1 et 4
    !                       2 pour les composantes 2 et 3
    !                       3 pour les composantes 5 et 6
    !------------------------------------------------------------------------

    USE gauss_points

    IMPLICIT NONE


    INTEGER     ,                 INTENT(IN)    :: type_op, mode  
    REAL(KIND=8),                 INTENT(IN)    :: visco, mass, stab
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8), DIMENSION(:,:) , POINTER      :: rr

    INTEGER :: k, l, m, ni, nj, i, j, p, np, ib, jb, ki, kj, k_max
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d) :: gl
    REAL(KIND=8) :: xij, masslm, viscolm, viscorlm, hloc
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w,mesh%gauss%l_G) :: wwprod
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w) :: aij, bij   
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w) :: y
    REAL(KIND=8)    :: ray, eps1, eps2

    CALL gauss(mesh)

    jj => mesh%jj
    me => mesh%me
    rr => mesh%rr
    np = SIZE(mesh%rr,2)

    IF (type_op == 1) THEN
       eps1 = 1.d0
       eps2 = 1.d0
       k_max = 2 !Structure vectorielle
    ELSEIF (type_op == 2) THEN
       eps1 = 1.d0
       eps2 = -1.d0
       k_max = 2 !Structure vectorielle
    ELSEIF (type_op == 3) THEN  
       !cas du laplacien scalaire   
       eps1 = 0.d0
       eps2 = 0.d0
       k_max = 1 !Structure scalaire
    ELSE
       WRITE(*,*) 'probleme de type d''operateur'
       STOP
    ENDIF


    DO l = 1, l_G 
       DO ni = 1, n_w 
          DO nj = 1, n_w 
             wwprod(ni,nj,l) = ww(ni,l)*ww(nj,l)
          END DO
       END DO
    END DO

    DO m = 1, me
       hloc = SQRT(SUM(rj(:,m)))
       aij = 0.d0
       bij = 0.d0

       DO l = 1, l_G


          !--------On calcule le rayon du point gauss
          ray = 0
          DO ni = 1, n_w;  i = jj(ni,m)
             ray = ray + rr(1,i)*ww(ni,l)
          END DO

          viscolm  = (visco + stab*hloc)*rj(l,m) ! We add the constant stabilization here

          DO nj = 1, n_w; j = jj(nj, m)

             DO ni = 1, n_w;  i = jj(ni, m)

                !grad(u).grad(v) en r et z
                xij = 0.
                DO k = 1, k_d
                   xij =  xij + dw(k,nj,l,m) * dw(k,ni,l,m)
                END DO

                !blocs diagonaux
                aij(ni,nj) =  aij(ni,nj) + ray * viscolm* xij    &
                     + viscolm*eps1*wwprod(ni,nj,l)/ray &
                     + mass*ray*wwprod(ni,nj,l)*rj(l,m) &
                     + viscolm*mode**2*wwprod(ni,nj,l)/ray
                !blocs couplant              
                bij(ni,nj) = bij(ni,nj) + eps2*viscolm*2*mode*wwprod(ni,nj,l)/ray

             ENDDO
          ENDDO


       ENDDO

       DO ki= 1, k_max
          DO ni = 1, n_w 
             i = jj(ni, m)
             ib = i + (ki-1)*np 
             DO kj = 1, k_max
                DO nj = 1, n_w
                   j = jj(nj, m)
                   jb = j + (kj-1)*np 
                   DO p = ia(ib),  ia(ib+1) - 1
                      IF (ja(p) == jb) THEN
                         IF (ki==kj) THEN
                            a0(p) = a0(p) + aij(ni,nj)
                         ELSE
                            a0(p) = a0(p) + bij(ni,nj)
                         END IF
                         EXIT
                      ENDIF
                   END DO
                END DO
             END DO
          END DO
       END DO

    ENDDO



  END SUBROUTINE qs_mass_diff_vect

  SUBROUTINE qs_mass_diff_vect_b(type_op, mesh, visco, mass, mode, ia, ja, a0)
    !========================================================================
    ! laplacien vectoriel qui renvoie une matrice 2np * 2np, pour deux composantes
    ! (V1,V4), (V2,V3) ou (V5,V6)
    ! calcul de l'operateur mass-visco(lap-eps1*1/r^2) et eps2*2m*visco/r^2
    ! eps1 et eps2 sont des parametres qui definissent le type d'operateur
    ! type de l'operateur : 1 pour les composantes 1 et 4
    !                       2 pour les composantes 2 et 3
    !                       3 pour les composantes 5 et 6
    !------------------------------------------------------------------------

    USE gauss_points

    IMPLICIT NONE


    INTEGER     ,                 INTENT(IN)    :: type_op, mode  
    REAL(KIND=8),                 INTENT(IN)    :: visco, mass
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8), DIMENSION(:,:) , POINTER      :: rr

    INTEGER :: k, l, m, ni, nj, i, j, p, np, ib, jb, ki, kj
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d) :: gl
    REAL(KIND=8) :: xij, masslm, viscolm, viscorlm
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w,mesh%gauss%l_G) :: wwprod
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w) :: aij, bij   
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w) :: y
    REAL(KIND=8)    :: ray, eps1, eps2

    CALL gauss(mesh)

    jj => mesh%jj
    me => mesh%me
    rr => mesh%rr
    np = SIZE(mesh%rr,2)

    IF (type_op == 1) THEN
       eps1 = 1.d0
       eps2 = 1.d0
    ELSEIF (type_op == 2) THEN
       eps1 = 1.d0
       eps2 = -1.d0
    ELSEIF (type_op == 3) THEN  
       !cas du laplacien scalaire   
       eps1 = 0.d0
       eps2 = 0.d0
    ELSE
       WRITE(*,*) 'probleme de type d''operateur'
       STOP
    ENDIF


    DO l = 1, l_G 
       DO ni = 1, n_w 
          DO nj = 1, n_w 
             wwprod(ni,nj,l) = ww(ni,l)*ww(nj,l)
          END DO
       END DO
    END DO

    DO m = 1, me

       aij = 0.d0
       bij = 0.d0

       DO l = 1, l_G


          !--------On calcule le rayon du point gauss
          ray = 0
          DO ni = 1, n_w;  i = jj(ni,m)
             ray = ray + rr(1,i)*ww(ni,l)
          END DO

          viscolm  = visco*rj(l,m)
          !   masslm   = (mass + visco*mode**2/ray**2)*rj(l,m)
          !   viscorlm = 2*viscolm*mode/ray**2
          !
          !        DO nj = 1, n_w; j = jj(nj, m)
          !       
          !           DO ni = 1, n_w;  i = jj(ni, m)

          !               xij = 0.
          !!               DO k = 1, k_d
          !                xij =  xij + dw(k,nj,l,m) * dw(k,ni,l,m)
          !!             END DO
          !             
          !             aij(ni,nj) =  aij(ni,nj) + ray * (viscolm*(xij &
          !                 + eps1*wwprod(ni,nj,l)/(ray**2)) + masslm*wwprod(ni,nj,l))
          !             bij(ni,nj) = bij(ni,nj) + eps2*ray*viscorlm*wwprod(ni,nj,l)
          !          ENDDO
          !       ENDDO

          DO nj = 1, n_w; j = jj(nj, m)

             DO ni = 1, n_w;  i = jj(ni, m)

                !grad(u).grad(v) en r et z
                xij = 0.
                DO k = 1, k_d
                   xij =  xij + dw(k,nj,l,m) * dw(k,ni,l,m)
                END DO

                !blocs diagonaux
                aij(ni,nj) =  aij(ni,nj) + ray * viscolm* xij    &
                     + viscolm*eps1*wwprod(ni,nj,l)/ray &
                     + mass*ray*wwprod(ni,nj,l)*rj(l,m) &
                     + viscolm*mode**2*wwprod(ni,nj,l)/ray
                !blocs couplant              
                bij(ni,nj) = bij(ni,nj) + eps2*viscolm*2*mode*wwprod(ni,nj,l)/ray

             ENDDO
          ENDDO


       ENDDO

       DO ki= 1, 2  
          DO ni = 1, n_w 
             i = jj(ni, m)
             ib = i + (ki-1)*np 
             DO kj = 1, 2
                DO nj = 1, n_w
                   j = jj(nj, m)
                   jb = j + (kj-1)*np 
                   DO p = ia(ib),  ia(ib+1) - 1
                      IF (ja(p) == jb) THEN
                         IF (ki==kj) THEN
                            a0(p) = a0(p) + aij(ni,nj)
                         ELSE
                            a0(p) = a0(p) + bij(ni,nj)
                         END IF
                         EXIT
                      ENDIF
                   END DO
                END DO
             END DO
          END DO
       END DO

    ENDDO



  END SUBROUTINE qs_mass_diff_vect_b

  SUBROUTINE qs_inst_adv_diff_vect(mesh, visco,beta,gg, ia, ja, a0)
    !=========================================================
    !On calcule l'operateur (beta - visco*(lap-1/r^2)+Vrdr+Vzdz)(.)
    !lap correspond ici a la premiere partie du laplacien, qui ne 
    !varie pas avec le mode. L'autre partie est calculee avec la procedure
    !qs_11_mtheta

    USE gauss_points

    IMPLICIT NONE

    REAL(KIND=8),                 INTENT(IN)    :: visco, beta
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0
    REAL(KIND=8), DIMENSION(:,:) , INTENT(IN)   :: gg
    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8), DIMENSION(:,:) , POINTER      :: rr

    INTEGER :: k, l, m, ni, nj, i, j, p
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d) :: gl
    REAL(KIND=8) :: xij, masslm, viscolm
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w,mesh%gauss%l_G) :: wwprod
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w) :: aij
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w) :: y
    REAL(KIND=8)    :: ray

    CALL gauss(mesh)

    jj => mesh%jj
    me => mesh%me
    rr => mesh%rr

    DO l = 1, l_G 
       DO ni = 1, n_w 
          DO nj = 1, n_w 
             wwprod(ni,nj,l) = ww(ni,l)*ww(nj,l)
          END DO
       END DO
    END DO

    DO m = 1, me

       aij = 0
       DO l = 1, l_G


          !--------On calcule le rayon du point gauss
          ray = 0
          DO ni = 1, n_w;  i = jj(ni,m)
             ray = ray + rr(1,i)*ww(ni,l)
          END DO

          viscolm = visco*rj(l,m)
          masslm  = beta*rj(l,m)
          gl = 0.
          DO k = 1, k_d
             DO ni =1 ,n_w
                gl(k) = gl(k) + gg(k, jj(ni,m)) * ww(ni,l)
             END DO
          ENDDO

          y = 0.
          DO ni = 1, n_w; 
             DO k = 1, k_d
                y(ni) =  y(ni) + gl(k) * dw(k,ni,l,m)
             END DO
          END DO
          y = rj(l,m)*y

          DO nj = 1, n_w; j = jj(nj, m)
             DO ni = 1, n_w;  i = jj(ni, m)

                xij = 0.
                DO k = 1, k_d
                   xij =  xij + dw(k,nj,l,m) * dw(k,ni,l,m)
                END DO

                aij(ni,nj) =  aij(ni,nj) + ray * (viscolm*(xij + wwprod(ni,nj,l)/(ray**2)) &
                     + masslm*wwprod(ni,nj,l) + y(nj) * ww(ni,l))

             ENDDO
          ENDDO

       ENDDO

       DO ni = 1, n_w; i = jj(ni, m)
          DO nj = 1, n_w;  j = jj(nj, m)
             DO p = ia(i),  ia(i+1) - 1
                IF (ja(p) == j) THEN;  a0(p) = a0(p) + aij(ni,nj);  EXIT
                ENDIF
             ENDDO
          END DO
       END DO

    ENDDO

  END SUBROUTINE qs_inst_adv_diff_vect



END MODULE fem_s_axi_M














