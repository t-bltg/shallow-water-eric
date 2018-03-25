!  =================================
!  Entry points for SCALAR equations
!  =================================

!  ff  and  u0  are Scalar functions  (fs is defined only on the boundary)
!  gg  and  v0  are Vector fields     (gs is defined only on the boundary)

!  w   denotes either a scalar of a vector weigthing function
!  ws  denotes either a scalar of a vector weigthing function
!      defined only on the boundary

!  < _ , _ >     means integration over elements in the list  m0
! << _ , _ >>    means scalar product and integration over elements in  m0
!  < _ , _ >_s   means surface integration over boundary elements in  ms0
! << _ , _ >>_s  means scalar product and surface integration over b.e. ms0

!     D  is Nabla spatial derivative operator (grad, div, curl)

MODULE  fem_s

CONTAINS

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

SUBROUTINE qs_00 (mesh, ff,  u0)
!=================================

!  < w, f >   ===>   u0

   USE Gauss_points

   IMPLICIT NONE
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0

   INTEGER ::  m, l
   REAL(KIND=8) :: fl

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   u0 = 0

   DO m = 1, me
      DO l = 1, l_G

         fl = SUM(ff(jj(:,m)) * ww(:,l)) * rj(l,m)

         u0(jj(:,m)) = u0(jj(:,m)) + ww(:,l) * fl

      ENDDO
   ENDDO

END SUBROUTINE qs_00


SUBROUTINE qs_000 (mesh, ff1, ff2, u0)
!=================================

!  < w, f g >   ===>   u0

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff1, ff2 
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0

   INTEGER :: m, l, n
   REAL(KIND=8) :: fl, gl

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   u0 = 0

   DO m = 1, me
      DO l = 1, l_G

         fl = 0 
         gl = 0 
         DO n = 1, n_w
            fl = fl + ff1(jj(n,m)) * ww(n,l)
            gl = gl + ff2(jj(n,m)) * ww(n,l)
         END DO

         fl = fl * gl * rj(l,m)

         u0(jj(:,m)) = u0(jj(:,m)) + ww(:,l) * fl

      ENDDO
   ENDDO

END SUBROUTINE qs_000


!------------------------------------------------------------------------------

SUBROUTINE qs_01 (mesh, gg,  u0)
!=================================

!  < w, D.g >   ===>   u0

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0

   INTEGER ::  m, l
   REAL(KIND=8) :: dgl

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   u0 = 0

   DO m = 1, me
      DO l = 1, l_G

         dgl = SUM(gg(:,jj(:,m)) * dw(:,:,l,m)) * rj(l,m)

         u0(jj(:,m)) = u0(jj(:,m)) + ww(:,l) * dgl

      ENDDO
   ENDDO

END SUBROUTINE qs_01

SUBROUTINE qs_010 (mesh, gg, uu,  u0)
!=================================

!  < w, u D.g >   ===>   u0

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: uu 
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0

   INTEGER :: m, l, n, k
   REAL(KIND=8) :: dgl, ul

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   u0 = 0

   DO m = 1, me
      DO l = 1, l_G

         dgl = 0
         ul = 0 
         DO n = 1, n_w
            ul = ul + uu(jj(n,m)) * ww(n,l)
            DO k =1 , k_d
               dgl = dgl + gg(k,jj(n,m)) * dw(k,n,l,m)
            END DO
         END DO
         ul = ul * dgl * rj(l,m) 

         u0(jj(:,m)) = u0(jj(:,m)) + ww(:,l) * ul 

      ENDDO
   ENDDO

END SUBROUTINE qs_010

!------------------------------------------------------------------------------

SUBROUTINE qs_10 (mesh, gg,  u0)
!=================================

!  << (Dw), g >>   ===>   u0

   USE Gauss_points

   IMPLICIT NONE
 
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   INTEGER :: m, l, k, n
   REAL(KIND=8), DIMENSION(mesh%gauss%k_d) :: gl

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   u0 = 0

   DO m = 1, me
      DO l = 1, l_G

         gl =0.
         DO k = 1, k_d
            DO n =1 ,n_w
!            gl(k) = SUM(gg(k,jj(:,m)) * ww(:,l)) * rj(l,m)
               gl(k) =gl(k) + gg(k,jj(n,m)) * ww(n,l)
            END DO
         ENDDO
         gl = gl *rj(l,m)

         DO n = 1, n_w
            DO k = 1, k_d
               u0(jj(n,m)) = u0(jj(n,m)) + dw(k,n,l,m) * gl(k)
            END DO
         ENDDO

      ENDDO
   ENDDO

END SUBROUTINE qs_10

SUBROUTINE qs_01_dx (mesh, ff,  u0)
!=================================

!  << w, d_x g >>   ===>   u0

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff 
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0

   INTEGER :: m, l, n
   REAL(KIND=8)  :: dfl

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   u0 = 0

   DO m = 1, me
      DO l = 1, l_G

         dfl =0.
         DO n =1 ,n_w
            dfl = dfl + ff(jj(n,m)) * dw(1,n,l,m) 
         ENDDO
         dfl = dfl * rj(l,m)

         DO n = 1, n_w
            u0(jj(n,m)) = u0(jj(n,m)) + ww(n,l) * dfl 
         ENDDO

      ENDDO
   ENDDO

END SUBROUTINE qs_01_dx

    !=====================================================
  SUBROUTINE qs_h_f_v_grad ( mesh, gg, ff, s0)
    !=====================================================

    !  << f , (g.D)w, >>

    USE Gauss_points

    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! Declaration des variables globales
 
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)    :: gg
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)    :: ff
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: s0

    !--------------------------------------------------------------------------
    ! Declaration des variables locales
    INTEGER                        :: l, k, ni, m

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me

    REAL(KIND=8), DIMENSION(mesh%gauss%k_d)   :: gl
    REAL(KIND=8)                   :: fl, mest, hloc
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w)   :: y

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

    !--------------------------------------------------------------------------
    boucle_mm : DO m = 1, me

      mest = 0
      DO l= 1, l_G
         mest = mest + rj(l,m)
      END DO
      hloc = mest**(1.d0/k_d)

       boucle_l : DO l = 1, l_G

          gl = 0
          boucle_k : DO k = 1, k_d
             boucle_ni : DO ni =1 ,n_w
                gl(k) = gl(k) + gg(k, jj(ni,m)) * ww(ni,l)
             END DO boucle_ni
          ENDDO  boucle_k

          y = 0.
          boucle_ni_2 : DO ni = 1, n_w
             boucle_k_2 : DO k = 1, k_d
                y(ni) =  y(ni) + gl(k) * dw(k,ni,l,m)
             END DO boucle_k_2
          END DO  boucle_ni_2
          y = y*hloc
          ! y "contient" v.grad(w)

          fl = SUM(ff(jj(:,m)) * ww(:,l))
          s0(jj(:,m)) = s0(jj(:,m)) + y(:)*fl*rj(l,m)

       ENDDO  boucle_l

    ENDDO boucle_mm

  END SUBROUTINE qs_h_f_v_grad


!------------------------------------------------------------------------------

SUBROUTINE qs_11 (mesh, uu,  u0)
!=================================

!  << (Dw), Du) >>   ===>   u0

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: uu
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   INTEGER ::  m, l, k, n
   REAL(KIND=8), DIMENSION(mesh%gauss%k_d) :: gl

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   u0 = 0

   DO m = 1, me
      DO l = 1, l_G

         DO k = 1, k_d
            gl(k) = SUM(uu(jj(:,m)) * dw(k,:,l,m)) * rj(l,m)
         ENDDO

         DO n = 1, n_w
            u0(jj(n,m)) = u0(jj(n,m)) + SUM(dw(:,n,l,m) * gl)
         ENDDO

      ENDDO
   ENDDO

END SUBROUTINE qs_11

!------------------------------------------------------------------------------

SUBROUTINE qs_001 (mesh, gg, uu,  u0)
!======================================

!  < w, (g.D)u >   ===>   u0

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: uu
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0

   INTEGER :: m, l, k, n
   REAL(KIND=8) :: s, fl, dul, gl

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   u0 = 0

   DO m = 1, me
      DO l = 1, l_G

         s = 0
         DO k = 1, k_d
            gl = 0
            dul = 0
            DO n = 1, n_w
               gl = gl + gg(k,jj(n,m)) * ww(n,l)
               dul = dul + uu(jj(n,m)) * dw(k,n,l,m)
            END DO
            s = s  +  gl * dul
         ENDDO
         fl = s * rj(l,m)

         u0(jj(:,m)) = u0(jj(:,m)) + ww(:,l) * fl

      ENDDO
   ENDDO

END SUBROUTINE qs_001

!------------------------------------------------------------------------------

SUBROUTINE qs_100 (mesh, gg, uu,  u0)
!======================================

!  < (g.D)w, u >   ===>   u0

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: uu
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   INTEGER :: m, l, k, n
   REAL(KIND=8), DIMENSION(mesh%gauss%k_d) :: gul

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   u0 = 0

   DO m = 1, me
      DO l = 1, l_G

         DO k = 1, k_d
            gul(k) = SUM(gg(k,jj(:,m)) * ww(:,l))   &
                   * SUM(uu(jj(:,m)) * ww(:,l)) * rj(l,m)
         ENDDO

         DO n = 1, n_w
            u0(jj(n,m)) = u0(jj(n,m)) + SUM(dw(:,n,l,m) * gul)
         ENDDO

      ENDDO
   ENDDO

END SUBROUTINE qs_100

SUBROUTINE qs_100_anal (mesh, gg, uu,  u0)
!======================================

!  < (g.D)w, u >   ===>   u0

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: uu
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   INTEGER :: m, l, k, n, index
   REAL(KIND=8), DIMENSION(mesh%gauss%k_d) :: gul

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   u0 = 0

   DO m = 1, me
      index = (m-1)*l_G
      DO l = 1, l_G
         index = index + 1
         DO k = 1, k_d
            gul(k) = SUM(gg(k,jj(:,m)) * ww(:,l))   &
                   * uu(index) * rj(l,m)
         ENDDO

         DO n = 1, n_w
            u0(jj(n,m)) = u0(jj(n,m)) + SUM(dw(:,n,l,m) * gul)
         ENDDO

      ENDDO
   ENDDO

END SUBROUTINE qs_100_anal

!------------------------------------------------------------------------------

SUBROUTINE qs_GALS_mass_adv (mesh, param, alpha, gg, uu,  u0)
!======================================

!  < w+ h(alpha w + (g.D)w), u >   ===>   u0

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8),                 INTENT(IN)  :: alpha,param
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: uu
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   INTEGER :: m, l, k, n
   REAL(KIND=8), DIMENSION(mesh%gauss%k_d) :: gul
   REAL(KIND=8) :: ul, mest, hloc

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   u0 = 0

   DO m = 1, me

       mest = 0
       DO l = 1, l_G
          mest = mest + rj(l,m)
       END DO
       hloc = param*mest**(1.d0/k_d)
 
      DO l = 1, l_G
 
         ul = SUM(uu(jj(:,m)) * ww(:,l)) * rj(l,m)

         DO k = 1, k_d
            gul(k) = SUM(gg(k,jj(:,m)) * ww(:,l))
         ENDDO

         DO n = 1, n_w
            u0(jj(n,m)) = u0(jj(n,m)) + &
             ul * (ww(n,l) + hloc*(alpha*ww(n,l)+SUM(dw(:,n,l,m) * gul)))
         ENDDO

      ENDDO
   ENDDO

END SUBROUTINE qs_GALS_mass_adv

!------------------------------------------------------------------------------

SUBROUTINE qs_LS_mass_adv_anal (mesh, alpha, gg, uu,  u0)
!======================================

!  < alpha w + (g.D)w, u >   ===>   u0

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8),                 INTENT(IN)  :: alpha
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: uu
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   INTEGER :: m, l, k, n, index
   REAL(KIND=8), DIMENSION(mesh%gauss%k_d) :: gul
   REAL(KIND=8) :: ul

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   u0 = 0

   DO m = 1, me
      index = (m-1)*l_G
      DO l = 1, l_G
         index = index + 1
         DO k = 1, k_d
            gul(k) = SUM(gg(k,jj(:,m)) * ww(:,l))   &
                   * uu(index) * rj(l,m)
         ENDDO

         ul = alpha*rj(l,m) * uu(index)

         DO n = 1, n_w
            u0(jj(n,m)) = u0(jj(n,m)) + SUM(dw(:,n,l,m) * gul) &
                          + ul*ww(n,l)
         ENDDO

      ENDDO
   ENDDO

END SUBROUTINE qs_LS_mass_adv_anal
!------------------------------------------------------------------------------

SUBROUTINE qs_00_lump (mesh, alpha,  a_l)
!==========================================

!  alpha < w, _ >   ===>   a_l   diagonally lumped mass matrix 

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8),                 INTENT(IN)  :: alpha
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: a_l

   INTEGER :: m, l, n, i

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   a_l = 0

   DO m = 1, me
      DO l = 1, l_G

         DO n = 1, n_w;  i = jj(n, m)
            a_l(i) = a_l(i) + ww(n,l) * rj(l,m)
         ENDDO

      ENDDO
   ENDDO

   a_l = alpha * a_l 

END SUBROUTINE qs_00_lump

!------------------------------------------------------------------------------

SUBROUTINE qs_00_diagonal (mesh, alpha,  a_d)
!==============================================

!  alpha < w, _ >   ===>   a_d   diagonal of mass matrix 

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8),                 INTENT(IN)  :: alpha
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: a_d

   INTEGER :: m, l, n, i

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   a_d = 0

   DO m = 1, me
      DO l = 1, l_G

         DO n = 1, n_w;  i = jj(n, m)
            a_d(i) = a_d(i) + ww(n,l)**2 * rj(l,m)
         ENDDO

      ENDDO
   ENDDO

   a_d = alpha * a_d 

END SUBROUTINE qs_00_diagonal

!------------------------------------------------------------------------------

SUBROUTINE qs_diagonal (mesh, alpha, gg,  a_d)
!==============================================

!   < Dw, D_ > + alpha < w, _ > + < w, g.D_ > + < w, _ (D.g) >/2 
!                         ===>   a_d   diagonal of mass matrix 

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8),                 INTENT(IN)  :: alpha
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: a_d

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   INTEGER :: m, k, l, n, i
   REAL(KIND=8) :: dg, x, y, z
   REAL(KIND=8), DIMENSION(mesh%gauss%k_d) ::  gl


   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   a_d = 0

   DO m = 1, me
      DO l = 1, l_G

        dg = SUM(gg(:, jj(:,m)) * dw(:,:,l,m)) * rj(l,m)
	DO k = 1, k_d
           gl(k) = SUM(gg(k, jj(:,m)) * ww(:,l)) * rj(l,m)
        ENDDO

         DO n = 1, n_w;  i = jj(n, m)
            x = SUM(dw(:,n,l,m)**2) * rj(l,m)
            y = alpha * ww(n,l)**2 * rj(l,m)
            z = ww(n,l) * ( SUM(gl*dw(:,n,l,m)) + dg*ww(n,l)/2 )
            a_d(i) = a_d(i) + x + y + z
         ENDDO

      ENDDO
   ENDDO

END SUBROUTINE qs_diagonal

!------------------------------------------------------------------------------

SUBROUTINE qs_00_s (mesh, fs,  u0)
!========================================

!  < ws, fs >_s   ===>   u0    incremental accumulation of boundary terms

   USE Gauss_points

   IMPLICIT NONE
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: fs
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0

   INTEGER :: ms, ls
   REAL(KIND=8) :: fls

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: js, is
   INTEGER,                      POINTER       :: mes

   CALL gauss(mesh)
   js => mesh%jjs
   is => mesh%iis
   mes => mesh%me

   DO ms = 1, mes
      DO ls = 1, l_Gs

         fls = SUM(fs(is(:,ms)) * wws(:,ls)) * rjs(ls,ms)

         u0(js(:,ms)) = u0(js(:,ms)) + wws(:,ls) * fls

      ENDDO
   ENDDO

END SUBROUTINE qs_00_s

!------------------------------------------------------------------------------

SUBROUTINE qs_01_s (mesh, gs,  u0)
!========================================

!  < ws, n.gs >_s   ===>   u0    incremental accumulation of boundary terms

   USE Gauss_points

   IMPLICIT NONE
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gs
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: js, is
   INTEGER,                      POINTER       :: mes

   INTEGER :: ms, ls, k, ns
   REAL(KIND=8) :: x
   REAL(KIND=8), DIMENSION(mesh%gauss%k_d) :: gls

   CALL gauss(mesh)
   js => mesh%jjs
   is => mesh%iis
   mes => mesh%mes

   DO ms = 1, mes
      DO ls = 1, l_Gs

         DO k = 1, k_d
            gls(k) = SUM(gs(k, is(:,ms)) * wws(:,ls)) * rjs(ls,ms)
         ENDDO

         DO ns = 1, n_ws
            x = wws(ns,ls) * SUM(gls * rnorms(:,ls,ms))
            u0(js(ns,ms)) = u0(js(ns,ms)) + x
         ENDDO

      ENDDO
   ENDDO

END SUBROUTINE qs_01_s

!------------------------------------------------------------------------------

SUBROUTINE qs_011 (mesh, f1, f2,  u0)
!======================================

!  < w, (Df1).(Df2) >   ===>   u0

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: f1, f2
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0

   INTEGER :: m, l, k
   REAL(KIND=8) :: dfl

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   u0 = 0

   DO m = 1, me
      DO l = 1, l_G

         dfl = 0
         DO k = 1, k_d
            dfl = dfl + SUM(f1(jj(:,m)) * dw(k,:,l,m))   &
                      * SUM(f2(jj(:,m)) * dw(k,:,l,m))
         ENDDO
         dfl = dfl * rj(l,m)

         u0(jj(:,m)) = u0(jj(:,m)) + ww(:,l) * dfl

      ENDDO
   ENDDO

END SUBROUTINE qs_011

!------------------------------------------------------------------------------

SUBROUTINE qs_101 (mesh, f1, f2,  u0)
!======================================

!  << (Dw), f1 Df2 >>   ===>   u0

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: f1, f2
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   INTEGER :: m, l, k, n
   REAL(KIND=8), DIMENSION(mesh%gauss%k_d) :: df, fdfl

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   u0 = 0

   DO m = 1, me
      DO l = 1, l_G

         DO k = 1, k_d
            df(k) = SUM(f2(jj(:,m)) * dw(k,:,l,m))
         ENDDO
         fdfl = SUM(f1(jj(:,m)) * ww(:,l)) * df * rj(l,m)

         DO n = 1, n_w
            u0(jj(n,m)) = u0(jj(n,m)) + SUM(dw(:,n,l,m) * fdfl)
         ENDDO

      ENDDO
   ENDDO

END SUBROUTINE qs_101

!------------------------------------------------------------------------------

SUBROUTINE qs_01_st (mesh, ff,  u0)
!=====================================

!  < ws, (t.Df) >_s   ===>   u0    (incremental accumulation)

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0

   INTEGER :: ms, ls
   REAL(KIND=8) :: df

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: js
   INTEGER,                      POINTER       :: mes

   CALL gauss(mesh)
   js => mesh%jjs
   mes => mesh%me

   DO ms = 1, mes
      DO ls = 1, l_Gs

         df = SUM(ff(js(:,ms)) * dwps(1,:,ls,ms))

         u0(js(:,ms)) = u0(js(:,ms)) + wws(:,ls) * df

      ENDDO
   ENDDO

END SUBROUTINE qs_01_st

!------------------------------------------------------------------------------

SUBROUTINE qs_dif_mass_adv(mesh, alpha, gg, uu,  u0)
!====================================================

!  << (Dw), Du) >>  
!  +  alpha * < w, u >  
!  +  < w, (g.D)u >  +  < w, (D.g) u > / 2 
!  ===>   u0

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8),                 INTENT(IN)  :: alpha 
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: uu
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   INTEGER :: k, l, m, n
   REAL(KIND=8), DIMENSION(mesh%gauss%k_d) :: gl 
   REAL(KIND=8) :: dg, gdu, dgu, ul

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   u0 = 0

   DO m = 1, me
      DO l = 1, l_G

         ul = SUM(uu(jj(:,m)) * ww(:,l)) * rj(l,m)

         dg = SUM(gg(:, jj(:,m)) * dw(:,:,l,m))
         dgu = dg * ul 

         gdu = 0
         DO k = 1, k_d
            gl(k) = SUM(uu(jj(:,m)) * dw(k,:,l,m))
            gdu = gdu + SUM(gg(k, jj(:,m)) * ww(:,l)) * gl(k) 
         ENDDO

	 gl = gl * rj(l,m)
         gdu = gdu * rj(l,m)

         DO n = 1, n_w
            u0(jj(n,m)) = u0(jj(n,m)) + SUM(dw(:,n,l,m) * gl)
         ENDDO

         u0(jj(:,m)) = u0(jj(:,m)) + ww(:,l) * (alpha*ul  +  gdu + dgu/2)

      ENDDO
   ENDDO

END SUBROUTINE qs_dif_mass_adv

SUBROUTINE bs_101 (mesh, ff, uu,  u0)
!======================================

! -  < J(w,u), f) >   ===>   u0
!   < (Dw) x k, f (Du) >   ===>   u0

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff, uu
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0

   INTEGER ::  m, l
   REAL(KIND=8) :: f
   REAL(KIND=8), DIMENSION(2) :: du

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   u0 = 0

   DO m = 1, me
      DO l = 1, l_G

         du(1) = SUM(uu(jj(:,m)) * dw(1,:,l,m))
         du(2) = SUM(uu(jj(:,m)) * dw(2,:,l,m))
             f = SUM(ff(jj(:,m)) * ww(:,l)) * rj(l,m)

         u0(jj(:,m)) = u0(jj(:,m))   &
                     - (dw(1,:,l,m) * du(2) - dw(2,:,l,m) * du(1)) * f

      ENDDO
   ENDDO

END SUBROUTINE bs_101


SUBROUTINE qsv_wave (mesh, f1, g1, g2, f2, u0, v0)
!=================================

!   < w, f1 + D.g1 >    ===>   u0
!  << w, g2 + D f2 >>   ===>   v0

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: f1, f2 
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: g1, g2 
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0
   REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: v0

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   INTEGER :: m, l, k, n
   REAL(KIND=8) :: f1l, dg1l
   REAL(KIND=8), DIMENSION(mesh%gauss%k_d) :: df2l, g2l
   REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w) :: dw_loc 

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   u0 = 0;  v0 = 0

   DO m = 1, me
      DO l = 1, l_G

         dw_loc = dw(:,:,l,m)

         f1l = 0;  dg1l = 0
         g2l = 0;  df2l = 0

         DO n = 1, n_w
            f1l = f1l + f1(jj(n,m)) * ww(n,l)
            DO k = 1, k_d
               dg1l    = dg1l    + g1(k,jj(n,m)) * dw_loc(k,n) 
                g2l(k) =  g2l(k) + g2(k,jj(n,m)) * ww(n,l)  
               df2l(k) = df2l(k) + f2(jj(n,m))   * dw_loc(k,n) 
            END DO
         END DO
         f1l = (f1l + dg1l) * rj(l,m)
         g2l = (g2l + df2l) * rj(l,m)
         
         u0(jj(:,m)) = u0(jj(:,m)) + ww(:,l) * f1l
         DO k = 1, k_d
            v0(k,jj(:,m)) = v0(k,jj(:,m)) + ww(:,l) * g2l(k) 
         END DO

      ENDDO
   ENDDO

END SUBROUTINE qsv_wave

SUBROUTINE qs_00_anal (mesh, ff,  u0)
!=================================

!  < w, f >   ===>   u0

   USE Gauss_points

   IMPLICIT NONE
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0

   INTEGER :: m, l, index

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   u0 = 0

   index = 0
   DO m = 1, me
      DO l = 1, l_G; index = index + 1

         u0(jj(:,m)) = u0(jj(:,m)) + ww(:,l) * ff(index) * rj(l,m)

      ENDDO
   ENDDO

END SUBROUTINE qs_00_anal

SUBROUTINE qs_1_sgs_1 (mesh, f1, f2, u0)
!=======================================

!  << (Dw), h* f1 (Df2) >>   ===>   u0    incremental accumulation

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: f1, f2
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0 

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   INTEGER ::  m, l, ni, i, k
   REAL(KIND=8) :: fl, h, exp, mest
   REAL(KIND=8), DIMENSION(mesh%gauss%k_d) :: df 

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   m = SIZE(jj,1)
   SELECT CASE(m)
      CASE(3); exp = 1.d0/2
      CASE(6); exp = 1.d0/2
      CASE(4); exp = 1.d0/3
      CASE(10); exp = 1.d0/3
   END SELECT

   u0 =0

   DO m = 1, me

      mest = 0
      DO l = 1, l_G
         mest = mest + rj(l,m)
      END DO
      h = mest**exp
!TEST : Ca ne marche pas bien
!      fnorm = 0 
!      DO l = 1, l_G
!         fl = 0
!         DO ni =1 ,n_w
!            fl = fl + f1(jj(ni,m)) * ww(ni,l)
!         END DO
!         fnorm = fnorm + ABS(fl)*rj(l,m)
!      END DO
!      fnorm  = fnorm/mest
!TEST

      DO l = 1, l_G
!TEST
         fl = 0
         DO ni =1 ,n_w
            fl = fl + f1(jj(ni,m)) * ww(ni,l)
         END DO
!TEST
         fl  = h*fl*rj(l,m)

         df = 0
         DO k = 1, k_d
            DO ni =1 ,n_w
               df(k) = df(k) + f2(jj(ni,m)) * dw(k,ni,l,m)
            END DO
         END DO 

         DO ni = 1, n_w; i = jj(ni, m)
            DO k = 1, k_d
                  u0(i) = u0(i) + dw(k,ni,l,m) * df(k) * fl 
            END DO
         END DO

      ENDDO
   ENDDO

END SUBROUTINE qs_1_SGS_1

SUBROUTINE qs_1_sgs_1_norm (mesh, f1, f2, u0)
!=======================================

!  << (Dw), h* f1 (Df2) >>   ===>   u0    incremental accumulation

   USE Gauss_points

   IMPLICIT NONE
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: f1, f2
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0 

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   INTEGER :: m, l, ni, i, k
   REAL(KIND=8) :: fl, h, exp, mest, fnorm
   REAL(KIND=8), DIMENSION(mesh%gauss%k_d) :: df 

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   m = SIZE(jj,1)
   SELECT CASE(m)
      CASE(3); exp = 1.d0/2
      CASE(6); exp = 1.d0/2
      CASE(4); exp = 1.d0/3
      CASE(10); exp = 1.d0/3
   END SELECT

   u0 =0

   DO m = 1, me

      mest = 0
      DO l = 1, l_G
         mest = mest + rj(l,m)
      END DO
      h = mest**exp

      fnorm = 0 
      DO l = 1, l_G
         fl = 0
         DO ni =1 ,n_w
            fl = fl + f2(jj(ni,m)) * ww(ni,l)
         END DO
         fnorm = fnorm + ABS(fl)*rj(l,m)
      END DO
      fnorm  = fnorm/mest
      IF (fnorm .LE. 1d-10) fnorm = 1

      DO l = 1, l_G
         fl = 0
         DO ni =1 ,n_w
            fl = fl + f1(jj(ni,m)) * ww(ni,l)
         END DO
         fl  = h*fl*rj(l,m)/fnorm

         df = 0
         DO k = 1, k_d
            DO ni =1 ,n_w
               df(k) = df(k) + f2(jj(ni,m)) * dw(k,ni,l,m)
            END DO
         END DO 

         DO ni = 1, n_w; i = jj(ni, m)
            DO k = 1, k_d
                  u0(i) = u0(i) + dw(k,ni,l,m) * df(k) * fl 
            END DO
         END DO

      ENDDO
   ENDDO

END SUBROUTINE qs_1_SGS_1_norm

SUBROUTINE qs_1_sgs_1_test (mesh, expo, f1, f2, u0)
!=======================================

!  << (Dw), h* f1 (Df2) >>   ===>   u0    incremental accumulation

   USE Gauss_points

   IMPLICIT NONE
   REAL(KIND=8),                 INTENT(IN)  :: expo
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: f1, f2
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0 

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   INTEGER :: mm, m, mloc, l, ni, i, k
   REAL(KIND=8) :: fl, h, exp, mest, fnorm, grad
   REAL(KIND=8), DIMENSION(mesh%gauss%k_d) :: df 

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   m = SIZE(jj,1)
   SELECT CASE(m)
      CASE(3); exp = 1.d0/2
      CASE(6); exp = 1.d0/2
      CASE(4); exp = 1.d0/3
      CASE(10); exp = 1.d0/3
   END SELECT

   u0 =0

   DO mm = 1, me/4

!   mest = 0
!   fnorm = 0 
!   DO mloc = 1, 4
!      m = 4*(mm-1) + mloc !  mmg(mloc,mg)  !Global index of element
!      DO l = 1, l_G
!         fl = 0
!         DO ni =1 ,n_w
!            fl = fl + f2(jj(ni,m)) * ww(ni,l)
!         END DO
!         fnorm = fnorm + ABS(fl)*rj(l,m)
!         mest = mest + rj(l,m)
!      END DO
!   END DO
!   fnorm  = fnorm/mest
!   IF (fnorm .LE. 1d-10) fnorm = 1
!   H = mest**exp


   DO mloc = 1, 4
      m = 4*(mm-1) + mloc !  mmg(mloc,mg)  !Global index of element

      mest = 0
      fnorm = 0 
      DO l = 1, l_G
         fl = 0
         DO ni =1 ,n_w
            fl = fl + f2(jj(ni,m)) * ww(ni,l)
         END DO
         fnorm = fnorm + ABS(fl)*rj(l,m)
         mest = mest + rj(l,m)
      END DO
      fnorm  = fnorm/mest
      IF (fnorm .LE. 1d-10) fnorm = 1
      h = mest**exp

      DO l = 1, l_G
         fl = 0
         DO ni =1 ,n_w
            fl = fl + f1(jj(ni,m)) * ww(ni,l)
         END DO

         grad = 0
         df = 0
         DO k = 1, k_d
            DO ni =1 ,n_w
               df(k) = df(k) + f2(jj(ni,m)) * dw(k,ni,l,m)
            END DO
         grad = grad + df(k)**2
         END DO 
         grad = SQRT(grad)

!BON fl = min((ABS(fl)/fnorm)**expo, 1.d0)
!if((ABS(fl)/fnorm)**expo .ge. 1) then
!beta = 1.d0 
!if(ABS(fl)**beta*(grad)**expo .ge. 1) then
!     fl = 1.d0
!else 
!     fl = 0
!end if 
fl = ABS(fl)*expo
!  fl = ABS(fl/fnorm)**expo
! fl = ABS(fl)/fnorm**expo
!TEST
!MAUVAIS  fl = (ABS(fl)/fnorm)**expo
!BOF      fl = min((ABS(fl)*SQRT(grad)/fnorm)**expo, 1.d0)
!BOF     fl = min(((ABS(fl)/fnorm)*(h*grad/fnorm))**expo, 1.d0)
!   beta = 1.d0 
!   fl = min(((ABS(fl/norm))**beta)*((h*grad/norm)**expo),1.d0)
!        fl  = -h*fl**(1./expo)*rj(l,m)
!TEST

         fl  = -h*fl*rj(l,m)

         DO ni = 1, n_w; i = jj(ni, m)
            DO k = 1, k_d
                  u0(i) = u0(i) + dw(k,ni,l,m) * df(k) * fl 
            END DO
         END DO

      ENDDO
   ENDDO
   ENDDO

END SUBROUTINE qs_1_SGS_1_test

SUBROUTINE qs_mass_adv(mesh, gg, ff, uu,  v0)
!=================================

!  << w, g >>
!  << w, u.Df >>
!     ===>   v0

   USE Gauss_points

   IMPLICIT NONE
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: uu
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: gg, ff
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: v0

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   INTEGER ::  m, l, k, n
   INTEGER, DIMENSION(mesh%gauss%n_w)          :: j_loc
   REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w) :: uukn
   REAL(KIND=8), DIMENSION(mesh%gauss%n_w)     :: ffn, ggn

   REAL(KIND=8) :: u, df, udf, g, resul

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   v0 = 0

   DO m = 1, me

      j_loc(:)  = jj(:,m)
      uukn(:,:) = uu(:,jj(:,m))
      ggn(:)    = gg(jj(:,m))
      ffn(:)    = ff(jj(:,m))

      DO l = 1, l_G

         udf = 0
         DO k = 1, k_d
            u = 0
            df = 0
            DO n = 1, n_w
               u  = u  + uukn(k,n) * ww(n,l)
               df = df + ffn(n)     * dw(k,n,l,m)
            END DO
            udf = udf + u*df
         ENDDO


         g = 0
         DO n = 1, n_w
            g =  g + ggn(n) * ww(n,l)
         END DO

         resul = (g + udf) * rj(l,m)

         DO n = 1, n_w
            v0(j_loc(n)) = v0(j_loc(n)) + ww(n,l) * resul
         ENDDO


      ENDDO
   ENDDO

END SUBROUTINE qs_mass_adv


!------------------------------------------------------------------------------
SUBROUTINE qs_01_hybrid (uu_mesh, pp_mesh, gg,  u0_c)
!================================================

!  < w_c, D.g >   ===>   u0_c

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0_c

   TYPE(mesh_type), TARGET                     :: uu_mesh, pp_mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj, jj_c
   INTEGER,      POINTER       :: me, nwc

   REAL(KIND=8), DIMENSION(pp_mesh%gauss%n_w, uu_mesh%gauss%l_G) ::  w_c

   INTEGER :: m, l, n, k
   REAL(KIND=8) :: dgl
   REAL(KIND=8), DIMENSION(uu_mesh%gauss%k_d,uu_mesh%gauss%n_w) :: dwkn, ggkn 
   REAL(KIND=8), DIMENSION(uu_mesh%gauss%n_w) :: u0_cn

   CALL gauss(uu_mesh)
   jj => uu_mesh%jj
   jj_c => pp_mesh%jj
   me => uu_mesh%me
   nwc => pp_mesh%gauss%n_w

   SELECT CASE(k_d)

      CASE(2)
         DO l = 1, l_G       
!         w_c(1,:) = ww(1,:) + 0.5*(ww(5,:) + ww(6,:)) 
!         w_c(2,:) = ww(2,:) + 0.5*(ww(6,:) + ww(4,:)) 
!         w_c(3,:) = ww(3,:) + 0.5*(ww(4,:) + ww(5,:)) 
         w_c(1,l) = ww(1,l) + 0.5*(ww(n_w-1,l) + ww(n_w,l)) 
         w_c(2,l) = ww(2,l) + 0.5*(ww(n_w,l) + ww(4,l)) 
         w_c(3,l) = ww(3,l) + 0.5*(ww(4,l) + ww(n_w-1,l)) 
         END DO  

      CASE(3)
         DO l = 1, l_G       
!         w_c(1,:) = ww(1,:) + 0.5*(ww(n_w-2,:) + ww(n_w-1,:) + ww(n_w,:)) 
!         w_c(2,:) = ww(2,:) + 0.5*(ww(6,:)     + ww(n_w-3,:) + ww(n_w,:))
!         w_c(3,:) = ww(3,:) + 0.5*(ww(5,:)     + ww(n_w-3,:) + ww(n_w-1,:)) 
!         w_c(4,:) = ww(4,:) + 0.5*(ww(5,:)     + ww(6,:)     + ww(n_w-2,:)) 
         w_c(1,l) = ww(1,l) + 0.5*(ww(n_w-2,l)  + ww(n_w-1,l) + ww(n_w,l)) 
         w_c(2,l) = ww(2,l) + 0.5*(ww(n_w-4,l)  + ww(n_w-3,l) + ww(n_w,l))
         w_c(3,l) = ww(3,l) + 0.5*(ww(n_w-5,l)  + ww(n_w-3,l) + ww(n_w-1,l)) 
         w_c(4,l) = ww(4,l) + 0.5*(ww(n_w-5,l)  + ww(n_w-4,l) + ww(n_w-2,l)) 
         END DO  

   END SELECT

   u0_c = 0

   DO m = 1, me

      ggkn(:,:) = gg(:,jj(:,m))

      u0_cn = 0.

      DO l = 1, l_G
         dwkn(:,:) = dw(:,:,l,m)

         dgl = 0.
         DO n = 1, n_w
            DO k = 1, k_d
               dgl = dgl + ggkn(k,n) * dwkn(k,n)
            END DO
         END DO

! FINIR d'optimiser la boucle.

         dgl = dgl * rj(l,m)

         DO n = 1, nwc 
           u0_cn(n) = u0_cn(n) + w_c(n,l) * dgl
         END DO

      ENDDO

!         u0_c(jj_c(:,m)) = u0_c(jj_c(:,m)) + u0_cn(:)
         DO n = 1, nwc 
            u0_c(jj_c(n,m)) = u0_c(jj_c(n,m)) + u0_cn(n)
         END DO

   ENDDO

END SUBROUTINE qs_01_hybrid

END MODULE fem_s
