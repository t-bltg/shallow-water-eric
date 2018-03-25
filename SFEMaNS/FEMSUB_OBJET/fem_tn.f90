!   ===========================================
!   Entry points for TOTAL quantities and NORMS
!   ===========================================
!
!   ff  and  u0  are Scalar functions  (fs is defined only on the boundary)
!   gg  and  v0  are Vector fields     (gs is defined only on the boundary)
!
!   w   denotes either a scalar of a vector weigthing function
!   ws  denotes either a scalar of a vector weigthing function
!       defined only on the boundary
!
!   < _ , _ >     means integration over elements in the list  m0
!  << _ , _ >>    means scalar product and integration over elements in  m0
!   < _ , _ >_s   means surface integration over boundary elements in  ms0
!  << _ , _ >>_s  means scalar product and surface integration over b.e. ms0
!
!   D  is Nabla spatial derivative operator (grad, div, curl)
!
MODULE  fem_tn

CONTAINS

SUBROUTINE ts_0 (mesh, ff,  t)
!===============================

!  < f >   ===>   t

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
   REAL(KIND=8),                 INTENT(OUT) :: t

   INTEGER :: m, l, n
   REAL(KIND=8) ::  fl

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   t = 0

   DO m = 1, me
      DO l = 1, l_G
         fl = 0
         DO n = 1, n_w
            fl = fl + ff(jj(n,m)) * ww(n,l)
         END DO
         t = t + fl * rj(l,m)

      ENDDO
   ENDDO

END SUBROUTINE ts_0

!------------------------------------------------------------------------------

SUBROUTINE tv_1 (mesh, gg,  t)
!===============================

!  < D.g >   ===>   t

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8),                 INTENT(OUT) :: t

   INTEGER ::  m, l

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   t = 0

   DO m = 1, me
      DO l = 1, l_G

         t = t + SUM(gg(:,jj(:,m)) * dw(:,:,l,m)) * rj(l,m)

      ENDDO
   ENDDO

END SUBROUTINE tv_1

!------------------------------------------------------------------------------

SUBROUTINE ts_0_s (mesh, fs,  t)
!======================================

!  < fs >_s   ===>   t

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: fs
   REAL(KIND=8),                 INTENT(OUT) :: t

   INTEGER :: ms, ls

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: is
   INTEGER,                      POINTER       :: mes

   CALL gauss(mesh)
   is => mesh%iis
   mes => mesh%mes

   t = 0

   DO ms = 1, mes
      DO ls = 1, l_Gs

         t = t + SUM(fs(is(:,ms)) * wws(:,ls)) * rjs(ls,ms)

      ENDDO
   ENDDO

END SUBROUTINE ts_0_s

!------------------------------------------------------------------------------

SUBROUTINE ts_0_s_G (mesh, fs_G,  t)
!======================================

!  < fs >_s   ===>   t

!  fs_G is defined at Gauss points of the boundary elements

   USE Gauss_points

   IMPLICIT NONE
 
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: fs_G
   REAL(KIND=8),                 INTENT(OUT) :: t

   INTEGER ::  ms, ls

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,                      POINTER       :: mes

   CALL gauss(mesh)
   mes => mesh%mes

   t = 0

   DO ms = 1, mes
      DO ls = 1, l_Gs

         t = t + fs_G(ls,ms) * rjs(ls,ms)

      ENDDO
   ENDDO

END SUBROUTINE ts_0_s_G

!------------------------------------------------------------------------------

SUBROUTINE tv_0_s (mesh, gs,  t)
!======================================

!  < n.gs >_s   ===>   t

   USE Gauss_points

   IMPLICIT NONE
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gs
   REAL(KIND=8),                 INTENT(OUT) :: t

   INTEGER :: ms, ls, k

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: is
   INTEGER,                      POINTER       :: mes

   REAL(KIND=8), DIMENSION(k_d) :: gls

   CALL gauss(mesh)
   is => mesh%iis
   mes => mesh%mes

   t = 0

   DO ms = 1, mes
      DO ls = 1, l_Gs

         DO k = 1, k_d
            gls(k) = SUM(gs(k, is(:,ms)) * wws(:,ls)) * rjs(ls,ms)
         ENDDO

         t = t + SUM(gls * rnorms(:,ls,ms))

      ENDDO
   ENDDO

END SUBROUTINE tv_0_s

SUBROUTINE couple_p1(mesh, uu, t)   
!======================================

! <u.Du.n>_s ==> t 

  USE Gauss_points

  IMPLICIT NONE


   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: uu 
   REAL(KIND=8),                 INTENT(OUT) :: t

   REAL(KIND=8), DIMENSION(k_d,k_d) :: du 
   REAL(KIND=8), DIMENSION(k_d)     :: u 
   REAL(KIND=8)                     :: s
   INTEGER ::  ms, ls, m, n, i, k, kp, ns

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj, jjs
   INTEGER,      DIMENSION(:),   POINTER       :: neighs
   INTEGER,                      POINTER       :: mes

   CALL gauss(mesh)
   jj => mesh%jj
   jjs => mesh%jjs
   neighs => mesh%neighs
   mes => mesh%mes

   t = 0

   DO ms = 1, mes

      m = neighs(ms)

      DO ls = 1, l_Gs

         du = 0
         DO n = 1, n_w;  i = jj(n,m)
            DO k = 1, k_d
               DO kp = 1, k_d
                  du(k,kp) = du(k,kp) + uu(k,i) * dw(kp,n,1,m)
               END DO
            END DO
         END DO
 
         u = 0
         DO ns = 1, n_ws;  i = jjs(ns,ms)
            DO k = 1, k_d
               u(k) = u(k) + uu(k,i) * wws(ns,ls)
            END DO
         END DO

        s = 0
        DO k = 1, k_d
           DO kp = 1, k_d
              s = s + u(k)*du(k,kp)*rnorms(kp,ls,ms) 
           END DO
       END DO

       t = t + s * rjs(ls,ms)

      END DO
   END DO

END SUBROUTINE couple_p1

!------------------------------------------------------------------------------

SUBROUTINE tb_0_s (mesh, gs,  t)
!======================================

!  < tau_s.gs >_s   ===>   t       ( 2d only )

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gs
   REAL(KIND=8),                 INTENT(OUT) :: t

   INTEGER :: ms, ls, k
   REAL(KIND=8) :: s
   REAL(KIND=8), DIMENSION(2) :: gls

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: is
   INTEGER,                      POINTER       :: mes

   CALL gauss(mesh)
   is => mesh%iis
   mes => mesh%mes

   t = 0

   DO ms = 1, mes
      DO ls = 1, l_Gs

         DO k = 1, 2
            gls(k) = SUM(gs(k,is(:,ms)) * wws(:,ls)) * rjs(ls, ms)
         ENDDO
         s = - gls(1) * rnorms(2,ls,ms) + gls(2) * rnorms(1,ls,ms)

         t = t + s

      ENDDO
   ENDDO

END SUBROUTINE tb_0_s

!------------------------------------------------------------------------------

SUBROUTINE ns_0 (mesh, ff,  t)
!===============================

!  sqrt(< f^2 >)   ===>   t

   USE Gauss_points

   IMPLICIT NONE
 
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
   REAL(KIND=8),                 INTENT(OUT) :: t

   INTEGER ::  m, l

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   t = 0

   DO m = 1, me
      DO l = 1, l_G

         t = t + SUM(ff(jj(:,m)) * ww(:,l))**2 * rj(l,m)

      ENDDO
   ENDDO

   t = SQRT(t)

END SUBROUTINE ns_0

!------------------------------------------------------------------------------

SUBROUTINE ns_1 (mesh, ff,  t)
!===============================

!  sqrt(<< (Df).(Df) >>)   ===>   t

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
   REAL(KIND=8),                 INTENT(OUT) :: t

   INTEGER ::  m, l, k
   REAL(KIND=8) :: s

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   t = 0

   DO m = 1, me
      DO l = 1, l_G

         s = 0
         DO k = 1, k_d
            s = s + SUM(ff(jj(:,m)) * dw(k,:,l,m))**2
         ENDDO

         t = t + s * rj(l,m)

      ENDDO
   ENDDO

   t = SQRT(t)

END SUBROUTINE ns_1

SUBROUTINE nh_1 (mesh, ff,  t)
!===============================

! sqrt(<f.f> + << (Df).(Df) >>)   ===>   t

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
   REAL(KIND=8),                 INTENT(OUT) :: t

   INTEGER ::  m, l, k
   REAL(KIND=8) :: s

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   t = 0

   DO m = 1, me
      DO l = 1, l_G

         s = SUM(ff(jj(:,m)) * ww(:,l))**2
         DO k = 1, k_d
            s = s + SUM(ff(jj(:,m)) * dw(k,:,l,m))**2
         ENDDO

         t = t + s * rj(l,m)

      ENDDO
   ENDDO

   t = SQRT(t)

END SUBROUTINE nh_1

SUBROUTINE ns_w11 (mesh, ff,  t)
!===============================

! (<ABS(f)> + << ABS(Df)>>)   ===>   t

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
   REAL(KIND=8),                 INTENT(OUT) :: t

   INTEGER ::  m, l, k
   REAL(KIND=8) :: s

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   t = 0

   DO m = 1, me
      DO l = 1, l_G

         s = ABS(SUM(ff(jj(:,m)) * ww(:,l)))
         DO k = 1, k_d
            s = s + ABS(SUM(ff(jj(:,m)) * dw(k,:,l,m)))
         ENDDO

         t = t + s * rj(l,m)

      ENDDO
   ENDDO


END SUBROUTINE ns_w11

!------------------------------------------------------------------------------

SUBROUTINE nv_0 (mesh, gg,  t)
!===============================

!  sqrt(<< g.g >>)   ===>   t

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8),                 INTENT(OUT) :: t

   INTEGER ::  m, l, k
   REAL(KIND=8) :: s

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   t = 0

   DO m = 1, me
      DO l = 1, l_G

         s = 0
         DO k = 1, k_d
            s = s + SUM(gg(k,jj(:,m)) * ww(:,l))**2
         ENDDO

         t = t + s * rj(l,m)

      ENDDO
   ENDDO

   t = SQRT(t)

END SUBROUTINE nv_0

!------------------------------------------------------------------------------

SUBROUTINE nv_1 (mesh, gg,  t)
!===============================

!  sqrt(<< D.g, D.g >>)   ===>   t

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8),                 INTENT(OUT) :: t

   INTEGER ::  m, l, k, n
   REAL(KIND=8) :: s

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   t = 0

   DO m = 1, me
      DO l = 1, l_G

         s = 0
         DO k = 1, k_d
            DO n = 1, n_w
               s = s + gg(k,jj(n,m)) * dw(k,n,l,m)
            END DO
         END DO
         t = t + rj(l,m) * (s**2)

      ENDDO
   ENDDO

   t = SQRT(t)

END SUBROUTINE nv_1

SUBROUTINE nv_1_bloc (mesh, gg,  t)
!===============================

!  sqrt(<< D.g, D.g >>)   ===>   t

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: gg
   REAL(KIND=8),                 INTENT(OUT) :: t

   INTEGER ::  m, l, k, n, bloc_size
   REAL(KIND=8) :: s

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   bloc_size = SIZE(gg)/k_d

   t = 0

   DO m = 1, me
      DO l = 1, l_G

         s = 0
         DO k = 1, k_d
            DO n = 1, n_w
               s = s + gg((k-1)*bloc_size+jj(n,m)) * dw(k,n,l,m)
            END DO
         END DO
         t = t + rj(l,m) * (s**2)

      ENDDO
   ENDDO

   t = SQRT(t)

END SUBROUTINE nv_1_bloc

SUBROUTINE nv_l1 (mesh, gg,  t)
!===============================

!  << |D.g| >>   ===>   t

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8),                 INTENT(OUT) :: t

   INTEGER ::  m, l, k, n
   REAL(KIND=8) :: s

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   t = 0

   DO m = 1, me
      DO l = 1, l_G

         s = 0
         DO k = 1, k_d
            DO n = 1, n_w
               s = s + gg(k,jj(n,m)) * dw(k,n,l,m)
            END DO
         END DO
         t = t + rj(l,m) * ABS(s)

      ENDDO
   ENDDO

END SUBROUTINE nv_l1
!------------------------------------------------------------------------------

SUBROUTINE nc_1 (mesh, gg,  t)
!===============================

!  sqrt(<< (D x g).(D x g) >>)   ===>   t

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8),                 INTENT(OUT) :: t

   INTEGER ::  m, l, k, k1, k2
   REAL(KIND=8) :: s

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   t = 0

   DO m = 1, me
      DO l = 1, l_G

         s = 0
         DO k = 1, k_d;  k1 = MODULO(k, k_d) + 1;  k2 = MODULO(k + 1, k_d) + 1
            s = s + SUM(gg(k2,jj(:,m)) * dw(k1,:,l,m)   &
                      - gg(k1,jj(:,m)) * dw(k2,:,l,m))**2
         ENDDO

         t = t + s * rj(l,m)

      ENDDO
   ENDDO

   t = SQRT(t)

END SUBROUTINE nc_1

!------------------------------------------------------------------------------

SUBROUTINE nb_1 (mesh, gg,  t)
!===============================

!  sqrt(< (D x g . k) (D x g . k) >)   ===>   t     ( 2d only )

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8),                 INTENT(OUT) :: t

   INTEGER ::  m, l
   REAL(KIND=8) :: s

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   t = 0

   DO m = 1, me
      DO l = 1, l_G

         s = SUM(gg(2,jj(:,m)) * dw(1,:,l,m)   &
               - gg(1,jj(:,m)) * dw(2,:,l,m))**2

         t = t + s * rj(l,m)

      ENDDO
   ENDDO

   t = SQRT(t)

END SUBROUTINE nb_1

SUBROUTINE nb_1_bloc (mesh, gg,  t)
!===============================

!  sqrt(< (D x g . k) (D x g . k) >)   ===>   t     ( 2d only )

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: gg
   REAL(KIND=8),                 INTENT(OUT) :: t

   INTEGER ::  m, l, bloc_size
   REAL(KIND=8) :: s

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me
   bloc_size = SIZE(gg)/k_d
   t = 0

   DO m = 1, me
      DO l = 1, l_G

         s = SUM(gg(bloc_size+jj(:,m)) * dw(1,:,l,m)   &
               - gg(jj(:,m)) * dw(2,:,l,m))**2

         t = t + s * rj(l,m)

      ENDDO
   ENDDO

   t = SQRT(t)

END SUBROUTINE nb_1_bloc
!------------------------------------------------------------------------------

SUBROUTINE tb_00_s (mesh, fs, gs,  t)
!===========================================

!  < fs tau_s.gs >_s   ===>   t     ( 2d only )

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: fs
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gs
   REAL(KIND=8),                 INTENT(OUT) :: t

   INTEGER ::  ms, ls, k
   REAL(KIND=8) :: fls, s
   REAL(KIND=8), DIMENSION(2) :: gls

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: is
   INTEGER,                      POINTER       :: mes

   CALL gauss(mesh)
   is => mesh%iis
   mes => mesh%mes

   t = 0

   DO ms = 1, mes
      DO ls = 1, l_Gs

         fls = SUM(fs(is(:,ms)) * wws(:,ls)) * rjs(ls, ms)
         DO k = 1, 2
            gls(k) = SUM(gs(k,is(:,ms)) * wws(:,ls))
         ENDDO
         s = fls * ( - gls(1) * rnorms(2,ls,ms) + gls(2) * rnorms(1,ls,ms) )

         t = t + s

      ENDDO
   ENDDO

END SUBROUTINE tb_00_s

!------------------------------------------------------------------------------

SUBROUTINE ts_00 (mesh, f1, f2,  t)
!====================================

!  < f1 f2 >   ===>   t

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: f1, f2
   REAL(KIND=8),                 INTENT(OUT) :: t

   INTEGER :: m, l

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   t = 0

   DO m = 1, me
      DO l = 1, l_G

         t = t + SUM(f1(jj(:,m)) * ww(:,l))   &
               * SUM(f2(jj(:,m)) * ww(:,l)) * rj(l,m)

      ENDDO
   ENDDO

END SUBROUTINE ts_00

!------------------------------------------------------------------------------

SUBROUTINE ts_11 (mesh, f1, f2,  t)
!====================================

!  << (Df1).(Df2) >>   ===>   t

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: f1, f2
   REAL(KIND=8),                 INTENT(OUT) :: t

   INTEGER ::  m, l, k
   REAL(KIND=8) :: s

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   t = 0

   DO m = 1, me
      DO l = 1, l_G

         s = 0
         DO k = 1, k_d
            s = s + SUM(f1(jj(:,m)) * dw(k,:,l,m))   &
                  * SUM(f2(jj(:,m)) * dw(k,:,l,m))
         ENDDO

         t = t + s * rj(l,m)

      ENDDO
   ENDDO

END SUBROUTINE ts_11

!------------------------------------------------------------------------------

SUBROUTINE ns_adv_l1 (mesh, gg, ff,  t)
!====================================

!  SQRT <<|g.(Df)| >>   ===>   t

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
   REAL(KIND=8),                 INTENT(OUT) :: t

   REAL(KIND=8), DIMENSION(k_d,n_w)          :: dw_loc, g_loc 
   REAL(KIND=8), DIMENSION(n_w)              :: f_loc 
   INTEGER ::  m, l, k, n
   REAL(KIND=8) :: gk, dfk, gdf 

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   t = 0

   DO m = 1, me

      g_loc = gg(:, jj(:,m))  
      f_loc = ff(jj(:,m))  

      DO l = 1, l_G

         dw_loc = dw(:,:,l,m)

         gdf = 0
         DO k = 1, k_d
            gk = 0; dfk =0
            DO n = 1, n_w 
               gk = gk + g_loc(k,n) * ww(n,l) 
               dfk = dfk + f_loc(n) * dw_loc(k,n)
            END DO
            gdf = gdf + gk*dfk 
         ENDDO

         t = t + ABS(gdf) * rj(l,m)

      ENDDO
   ENDDO

END SUBROUTINE ns_adv_l1

SUBROUTINE ns_adv (mesh, gg, ff,  t)
!====================================

!  SQRT << ff^2 + (g.(Df))^2 >>   ===>   t

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
   REAL(KIND=8),                 INTENT(OUT) :: t

   REAL(KIND=8), DIMENSION(k_d,n_w)          :: dw_loc, g_loc 
   REAL(KIND=8), DIMENSION(n_w)              :: f_loc 
   INTEGER ::  m, l, k, n
   REAL(KIND=8) :: fl, gk, dfk, gdf 

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   t = 0

   DO m = 1, me

      g_loc = gg(:, jj(:,m))  
      f_loc = ff(jj(:,m))  

      DO l = 1, l_G

         dw_loc = dw(:,:,l,m)


         fl = 0
         DO n = 1, n_w 
            fl = fl + f_loc(n) * ww(n,l) 
         END DO

         gdf = 0
         DO k = 1, k_d
            gk = 0; dfk =0
            DO n = 1, n_w 
               gk = gk + g_loc(k,n) * ww(n,l) 
               dfk = dfk + f_loc(n) * dw_loc(k,n)
            END DO
            gdf = gdf + gk*dfk 
         ENDDO

         t = t + (fl**2 + gdf**2) * rj(l,m)

      ENDDO
   ENDDO

   t = SQRT(t)

END SUBROUTINE ns_adv

!------------------------------------------------------------------------------

SUBROUTINE ns_01 (mesh, gg, ff,  t)
!====================================

!  SQRT << (g.(Df))^2 >>   ===>   t

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
   REAL(KIND=8),                 INTENT(OUT) :: t

   INTEGER ::  m, l, k
   REAL(KIND=8) :: s

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   t = 0

   DO m = 1, me
      DO l = 1, l_G

         s = 0
         DO k = 1, k_d
            s = s + (SUM(gg(k, jj(:,m)) * ww(:,l))   &
                  * SUM(ff(jj(:,m)) * dw(k,:,l,m)))**2
         ENDDO

         t = t + s * rj(l,m)

      ENDDO
   ENDDO
   t = SQRT(t)

END SUBROUTINE ns_01

SUBROUTINE ts_01 (mesh, gg, ff,  t)
!====================================

!  << g.(Df) >>   ===>   t

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
   REAL(KIND=8),                 INTENT(OUT) :: t

   INTEGER ::  m, l, k
   REAL(KIND=8) :: s

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   t = 0

   DO m = 1, me

      DO l = 1, l_G

         s = 0
         DO k = 1, k_d
            s = s + SUM(gg(k, jj(:,m)) * ww(:,l))   &
                  * SUM(ff(jj(:,m)) * dw(k,:,l,m))
         ENDDO

         t = t + s * rj(l,m)

      ENDDO
   ENDDO

END SUBROUTINE ts_01

!------------------------------------------------------------------------------

SUBROUTINE tb_01 (mesh, gg, ff,  t)
!====================================

!  << (kxg).(Df) >>   ===>   t    (two dimensions only)

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
   REAL(KIND=8),                 INTENT(OUT) :: t

   INTEGER ::  m, l
   REAL(KIND=8) :: s

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   t = 0

   DO m = 1, me

      DO l = 1, l_G

         s =   SUM(gg(2, jj(:,m)) * ww(:,l)) * SUM(ff(jj(:,m)) * dw(1,:,l,m)) &
             - SUM(gg(1, jj(:,m)) * ww(:,l)) * SUM(ff(jj(:,m)) * dw(2,:,l,m))

         t = t + s * rj(l,m)

      ENDDO
   ENDDO

END SUBROUTINE tb_01

!------------------------------------------------------------------------------

SUBROUTINE tb_10 (mesh, gg, ff,  t)
!====================================

!  << (k.Dxg) f >>   ===>   t    (two dimensions only)

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
   REAL(KIND=8),                 INTENT(OUT) :: t

   INTEGER ::  m, l
   REAL(KIND=8) :: s

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   t = 0

   DO m = 1, me

      DO l = 1, l_G

         s =   SUM(gg(2, jj(:,m)) * dw(1,:,l,m))  &
             - SUM(gg(1, jj(:,m)) * dw(2,:,l,m))

         t = t + s * SUM(ff(jj(:,m)) * ww(:,l)) * rj(l,m)

      ENDDO
   ENDDO

END SUBROUTINE tb_10

!------------------------------------------------------------------------------

SUBROUTINE tv_11 (mesh, g1, g2,  t)
!====================================

!  << (D.g1) (D.g2) >>   ===>   t

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: g1, g2
   REAL(KIND=8),                 INTENT(OUT) :: t

   INTEGER ::  m, l
   REAL(KIND=8) :: s

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   t = 0

   DO m = 1, me

      DO l = 1, l_G

         s = SUM(g1(:, jj(:,m)) * dw(:,:,l,m)) * SUM(g2(:, jj(:,m)) * dw(:,:,l,m))

         t = t + s * rj(l,m)

      ENDDO
   ENDDO

END SUBROUTINE tv_11

!------------------------------------------------------------------------------

SUBROUTINE tc_11 (mesh, g1, g2,  t)
!====================================

!  << k.(Dxg1) k.(Dxg2) >>   ===>   t    (two dimensions only)

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: g1, g2
   REAL(KIND=8),                 INTENT(OUT) :: t

   INTEGER ::  m, l
   REAL(KIND=8) :: s

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   t = 0

   DO m = 1, me
      DO l = 1, l_G

         s = (   SUM(g1(2, jj(:,m)) * dw(1,:,l,m))    &
               - SUM(g1(1, jj(:,m)) * dw(2,:,l,m)) )  &
           * (   SUM(g2(2, jj(:,m)) * dw(1,:,l,m))    &
               - SUM(g2(1, jj(:,m)) * dw(2,:,l,m)) )

         t = t + s * rj(l,m)

      ENDDO
   ENDDO

END SUBROUTINE tc_11

!------------------------------------------------------------------------------

SUBROUTINE tb_Jac (mesh, pp, uu, ff,  t)
!=========================================

!  < J(u,p) f >   ===>   t    (two dimensions only)

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: pp, uu, ff
   REAL(KIND=8),                 INTENT(OUT) :: t

   INTEGER ::  m, l
   REAL(KIND=8) :: s

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   t = 0

   DO m = 1, me
      DO l = 1, l_G

         s = SUM(uu(jj(:,m)) * dw(1,:,l,m)) * SUM(pp(jj(:,m)) * dw(2,:,l,m)) &
           - SUM(uu(jj(:,m)) * dw(2,:,l,m)) * SUM(pp(jj(:,m)) * dw(1,:,l,m))

         t = t + s * SUM(ff(jj(:,m)) * ww(:,l)) * rj(l,m)

      ENDDO
   ENDDO

END SUBROUTINE tb_Jac

!------------------------------------------------------------------------------

SUBROUTINE tb_101 (mesh, pp, uu, ff,  t)
!=========================================

!  < u k x (Dp) . Df >   ===>   t    (two dimensions only)

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: pp, uu, ff
   REAL(KIND=8),                 INTENT(OUT) :: t

   INTEGER :: m, l
   REAL(KIND=8) :: s

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   t = 0

   DO m = 1, me
      DO l = 1, l_G

         s = - SUM(pp(jj(:,m)) * dw(2,:,l,m)) * SUM(ff(jj(:,m)) * dw(1,:,l,m)) &
             + SUM(pp(jj(:,m)) * dw(1,:,l,m)) * SUM(ff(jj(:,m)) * dw(2,:,l,m))

         t = t + SUM(uu(jj(:,m)) * ww(:,l)) * s * rj(l,m)

      ENDDO
   ENDDO

END SUBROUTINE tb_101

!------------------------------------------------------------------------------

SUBROUTINE tb_Jac_s (mesh, pp, uu, ff,  t)
!============================================

!  < u (T.Dp) f >_s   ===>   t    (two dimensions only)

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: pp, uu, ff
   REAL(KIND=8),                 INTENT(OUT) :: t

   INTEGER ::  ms, ls
   REAL(KIND=8) :: dp, u, dpu

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: js
   INTEGER,                      POINTER       :: mes

   CALL gauss(mesh)
   js => mesh%jjs
   mes => mesh%mes

   t = 0

   DO ms = 1, mes
      DO ls = 1, l_Gs

         dp = SUM(pp(js(:,ms)) * dwps(1,:,ls,ms))
         u  = SUM(uu(js(:,ms)) * wws(:,ls))
         dpu = dp * u

         t = t + dpu * SUM(ff(js(:,ms)) * wws(:,ls))

      ENDDO
   ENDDO

END SUBROUTINE tb_Jac_s

!------------------------------------------------------------------------------

SUBROUTINE tb_11s_s (mesh, pp, ff,  t)
!========================================

!  < (T.D_s f) |D_s p|^2 >_s   ===>   t    (two dimensions only)

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: pp, ff
   REAL(KIND=8),                 INTENT(OUT) :: t

   INTEGER ::  ms, ls
   REAL(KIND=8) :: dp

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: js
   INTEGER,                      POINTER       :: mes

   CALL gauss(mesh)
   js => mesh%jjs
   mes => mesh%mes

   t = 0

   DO ms = 1, mes
      DO ls = 1, l_Gs

         dp = SUM(pp(js(:,ms)) * dws(1,:,ls,ms))

         t = t + dp * dp * SUM(ff(js(:,ms)) * dwps(1,:,ls,ms))

      ENDDO
   ENDDO

END SUBROUTINE tb_11s_s

!------------------------------------------------------------------------------

SUBROUTINE ts_00_s (mesh, fs, ff, t)
!==========================================

!  < fs f >_s   ===>   t

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: fs, ff
   REAL(KIND=8),                 INTENT(OUT) :: t

   INTEGER ::  ms, ls

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: js, is
   INTEGER,                      POINTER       :: mes

   CALL gauss(mesh)
   js => mesh%jjs
   is => mesh%iis
   mes => mesh%mes

   t = 0

   DO ms = 1, mes
      DO ls = 1, l_Gs

         t = t + SUM(fs(is(:,ms)) * wws(:,ls))   &
               * SUM(ff(js(:,ms)) * wws(:,ls)) * rjs(ls,ms)

      ENDDO
   ENDDO

END SUBROUTINE ts_00_s

SUBROUTINE ns_anal (mesh, ff, f_anal,  t)
!===============================================
   USE Gauss_points
   IMPLICIT NONE
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: f_anal
   REAL(KIND=8),                 INTENT(OUT) :: t
   INTEGER ::  m, l, n, index
   REAL(KIND=8) :: fl 
   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   t = 0

   DO m = 1, me
      index = (m-1)*l_G
      DO l = 1, l_G; index = index + 1
         fl = 0
         DO n = 1, n_w
            fl = fl + ff(jj(n,m)) * ww(n,l)
         END DO
         t = t + (fl - f_anal(index)) * rj(l,m)
      ENDDO
   ENDDO

   t = ABS(t)
END SUBROUTINE ns_anal

SUBROUTINE ns_anal_0 (mesh, ff, f_anal,  t)
!===============================================

!  sqrt(< f^2 >)   ===>   t

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: f_anal
   REAL(KIND=8),                 INTENT(OUT) :: t

   INTEGER ::  m, l, n, index
   REAL(KIND=8) :: fl 

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   t = 0

   DO m = 1, me
      index = (m-1)*l_G
      DO l = 1, l_G; index = index + 1
         fl = 0
         DO n = 1, n_w
            fl = fl + ff(jj(n,m)) * ww(n,l)
         END DO
         t = t + (fl - f_anal(index))**2 * rj(l,m)
      ENDDO
   ENDDO

   t = SQRT(t)


END SUBROUTINE ns_anal_0

SUBROUTINE ns_anal_infty (mesh, ff, f_anal,  t)
!===============================================

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: f_anal
   REAL(KIND=8),                 INTENT(OUT) :: t

   INTEGER ::  m, l, n, index
   REAL(KIND=8) :: fl 

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   t = 0

   DO m = 1, me
      index = (m-1)*l_G
      DO l = 1, l_G; index = index + 1
         fl = 0
         DO n = 1, n_w
            fl = fl + ff(jj(n,m)) * ww(n,l)
         END DO
         t = MAX(t,ABS(fl - f_anal(index)))
      ENDDO
   ENDDO


END SUBROUTINE ns_anal_infty

SUBROUTINE nv_anal_0 (mesh, ff, f_anal,  t)
!===============================================

!  sqrt(< f^2 >)   ===>   t

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: ff
   REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: f_anal
   REAL(KIND=8),                   INTENT(OUT) :: t

   INTEGER ::  m, l, n, index, k
   REAL(KIND=8) :: fl 

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   t = 0

   DO m = 1, me
      index = (m-1)*l_G
      DO l = 1, l_G; index = index + 1
         DO k =1 , SIZE(ff,1)
            fl = 0
            DO n = 1, n_w
               fl = fl + ff(k,jj(n,m)) * ww(n,l)
            END DO
            t = t + (fl - f_anal(k,index))**2 * rj(l,m)
         END DO
      ENDDO
   ENDDO

   t = SQRT(t)


END SUBROUTINE nv_anal_0

SUBROUTINE nv_anal_l1 (mesh, ff, f_anal,  t)
!===============================================

!  (< |f| >)   ===>   t

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: ff
   REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: f_anal
   REAL(KIND=8),                   INTENT(OUT) :: t

   INTEGER ::  m, l, n, index, k
   REAL(KIND=8) :: fl 

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   t = 0
 
   DO m = 1, me
      index = (m-1)*l_G
      DO l = 1, l_G; index = index + 1
         DO k =1 , SIZE(ff,1)
            fl = 0
            DO n = 1, n_w
               fl = fl + ff(k,jj(n,m)) * ww(n,l)
            END DO
            t = t + ABS(fl - f_anal(k,index)) * rj(l,m)
         END DO
      ENDDO
   ENDDO

END SUBROUTINE nv_anal_l1


SUBROUTINE ns_anal_l1 (mesh, ff, f_anal,  t)
!===============================================

!  < |f| >)   ===>   t

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: f_anal
   REAL(KIND=8),                 INTENT(OUT) :: t

   INTEGER ::  m, l, n, index
   REAL(KIND=8) :: fl

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   t = 0

   DO m = 1, me
      index = (m-1)*l_G
      DO l = 1, l_G; index = index + 1
         fl = 0
         DO n = 1, n_w
            fl = fl + ff(jj(n,m)) * ww(n,l)
         END DO
         t = t + ABS(fl - f_anal(index)) * rj(l,m)
      ENDDO
   ENDDO

END SUBROUTINE ns_anal_l1

SUBROUTINE nb_anal_1 (mesh, f_anal, gg,  t)
!===============================

!  sqrt(< (D x g . k) (D x g . k) >)   ===>   t     ( 2d only )

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: f_anal
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8),                 INTENT(OUT) :: t

   INTEGER ::  m, l, index
   REAL(KIND=8) :: s

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   t = 0

   DO m = 1, me
      index = (m-1)*l_G
      DO l = 1, l_G; index = index + 1

         s = (f_anal(index) &
               -SUM(gg(2,jj(:,m)) * dw(1,:,l,m)   &
                  - gg(1,jj(:,m)) * dw(2,:,l,m)))**2

         t = t + s * rj(l,m)

      ENDDO
   ENDDO

   t = SQRT(t)

END SUBROUTINE nb_anal_1

SUBROUTINE ns_anal_1 (mesh, ff, f_anal, gdf_anal,  t)
!====================================

!  SQRT << (ff-f_anal)^2 + ((Df-gdf_anal))^2 >>   ===>   t

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff, f_anal
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gdf_anal
   REAL(KIND=8),                 INTENT(OUT) :: t

   REAL(KIND=8), DIMENSION(k_d,n_w)          :: dw_loc
   REAL(KIND=8), DIMENSION(n_w)              :: f_loc 
   INTEGER ::  m, l, k, n, index
   REAL(KIND=8) :: fl, dfk, gdf 

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   t = 0

   index = 0
   DO m = 1, me

      f_loc = ff(jj(:,m))  

      DO l = 1, l_G; index = index + 1

         dw_loc = dw(:,:,l,m)


         fl = 0
         DO n = 1, n_w 
            fl = fl + f_loc(n) * ww(n,l) 
         END DO

         gdf = 0
         DO k = 1, k_d
            dfk =0
            DO n = 1, n_w 
               dfk = dfk + f_loc(n) * dw_loc(k,n)
            END DO
            gdf = gdf + (dfk - gdf_anal(k,index))**2 
         ENDDO

         t = t + ((fl-f_anal(index))**2 + gdf) * rj(l,m)

      ENDDO
   ENDDO

   t = SQRT(t)

END SUBROUTINE ns_anal_1

SUBROUTINE ns_anal_w11 (mesh, ff, f_anal, gdf_anal,  t)
!====================================

!   << |ff-f_anal| + ||Df-gdf_anal|| >>   ===>   t

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff, f_anal
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gdf_anal
   REAL(KIND=8),                 INTENT(OUT) :: t

   REAL(KIND=8), DIMENSION(k_d,n_w)          :: dw_loc
   REAL(KIND=8), DIMENSION(n_w)              :: f_loc 
   INTEGER ::  m, l, k, n, index
   REAL(KIND=8) :: fl, dfk, gdf 

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   t = 0

   index = 0
   DO m = 1, me

      f_loc = ff(jj(:,m))  

      DO l = 1, l_G; index = index + 1

         dw_loc = dw(:,:,l,m)


         fl = 0
         DO n = 1, n_w 
            fl = fl + f_loc(n) * ww(n,l) 
         END DO

         gdf = 0
         DO k = 1, k_d
            dfk =0
            DO n = 1, n_w 
               dfk = dfk + f_loc(n) * dw_loc(k,n)
            END DO
            gdf = gdf + ABS(dfk - gdf_anal(k,index))
         ENDDO

         t = t + (ABS(fl-f_anal(index)) + gdf) * rj(l,m)

      ENDDO
   ENDDO

END SUBROUTINE ns_anal_w11

!------------------------------------------------------------------------------

SUBROUTINE ns_anal_adv (mesh, gg, ff, f_anal, gdf_anal,  t)
!====================================

!  SQRT << (ff-f_anal)^2 + (g.(Df-gdf_anal))^2 >>   ===>   t

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff, f_anal
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gdf_anal
   REAL(KIND=8),                 INTENT(OUT) :: t

   REAL(KIND=8), DIMENSION(k_d,n_w)          :: dw_loc, g_loc 
   REAL(KIND=8), DIMENSION(n_w)              :: f_loc 
   INTEGER ::  m, l, k, n, index
   REAL(KIND=8) :: fl, gk, dfk, gdf 

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   t = 0

   index = 0
   DO m = 1, me

      g_loc = gg(:, jj(:,m))  
      f_loc = ff(jj(:,m))  

      DO l = 1, l_G; index = index + 1

         dw_loc = dw(:,:,l,m)


         fl = 0
         DO n = 1, n_w 
            fl = fl + f_loc(n) * ww(n,l) 
         END DO

         gdf = 0
         DO k = 1, k_d
            gk = 0; dfk =0
            DO n = 1, n_w 
               gk = gk + g_loc(k,n) * ww(n,l) 
               dfk = dfk + f_loc(n) * dw_loc(k,n)
            END DO
            gdf = gdf + gk*(dfk - gdf_anal(k,index)) 
         ENDDO

         t = t + ((fl - f_anal(index))**2 + gdf**2) * rj(l,m)

      ENDDO
   ENDDO

   t = SQRT(t)

END SUBROUTINE ns_anal_adv

SUBROUTINE ns_anal_adv_l1 (mesh, gg, ff, f_anal, gdf_anal,  t)
!====================================

!  SQRT << (ff-f_anal)^2 + (g.(Df-gdf_anal))^2 >>   ===>   t

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff, f_anal
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gdf_anal
   REAL(KIND=8),                 INTENT(OUT) :: t

   REAL(KIND=8), DIMENSION(k_d,n_w)          :: dw_loc, g_loc 
   REAL(KIND=8), DIMENSION(n_w)              :: f_loc 
   INTEGER ::  m, l, k, n, index
   REAL(KIND=8) :: fl, gk, dfk, gdf 

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   t = 0

   index = 0
   DO m = 1, me

      g_loc = gg(:, jj(:,m))  
      f_loc = ff(jj(:,m))  

      DO l = 1, l_G; index = index + 1

         dw_loc = dw(:,:,l,m)


         fl = 0
         DO n = 1, n_w 
            fl = fl + f_loc(n) * ww(n,l) 
         END DO

         gdf = 0
         DO k = 1, k_d
            gk = 0; dfk =0
            DO n = 1, n_w 
               gk = gk + g_loc(k,n) * ww(n,l) 
               dfk = dfk + f_loc(n) * dw_loc(k,n)
            END DO
            gdf = gdf + gk*(dfk - gdf_anal(k,index)) 
         ENDDO

         t = t + (ABS(fl - f_anal(index)) + ABS(gdf)) * rj(l,m)

      ENDDO
   ENDDO


END SUBROUTINE ns_anal_adv_l1

!------------------------------------------------------------------------------

SUBROUTINE ns_test (mesh, gg, ff,  t)
!====================================

!  SQRT << (g.Df,f >>   ===>   t

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
   REAL(KIND=8),                 INTENT(OUT) :: t

   REAL(KIND=8), DIMENSION(k_d,n_w)          :: dw_loc, g_loc 
   REAL(KIND=8), DIMENSION(n_w)              :: f_loc 
   INTEGER ::  m, l, k, n
   REAL(KIND=8) :: fl, gk, dfk, gdf 

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   t = 0

   DO m = 1, me

      g_loc = gg(:, jj(:,m))  
      f_loc = ff(jj(:,m))  

      DO l = 1, l_G;  !index = index + 1

         dw_loc = dw(:,:,l,m)


         fl = 0
         DO n = 1, n_w 
            fl = fl + f_loc(n) * ww(n,l) 
         END DO

         gdf = 0
         DO k = 1, k_d
            gk = 0; dfk =0
            DO n = 1, n_w 
               gk = gk + g_loc(k,n) * ww(n,l) 
               dfk = dfk + f_loc(n) * dw_loc(k,n)
            END DO
            gdf = gdf + gk*dfk 
         ENDDO

         t = t +  gdf*fl* rj(l,m)

      ENDDO
   ENDDO

END SUBROUTINE ns_test

SUBROUTINE n_bloc_0 (mesh, gg, nb_bloc, bloc_size, t)
!===============================

!  sqrt(<< g.g >>)   ===>   t

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: gg
   INTEGER,                    INTENT(IN)  :: nb_bloc, bloc_size
   REAL(KIND=8),                 INTENT(OUT) :: t

   INTEGER ::  m, l, k, n
   REAL(KIND=8) :: s, gl

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   t = 0

   DO m = 1, me
      DO l = 1, l_G

         s = 0
         DO k = 1, nb_bloc 
            gl = 0
            DO n = 1, n_w
               gl = gl + gg((k-1)*bloc_size+jj(n,m)) * ww(n,l)
            END DO
            s = s + gl**2
         ENDDO

         t = t + s * rj(l,m)

      ENDDO
   ENDDO

   t = SQRT(t)

END SUBROUTINE n_bloc_0

SUBROUTINE ns_l1 (mesh, ff, t)
!===============================================

!  < |f| >)   ===>   t

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
   REAL(KIND=8),                 INTENT(OUT) :: t

   INTEGER ::  m, l, n
   REAL(KIND=8) :: fl

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me
   t = 0

   DO m = 1, me
      DO l = 1, l_G
         fl = 0
         DO n = 1, n_w
            fl = fl + ff(jj(n,m)) * ww(n,l)
         END DO
         t = t + ABS(fl) * rj(l,m)
      ENDDO
   ENDDO

END SUBROUTINE ns_l1


!------------------------------------------------------------------------------
SUBROUTINE ns_01_hybrid (uu_mesh, pp_mesh, ff, gg,  u0_c)
!================================================

!  < f, D.g >   ===>   u0_c

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff 
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8),                 INTENT(OUT) :: u0_c

   TYPE(mesh_type), TARGET                     :: uu_mesh, pp_mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj, jj_c
   INTEGER,      POINTER       :: me

   REAL(KIND=8), DIMENSION(pp_mesh%gauss%n_w, uu_mesh%gauss%l_G) ::  w_c

   INTEGER :: m, l, n, k
   REAL(KIND=8) :: dgl
   REAL(KIND=8), DIMENSION(uu_mesh%gauss%k_d,uu_mesh%gauss%n_w) :: dwkn, ggkn 
   REAL(KIND=8), DIMENSION(uu_mesh%gauss%n_w) :: u0_cn

   CALL gauss(uu_mesh)
   jj => uu_mesh%jj
   jj_c => pp_mesh%jj
   me => uu_mesh%me


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

      u0_cn = 0

      DO l = 1, l_G
         dwkn(:,:) = dw(:,:,l,m)

         dgl = 0
         DO n = 1, n_w
            DO k = 1, k_d
               dgl = dgl + ggkn(k,n) * dwkn(k,n)
            END DO
         END DO

! FINIR d'optimiser la boucle.

         u0_c = u0_c + dgl * SUM(ff(jj_c(:,m))*w_c(:,l)) * rj(l,m)

      ENDDO

   ENDDO

END SUBROUTINE ns_01_hybrid

SUBROUTINE eval_h (mesh, hh)
!===============================

!  h  ===>   t

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:),   INTENT(OUT)  :: hh

   INTEGER :: m
   REAL(KIND=8) :: h , inv_dim

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me
   inv_dim = 1.d0/SIZE(mesh%rr,1)

   hh  = 0

   DO m = 1, me
      h = SUM(rj(:,m))**inv_dim
      hh(jj(:,m)) =  hh(jj(:,m)) + h 
   ENDDO

END SUBROUTINE eval_h 

END MODULE fem_tn
