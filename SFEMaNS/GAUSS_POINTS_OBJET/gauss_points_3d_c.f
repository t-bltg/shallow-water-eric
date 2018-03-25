MODULE Gauss_points_c

   PRIVATE

   INTEGER, PARAMETER, PUBLIC :: k_d_c = 3,  n_w_c  = 4,  l_G_c  = 4,   &
                                           n_ws_c = 3,  l_Gs_c = 3

   REAL(KIND=8), DIMENSION(n_w_c,  l_G_c),            PUBLIC :: ww_c
   REAL(KIND=8), DIMENSION(n_ws_c, l_Gs_c),           PUBLIC :: wws_c

   REAL(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE,  PUBLIC :: dw_c
   REAL(KIND=8), DIMENSION(:,    :, :), ALLOCATABLE,  PUBLIC :: rnorms_c
   REAL(KIND=8), DIMENSION(:, :),       ALLOCATABLE,  PUBLIC :: rj_c
   REAL(KIND=8), DIMENSION(:, :),       ALLOCATABLE,  PUBLIC :: rjs_c
   REAL(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE,  PUBLIC :: dwps_c !special!
   REAL(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE,  PUBLIC :: dws_c  !SPECIAL!

!  REAL(KIND=8), DIMENSION(k_d_c,  n_w_c,  l_G_c,  me), PUBLIC :: dw_c
!  REAL(KIND=8), DIMENSION(k_d_c,        l_Gs_c,  mes), PUBLIC :: rnorms_c
!  REAL(KIND=8), DIMENSION(l_G_c,   me),              PUBLIC :: rj_c
!  REAL(KIND=8), DIMENSION(l_Gs_c,  mes),             PUBLIC :: rjs_c

   PUBLIC Gauss_gen_c

CONTAINS

SUBROUTINE Gauss_gen_c(np, me, nps, mes, jj, js, rr)

   IMPLICIT NONE

   INTEGER, INTENT(IN) :: np, me, nps, mes

   INTEGER,      DIMENSION(n_w_c,  me),  INTENT(IN) :: jj
   INTEGER,      DIMENSION(n_ws_c, mes), INTENT(IN) :: js
   REAL(KIND=8), DIMENSION(k_d_c,  np),    INTENT(IN) :: rr

   REAL(KIND=8), DIMENSION(k_d_c,   n_w_c,  l_G_c)  :: dd
   REAL(KIND=8), DIMENSION(k_d_c-1, n_ws_c, l_Gs_c) :: dds
   REAL(KIND=8), DIMENSION(l_G_c)                 :: pp
   REAL(KIND=8), DIMENSION(l_Gs_c)                :: pps

   REAL(KIND=8), DIMENSION(n_w_c)        :: r
   REAL(KIND=8), DIMENSION(n_ws_c)       :: rs
   REAL(KIND=8), DIMENSION(k_d_c,   k_d_c)   :: dr
   REAL(KIND=8), DIMENSION(k_d_c-1, k_d_c)   :: drs

   REAL(KIND=8), DIMENSION(k_d_c, k_d_c) :: mnr
   REAL(KIND=8), DIMENSION(k_d_c)      :: mnrs
   REAL(KIND=8) :: rjac_c, rjacs_c
   INTEGER      :: m, l, k, k1, k2, h, h1, h2, n,  ms, ls

   IF(ALLOCATED(dw_c)) THEN
      DEALLOCATE(dw_c,rnorms_c,rj_c,rjs_c,dwps_c)
   END IF

   ALLOCATE(dw_c(k_d_c, n_w_c, l_G_c,  me ))
   ALLOCATE(rnorms_c(k_d_c,  l_Gs_c, mes))
   ALLOCATE( rj_c(l_G_c,  me ))
   ALLOCATE(rjs_c(l_Gs_c, mes))
   ALLOCATE(dwps_c(k_d_c-1,  n_ws_c,  l_Gs_c,  mes))
   ALLOCATE(dws_c(k_d_c-1,  n_ws_c,  l_Gs_c,  mes))

!  evaluate and store the values of derivatives and of the
!  jacobian determinant at Gauss points of all volume elements

!  volume elements

   CALL element_3d(ww_c, dd, pp)

   DO m = 1, me

      DO l = 1, l_G_c

         DO k = 1, k_d_c
            r = rr(k, jj(:,m))
            DO k1 = 1, k_d_c
               dr(k, k1) = SUM(r * dd(k1,:,l))
            ENDDO
         ENDDO

         DO k = 1, k_d_c;  k1 = MODULO(k,k_d_c) + 1;  k2 = MODULO(k+1,k_d_c) + 1
            DO h = 1, k_d_c;  h1 = MODULO(h,k_d_c) + 1;  h2 = MODULO(h+1,k_d_c) + 1
               mnr(k,h) = dr(k1,h1) * dr(k2,h2)  -  dr(k1,h2) * dr(k2,h1)
            ENDDO
         ENDDO

         rjac_c = SUM(dr(1,:) * mnr(1,:))

         DO n = 1, n_w_c
            DO k = 1, k_d_c
               dw_c(k, n, l, m) = SUM(mnr(k,:) * dd(:,n,l))/rjac_c
            ENDDO
         ENDDO

         rj_c(l, m) = rjac_c * pp(l)

      ENDDO

   ENDDO

!  surface elements

   CALL element_2d(wws_c, dds, pps)

   DO ms = 1, mes

      DO ls = 1, l_Gs_c

         DO k = 1, k_d_c
            rs = rr(k, js(:,ms))
            DO h = 1, 2  ! = k_d_c - 1
               drs(h, k) = SUM(rs * dds(h,:,ls))
            ENDDO
         ENDDO

         DO k = 1, k_d_c;  k1 = MODULO(k,k_d_c) + 1;  k2 = MODULO(k+1,k_d_c) + 1
            mnrs(k) = drs(1,k1) * drs(2,k2)  -  drs(2,k1) * drs(1,k2)
         ENDDO

         rjacs_c = SQRT(SUM(mnrs**2))

         rnorms_c(:, ls, ms) = mnrs/rjacs_c   ! outward normal

         rjs_c(ls, ms) = rjacs_c * pps(ls)

      ENDDO

   ENDDO


   DO ms = 1, mes  ! necessary only for evaluating gradient
                   ! tangential to the surface (ex COMMON Gauss_tan)

      DO ls = 1, l_Gs_c

          dwps_c(:, :, ls, ms) = dds(:, :, ls) * pps(ls)
           dws_c(:, :, ls, ms) = dds(:, :, ls)

      ENDDO

   ENDDO

   print*, 'end of gen_Gauss_3d_c_p1'

   CONTAINS

   SUBROUTINE element_3d (w, d, p)

!     tetraedral  element with linear interpolation
!     and four Gauss integration points

!        w(n_w_c, l_G_c) : values of shape functions at Gauss points
!     d(3, n_w_c, l_G_c) : derivatives values of shape functions at Gauss points
!               p(l_G_c) : weight for Gaussian quadrature at Gauss points

      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(   n_w_c, l_G_c), INTENT(OUT) :: w
      REAL(KIND=8), DIMENSION(3, n_w_c, l_G_c), INTENT(OUT) :: d
      REAL(KIND=8), DIMENSION(l_G_c),           INTENT(OUT) :: p

      REAL(KIND=8), DIMENSION(l_G_c) :: xx, yy, zz
      INTEGER :: j

      REAL(KIND=8) :: zero = 0,  one = 1,  five = 5
      REAL(KIND=8) :: a, b

      REAL(KIND=8) :: f1, f2, f3, f4, x, y, z
      f1(x, y, z) = one - x - y - z
      f2(x, y, z) = x
      f3(x, y, z) = y
      f4(x, y, z) = z

      a = (five  -  SQRT(five))/20
      b = (five + 3*SQRT(five))/20

      xx    = a;      yy = a;      zz = a
      xx(3) = b;   yy(2) = b;   zz(1) = b

      DO j = 1, l_G_c

            w(1, j) = f1(xx(j), yy(j), zz(j))
         d(1, 1, j) = - one
         d(2, 1, j) = - one
         d(3, 1, j) = - one

            w(2, j) = f2(xx(j), yy(j), zz(j))
         d(1, 2, j) = one
         d(2, 2, j) = zero
         d(3, 2, j) = zero

            w(3, j) = f3(xx(j), yy(j), zz(j))
         d(1, 3, j) = zero
         d(2, 3, j) = one
         d(3, 3, j) = zero

            w(4, j) = f4(xx(j), yy(j), zz(j))
         d(1, 4, j) = zero
         d(2, 4, j) = zero
         d(3, 4, j) = one

               p(j) = one/24

      ENDDO

   END SUBROUTINE element_3d

!------------------------------------------------------------------------------

   SUBROUTINE element_2d (w, d, p)

!     triangular element with linear interpolation
!     and three Gauss integration points

!        w(n_ws_c, l_Gs_c) : values of shape functions at Gauss points
!     d(2, n_ws_c, l_Gs_c) : derivatives values of shape functions at Gauss points
!                p(l_Gs_c) : weight for Gaussian quadrature at Gauss points

      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(   n_ws_c, l_Gs_c), INTENT(OUT) :: w
      REAL(KIND=8), DIMENSION(2, n_ws_c, l_Gs_c), INTENT(OUT) :: d
      REAL(KIND=8), DIMENSION(l_Gs_c),            INTENT(OUT) :: p

      REAL(KIND=8), DIMENSION(l_Gs_c) :: xx, yy
      INTEGER :: j

      REAL(KIND=8) :: zero = 0,  one = 1,  two = 2,  three = 3,  six = 6

      REAL(KIND=8) :: f1, f2, f3, x, y
      f1(x, y) = one - x - y
      f2(x, y) = x
      f3(x, y) = y

      xx(1) = one/six;  xx(2) = two/three;  xx(3) = one/six
      yy(1) = one/six;  yy(2) = one/six;    yy(3) = two/three

      DO j = 1, l_Gs_c

            w(1, j) = f1(xx(j), yy(j))
         d(1, 1, j) = - one
         d(2, 1, j) = - one

            w(2, j) = f2(xx(j), yy(j))
         d(1, 2, j) = one
         d(2, 2, j) = zero

            w(3, j) = f3(xx(j), yy(j))
         d(1, 3, j) = zero
         d(2, 3, j) = one

               p(j) = one/six

      ENDDO

   END SUBROUTINE element_2d

END SUBROUTINE Gauss_gen_c

END MODULE Gauss_points_c
