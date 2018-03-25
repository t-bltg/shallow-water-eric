MODULE Gauss_points

   PRIVATE

   INTEGER, PARAMETER, PUBLIC :: k_d = 3,  n_w  = 4,  l_G  = 4,   &
                                           n_ws = 3,  l_Gs = 3

   REAL(KIND=8), DIMENSION(n_w,  l_G),               PUBLIC :: ww
   REAL(KIND=8), DIMENSION(n_ws, l_Gs),              PUBLIC :: wws

   REAL(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE, PUBLIC :: dw
   REAL(KIND=8), DIMENSION(:,    :, :), ALLOCATABLE, PUBLIC :: rnorms
   REAL(KIND=8), DIMENSION(:, :),       ALLOCATABLE, PUBLIC :: rj
   REAL(KIND=8), DIMENSION(:, :),       ALLOCATABLE, PUBLIC :: rjs
   REAL(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE, PUBLIC :: dwps !special!
   REAL(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE, PUBLIC :: dws  !SPECIAL!

!  REAL(KIND=8), DIMENSION(k_d,  n_w,  l_G,   me),   PUBLIC :: dw
!  REAL(KIND=8), DIMENSION(k_d,        l_Gs,  mes),  PUBLIC :: rnorms
!  REAL(KIND=8), DIMENSION(l_G,   me),               PUBLIC :: rj
!  REAL(KIND=8), DIMENSION(l_Gs,  mes),              PUBLIC :: rjs

   PUBLIC Gauss_gen

CONTAINS

SUBROUTINE Gauss_gen(np, me, nps, mes, jj, js, rr)

   IMPLICIT NONE

   INTEGER, INTENT(IN) :: np, me, nps, mes

   INTEGER,      DIMENSION(n_w,  me),  INTENT(IN) :: jj
   INTEGER,      DIMENSION(n_ws, mes), INTENT(IN) :: js
   REAL(KIND=8), DIMENSION(k_d,  np),  INTENT(IN) :: rr

   REAL(KIND=8), DIMENSION(k_d,   n_w,  l_G)  :: dd
   REAL(KIND=8), DIMENSION(k_d-1, n_ws, l_Gs) :: dds
   REAL(KIND=8), DIMENSION(l_G)               :: pp
   REAL(KIND=8), DIMENSION(l_Gs)              :: pps

   REAL(KIND=8), DIMENSION(n_w)        :: r
   REAL(KIND=8), DIMENSION(n_ws)       :: rs
   REAL(KIND=8), DIMENSION(k_d,   k_d) :: dr
   REAL(KIND=8), DIMENSION(k_d-1, k_d) :: drs

   REAL(KIND=8), DIMENSION(k_d, k_d) :: mnr
   REAL(KIND=8), DIMENSION(k_d)      :: mnrs
   REAL(KIND=8) :: rjac, rjacs
   INTEGER      :: m, l, k, k1, k2, h, h1, h2, n,  ms, ls

   IF(ALLOCATED(dw)) THEN
      DEALLOCATE(dw,rnorms,rj,rjs,dwps)
   END IF

   ALLOCATE(dw(k_d, n_w, l_G,  me ))
   ALLOCATE(rnorms(k_d,  l_Gs, mes))
   ALLOCATE( rj(l_G,  me ))
   ALLOCATE(rjs(l_Gs, mes))
   ALLOCATE(dwps(k_d-1,  n_ws,  l_Gs,  mes))
   ALLOCATE(dws(k_d-1,  n_ws,  l_Gs,  mes))

!  evaluate and store the values of derivatives and of the
!  jacobian determinant at Gauss points of all volume elements

!  volume elements

   CALL element_3d(ww, dd, pp)

   DO m = 1, me

      DO l = 1, l_G

         DO k = 1, k_d
            r = rr(k, jj(:,m))
            DO k1 = 1, k_d
               dr(k, k1) = SUM(r * dd(k1,:,l))
            ENDDO
         ENDDO

         DO k = 1, k_d;  k1 = MODULO(k,k_d) + 1;  k2 = MODULO(k+1,k_d) + 1
            DO h = 1, k_d;  h1 = MODULO(h,k_d) + 1;  h2 = MODULO(h+1,k_d) + 1
               mnr(k,h) = dr(k1,h1) * dr(k2,h2)  -  dr(k1,h2) * dr(k2,h1)
            ENDDO
         ENDDO

         rjac = SUM(dr(1,:) * mnr(1,:))

!   rjac = dr(1,1)*dr(2,2)*dr(3,3)+dr(1,2)*dr(2,3)*dr(3,1)+dr(1,3)*dr(2,1)*dr(3,2) &
!        - dr(3,1)*dr(2,2)*dr(1,3)-dr(3,2)*dr(2,3)*dr(1,1)+dr(3,3)*dr(2,1)*dr(1,2) 

         DO n = 1, n_w
            DO k = 1, k_d
               dw(k, n, l, m) = SUM(mnr(k,:) * dd(:,n,l))/rjac
            ENDDO
         ENDDO

         rj(l, m) = rjac * pp(l)

      ENDDO

   ENDDO

!  surface elements

   CALL element_2d(wws, dds, pps)

   DO ms = 1, mes

      DO ls = 1, l_Gs

         DO k = 1, k_d
            rs = rr(k, js(:,ms))
            DO h = 1, 2  ! = k_d - 1
               drs(h, k) = SUM(rs * dds(h,:,ls))
            ENDDO
         ENDDO

         DO k = 1, k_d;  k1 = MODULO(k,k_d) + 1;  k2 = MODULO(k+1,k_d) + 1
            mnrs(k) = drs(1,k1) * drs(2,k2)  -  drs(2,k1) * drs(1,k2)
         ENDDO

         rjacs = SQRT(SUM(mnrs**2))

         rnorms(:, ls, ms) = mnrs/rjacs   ! outward normal

         rjs(ls, ms) = rjacs * pps(ls)

      ENDDO

   ENDDO


   DO ms = 1, mes  ! necessary only for evaluating gradient
                   ! tangential to the surface (ex COMMON Gauss_tan)

      DO ls = 1, l_Gs

          dwps(:, :, ls, ms) = dds(:, :, ls) * pps(ls)
           dws(:, :, ls, ms) = dds(:, :, ls)

      ENDDO

   ENDDO

   print*, 'end of gen_Gauss_3d_p1'

   CONTAINS

   SUBROUTINE element_3d (w, d, p)

!     tetraedral element with linear interpolation
!     and four Gauss integration points

!        w(n_w, l_G) : values of shape functions at Gauss points
!     d(3, n_w, l_G) : derivatives values of shape functions at Gauss points
!             p(l_G) : weight for Gaussian quadrature at Gauss points

      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(   n_w, l_G), INTENT(OUT) :: w
      REAL(KIND=8), DIMENSION(3, n_w, l_G), INTENT(OUT) :: d
      REAL(KIND=8), DIMENSION(l_G),         INTENT(OUT) :: p

      REAL(KIND=8), DIMENSION(l_G) :: xx, yy, zz
      INTEGER :: j

      REAL(KIND=8) :: zero = 0,  one = 1.,  five = 5.
      REAL(KIND=8) :: a, b

      REAL(KIND=8) :: f1, f2, f3, f4, x, y, z
      f1(x, y, z) = one - x - y - z
      f2(x, y, z) = x
      f3(x, y, z) = y
      f4(x, y, z) = z

      a = (five  -  SQRT(five))/20
      b = (five + 3.*SQRT(five))/20

      xx    = a;      yy = a;      zz = a
      xx(3) = b;   yy(2) = b;   zz(1) = b

      DO j = 1, l_G

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

!        w(n_ws, l_Gs) : values of shape functions at Gauss points
!     d(2, n_ws, l_Gs) : derivatives values of shape functions at Gauss points
!              p(l_Gs) : weight for Gaussian quadrature at Gauss points

      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(   n_ws, l_Gs), INTENT(OUT) :: w
      REAL(KIND=8), DIMENSION(2, n_ws, l_Gs), INTENT(OUT) :: d
      REAL(KIND=8), DIMENSION(l_Gs),          INTENT(OUT) :: p

      REAL(KIND=8), DIMENSION(l_Gs) :: xx, yy
      INTEGER :: j

      REAL(KIND=8) :: zero = 0,  one = 1,  two = 2,  three = 3,  six = 6

      REAL(KIND=8) :: f1, f2, f3, x, y
      f1(x, y) = one - x - y
      f2(x, y) = x
      f3(x, y) = y

      xx(1) = one/six;  xx(2) = two/three;  xx(3) = one/six
      yy(1) = one/six;  yy(2) = one/six;    yy(3) = two/three

      DO j = 1, l_Gs

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


END SUBROUTINE Gauss_gen

END MODULE Gauss_points
