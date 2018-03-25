MODULE Gauss_points

   PRIVATE

   INTEGER, PARAMETER, PUBLIC :: k_d = 3,  n_w  = 10,  l_G  = 15,   &
                                           n_ws = 6,   l_Gs = 7

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

         DO n = 1, n_w
            DO k = 1, k_d
               dw(k, n, l, m) = SUM(mnr(k,:) * dd(:,n,l))/rjac
            ENDDO
         ENDDO

         rj(l, m) = rjac * pp(l)

      ENDDO

   ENDDO

!  surface elements

   CALL element_2d_p2(wws, dds, pps)

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

   print*, 'end of gen_Gauss_3d_p2'

   CONTAINS

   SUBROUTINE element_3d (w, d, p)

!     tetraedral element with quadratic interpolation
!     and fifteen Gauss integration points

!        w(n_w, l_G) : values of shape functions at Gauss points
!     d(3, n_w, l_G) : derivatives values of shape functions at Gauss points
!             p(l_G) : weight for Gaussian quadrature at Gauss points

      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(   n_w, l_G), INTENT(OUT) :: w
      REAL(KIND=8), DIMENSION(3, n_w, l_G), INTENT(OUT) :: d
      REAL(KIND=8), DIMENSION(l_G),         INTENT(OUT) :: p

      REAL(KIND=8), DIMENSION(l_G) :: xx, yy, zz
      INTEGER :: j

      REAL(KIND=8) :: zero = 0,  one = 1,  two = 2,  three = 3,  four = 4,  &
                      five = 5,  six = 6,                                   &
                      r, a,  s1, s2, t1, t2, b1, b2,  u, v, c,  vol, sq

      REAL(KIND=8) ::  f1,   f2,   f3,   f4,   f5,   f6,   f7,   f8,   f9,   f10,  & 
                      df1x, df2x, df3x, df4x, df5x, df6x, df7x, df8x, df9x, df10x, & 
                      df1y, df2y, df3y, df4y, df5y, df6y, df7y, df8y, df9y, df10y, & 
                      df1z, df2z, df3z, df4z, df5z, df6z, df7z, df8z, df9z, df10z, & 
                      x, y, z, t, tj

      f1(x, y, z, t) = t * (two*t - one)
      f2(x, y, z, t) = x * (two*x - one)
      f3(x, y, z, t) = y * (two*y - one) 
      f4(x, y, z, t) = z * (two*z - one)
      f5(x, y, z, t) = four * y * z 
      f6(x, y, z, t) = four * z * x
      f7(x, y, z, t) = four * x * y
      f8(x, y, z, t) = four * z * t
      f9(x, y, z, t) = four * y * t
     f10(x, y, z, t) = four * x * t 

    df1x(x, y, z, t) = one - four*t
    df1y(x, y, z, t) = one - four*t
    df1z(x, y, z, t) = one - four*t

    df2x(x, y, z, t) = -one + four*x
    df2y(x, y, z, t) = zero
    df2z(x, y, z, t) = zero

    df3x(x, y, z, t) = zero
    df3y(x, y, z, t) = -one + four*y 
    df3z(x, y, z, t) = zero

    df4x(x, y, z, t) = zero
    df4y(x, y, z, t) = zero 
    df4z(x, y, z, t) = -one + four*z

    df5x(x, y, z, t) = zero
    df5y(x, y, z, t) = four * z 
    df5z(x, y, z, t) = four * y

    df6x(x, y, z, t) = four * z
    df6y(x, y, z, t) = zero 
    df6z(x, y, z, t) = four * x

    df7x(x, y, z, t) = four * y
    df7y(x, y, z, t) = four * x 
    df7z(x, y, z, t) = zero

    df8x(x, y, z, t) = -four * z
    df8y(x, y, z, t) = -four * z 
    df8z(x, y, z, t) =  four * (t-z)

    df9x(x, y, z, t) = -four * y
    df9y(x, y, z, t) =  four * (t-y) 
    df9z(x, y, z, t) = -four * y

   df10x(x, y, z, t) =  four * (t-x)
   df10y(x, y, z, t) = -four * x 
   df10z(x, y, z, t) = -four * x

!     Degree 5; 15 Points;  Stroud: p. 315.

      vol = one/six

      sq = sqrt(three*five)

      r  = one/four;                              a = (vol * 16)/135

      s1 = (7 - sq)/34;    t1 = (13 + 3*sq)/34;  b1 = vol * (2665 + 14*sq)/37800

      s2 = (7 + sq)/34;    t2 = (13 - 3*sq)/34;  b2 = vol * (2665 - 14*sq)/37800

      u  = (10 - 2*sq)/40;  v = (10 + 2*sq)/40;   c = (vol * 20)/378


      xx(1) = r;    yy(1) = r;    zz(1) = r;    p(1) = a

      xx(2) = s1;   yy(2) = s1;   zz(2) = s1;   p(2) = b1
      xx(3) = s1;   yy(3) = s1;   zz(3) = t1;   p(3) = b1
      xx(4) = s1;   yy(4) = t1;   zz(4) = s1;   p(4) = b1
      xx(5) = t1;   yy(5) = s1;   zz(5) = s1;   p(5) = b1

      xx(6) = s2;   yy(6) = s2;   zz(6) = s2;   p(6) = b2
      xx(7) = s2;   yy(7) = s2;   zz(7) = t2;   p(7) = b2
      xx(8) = s2;   yy(8) = t2;   zz(8) = s2;   p(8) = b2
      xx(9) = t2;   yy(9) = s2;   zz(9) = s2;   p(9) = b2

      xx(10) = u;   yy(10) = u;   zz(10) = v;   p(10) = c
      xx(11) = u;   yy(11) = v;   zz(11) = u;   p(11) = c
      xx(12) = v;   yy(12) = u;   zz(12) = u;   p(12) = c
      xx(13) = v;   yy(13) = v;   zz(13) = u;   p(13) = c
      xx(14) = v;   yy(14) = u;   zz(14) = v;   p(14) = c
      xx(15) = u;   yy(15) = v;   zz(15) = v;   p(15) = c

      DO j = 1, l_G

         tj = one - xx(j) - yy(j) - zz(j) !This variable is a parameter of shape functions

            w(1, j) =  f1 (xx(j), yy(j), zz(j), tj)
         d(1, 1, j) = df1x(xx(j), yy(j), zz(j), tj)
         d(2, 1, j) = df1y(xx(j), yy(j), zz(j), tj)
         d(3, 1, j) = df1z(xx(j), yy(j), zz(j), tj)

            w(2, j) =  f2 (xx(j), yy(j), zz(j), tj)
         d(1, 2, j) = df2x(xx(j), yy(j), zz(j), tj)
         d(2, 2, j) = df2y(xx(j), yy(j), zz(j), tj)
         d(3, 2, j) = df2z(xx(j), yy(j), zz(j), tj)

            w(3, j) =  f3 (xx(j), yy(j), zz(j), tj)
         d(1, 3, j) = df3x(xx(j), yy(j), zz(j), tj)
         d(2, 3, j) = df3y(xx(j), yy(j), zz(j), tj)
         d(3, 3, j) = df3z(xx(j), yy(j), zz(j), tj)

            w(4, j) =  f4 (xx(j), yy(j), zz(j), tj)
         d(1, 4, j) = df4x(xx(j), yy(j), zz(j), tj)
         d(2, 4, j) = df4y(xx(j), yy(j), zz(j), tj)
         d(3, 4, j) = df4z(xx(j), yy(j), zz(j), tj)

            w(5, j) =  f5 (xx(j), yy(j), zz(j), tj)
         d(1, 5, j) = df5x(xx(j), yy(j), zz(j), tj)
         d(2, 5, j) = df5y(xx(j), yy(j), zz(j), tj)
         d(3, 5, j) = df5z(xx(j), yy(j), zz(j), tj)

            w(6, j) =  f6 (xx(j), yy(j), zz(j), tj)
         d(1, 6, j) = df6x(xx(j), yy(j), zz(j), tj)
         d(2, 6, j) = df6y(xx(j), yy(j), zz(j), tj)
         d(3, 6, j) = df6z(xx(j), yy(j), zz(j), tj)

            w(7, j) =  f7 (xx(j), yy(j), zz(j), tj)
         d(1, 7, j) = df7x(xx(j), yy(j), zz(j), tj)
         d(2, 7, j) = df7y(xx(j), yy(j), zz(j), tj)
         d(3, 7, j) = df7z(xx(j), yy(j), zz(j), tj)

            w(8, j) =  f8 (xx(j), yy(j), zz(j), tj)
         d(1, 8, j) = df8x(xx(j), yy(j), zz(j), tj)
         d(2, 8, j) = df8y(xx(j), yy(j), zz(j), tj)
         d(3, 8, j) = df8z(xx(j), yy(j), zz(j), tj)

            w(9, j) =  f9 (xx(j), yy(j), zz(j), tj)
         d(1, 9, j) = df9x(xx(j), yy(j), zz(j), tj)
         d(2, 9, j) = df9y(xx(j), yy(j), zz(j), tj)
         d(3, 9, j) = df9z(xx(j), yy(j), zz(j), tj)

           w(10, j) =  f10 (xx(j), yy(j), zz(j), tj)
        d(1, 10, j) = df10x(xx(j), yy(j), zz(j), tj)
        d(2, 10, j) = df10y(xx(j), yy(j), zz(j), tj)
        d(3, 10, j) = df10z(xx(j), yy(j), zz(j), tj)

      ENDDO

   END SUBROUTINE element_3d

!------------------------------------------------------------------------------

   SUBROUTINE element_2d_p2 (w, d, p)

!     triangular element with quadratic interpolation
!     and seven Gauss integration points

!        w(n_ws, l_Gs) : values of shape functions at Gauss points
!     d(2, n_ws, l_Gs) : derivatives values of shape functions at Gauss points
!             p(l_Gs)  : weight for Gaussian quadrature at Gauss points

      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(   n_ws, l_Gs), INTENT(OUT) :: w
      REAL(KIND=8), DIMENSION(2, n_ws, l_Gs), INTENT(OUT) :: d
      REAL(KIND=8), DIMENSION(l_Gs),          INTENT(OUT) :: p

      REAL(KIND=8), DIMENSION(l_Gs) :: xx, yy
      INTEGER :: j

      REAL(KIND=8) :: zero = 0,  half  = 0.5,  one  = 1,  &
		      two  = 2,  three = 3,    four = 4,  &
		      five = 5,   nine = 9

      REAL(KIND=8) ::  f1,   f2,   f3,   f4,   f5,   f6,  &
                      df1x, df2x, df3x, df4x, df5x, df6x, &
                      df1y, df2y, df3y, df4y, df5y, df6y, &
                      x, y, r, a,  s1, s2, t1, t2, b1, b2,  area, sq


      f1(x, y) = (half - x - y) * (one - x - y) * two
      f2(x, y) = x * (x - half) * two
      f3(x, y) = y * (y - half) * two
      f4(x, y) = x * y * four
      f5(x, y) = y * (one - x - y) * four
      f6(x, y) = x * (one - x - y) * four

      df1x(x, y) = -three + four * (x + y)
      df2x(x, y) = (two*x - half) * two
      df3x(x, y) = zero
      df4x(x, y) =  y * four
      df5x(x, y) = -y * four
      df6x(x, y) = (one - two*x - y) * four

      df1y(x, y) = -three + four * (x + y)
      df2y(x, y) = zero
      df3y(x, y) = (two*y - half) * two
      df4y(x, y) =  x * four
      df5y(x, y) = (one - x - two*y) * four
      df6y(x, y) = -x * four

!     Degree 5; 7 Points;  Stroud: p. 314, Approximate calculation of
!                          Multiple integrals (Prentice--Hall), 1971.

      area = one/two

      sq = sqrt(three*five)

      r  = one/three;                          a = area * nine/40

      s1 = (6 - sq)/21;  t1 = (9 + 2*sq)/21;  b1 = area * (155 - sq)/1200

      s2 = (6 + sq)/21;  t2 = (9 - 2*sq)/21;  b2 = area * (155 + sq)/1200

      xx(1) = r;    yy(1) = r;    p(1) = a

      xx(2) = s1;   yy(2) = s1;   p(2) = b1
      xx(3) = s1;   yy(3) = t1;   p(3) = b1
      xx(4) = t1;   yy(4) = s1;   p(4) = b1

      xx(5) = s2;   yy(5) = s2;   p(5) = b2
      xx(6) = s2;   yy(6) = t2;   p(6) = b2
      xx(7) = t2;   yy(7) = s2;   p(7) = b2


      DO j = 1, l_Gs

            w(1, j) =  f1 (xx(j), yy(j))
         d(1, 1, j) = df1x(xx(j), yy(j))
         d(2, 1, j) = df1y(xx(j), yy(j))

            w(2, j) =  f2 (xx(j), yy(j))
         d(1, 2, j) = df2x(xx(j), yy(j))
         d(2, 2, j) = df2y(xx(j), yy(j))

            w(3, j) =  f3 (xx(j), yy(j))
         d(1, 3, j) = df3x(xx(j), yy(j))
         d(2, 3, j) = df3y(xx(j), yy(j))

            w(4, j) =  f4 (xx(j), yy(j))
         d(1, 4, j) = df4x(xx(j), yy(j))
         d(2, 4, j) = df4y(xx(j), yy(j))

            w(5, j) =  f5 (xx(j), yy(j))
         d(1, 5, j) = df5x(xx(j), yy(j))
         d(2, 5, j) = df5y(xx(j), yy(j))

            w(6, j) =  f6 (xx(j), yy(j))
         d(1, 6, j) = df6x(xx(j), yy(j))
         d(2, 6, j) = df6y(xx(j), yy(j))

      ENDDO

   END SUBROUTINE element_2d_p2

END SUBROUTINE Gauss_gen

END MODULE Gauss_points
