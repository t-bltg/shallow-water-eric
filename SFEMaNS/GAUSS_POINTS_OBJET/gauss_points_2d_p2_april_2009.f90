!
!Authors: Jean-Luc Guermond, Lugi Quartapelle, Copyright 1994
!
MODULE mod_gauss_points_2d_p2

  PRIVATE

  PUBLIC gauss_points_2d_p2

CONTAINS

  SUBROUTINE gauss_points_2d_p2(mesh)

    USE def_type_mesh

    IMPLICIT NONE

    TYPE(mesh_type), TARGET :: mesh 

    INTEGER,                       POINTER :: me, mes
    INTEGER, DIMENSION(:, :),      POINTER :: jj
    INTEGER, DIMENSION(:, :),      POINTER :: js
    REAL(KIND=8), DIMENSION(:, :), POINTER :: rr

    INTEGER, PARAMETER :: k_d = 2,  n_w  = 6,  l_G  = 7,   &
         n_ws = 3,  l_Gs = 3

    REAL(KIND=8), DIMENSION(:, :),       POINTER :: ww
    REAL(KIND=8), DIMENSION(:, :),       POINTER :: wws
    REAL(KIND=8), DIMENSION(:, :, :, :), POINTER :: dw
    REAL(KIND=8), DIMENSION(:,    :, :), POINTER :: rnorms
    REAL(KIND=8), DIMENSION(:, :),       POINTER :: rj
    REAL(KIND=8), DIMENSION(:, :),       POINTER :: rjs
    REAL(KIND=8), DIMENSION(:, :, :, :), POINTER :: dw_s
    REAL(KIND=8), DIMENSION(:, :, :, :), POINTER :: dwps !special!
    REAL(KIND=8), DIMENSION(:, :, :, :), POINTER :: dws  !SPECIAL!

    REAL(KIND=8), DIMENSION(k_d, n_w,  l_G)  :: dd
    REAL(KIND=8), DIMENSION(1,   n_ws, l_Gs) :: dds
    REAL(KIND=8), DIMENSION(l_G)             :: pp
    REAL(KIND=8), DIMENSION(l_Gs)            :: pps
    REAL(KIND=8), DIMENSION(n_w)             :: r
    REAL(KIND=8), DIMENSION(n_ws)            :: rs
    REAL(KIND=8), DIMENSION(k_d, k_d)        :: dr
    REAL(KIND=8), DIMENSION(1,   k_d)        :: drs

    REAL(KIND=8) :: rjac, rjacs, x
    INTEGER :: m, l, k, k1, n,  ms, ls, face, n1, n2, sign
    !TEST
    REAL(KIND=8), DIMENSION(n_w,  l_G)  :: ww_s
    !TEST
    !NEW Mai 26, 2004
    INTEGER,      DIMENSION(n_ws)            :: nface
    INTEGER,      DIMENSION(l_Gs)            :: ngauss
    !NEW Mai 26, 2004
    me => mesh%me
    mes => mesh%mes
    jj => mesh%jj
    js => mesh%jjs
    rr => mesh%rr

    ALLOCATE(mesh%gauss%ww(n_w, l_g))
    ALLOCATE(mesh%gauss%wws(n_ws, l_gs))
    ALLOCATE(mesh%gauss%dw(k_d, n_w, l_G,  me ))
    ALLOCATE(mesh%gauss%rnorms(k_d,  l_Gs, mes))
    ALLOCATE(mesh%gauss%rj(l_G, me ))
    ALLOCATE(mesh%gauss%rjs(l_Gs, mes))
    !NEW Juin 8, 2004
    ALLOCATE(mesh%gauss%dw_s(k_d, n_w,  l_Gs,  mesh%mes))
    !NEW Juin 8, 2004
    ALLOCATE(mesh%gauss%dwps(1, n_ws, l_Gs, mes))
    ALLOCATE(mesh%gauss%dws(1, n_ws, l_Gs, mes))

    ww => mesh%gauss%ww
    wws => mesh%gauss%wws
    dw => mesh%gauss%dw
    rnorms => mesh%gauss%rnorms
    rj => mesh%gauss%rj
    rjs => mesh%gauss%rjs
    !NEW Mai 26, 2004
    dw_s => mesh%gauss%dw_s
    !NEW Mai 26, 2004
    dwps => mesh%gauss%dwps
    dws => mesh%gauss%dws

    mesh%gauss%k_d  = k_d
    mesh%gauss%n_w  = n_w
    mesh%gauss%l_G  = l_G
    mesh%gauss%n_ws = n_ws
    mesh%gauss%l_Gs = l_Gs

    !  evaluate and store the values of derivatives and of the
    !  jacobian determinant at Gauss points of all volume elements

    !  volume elements

    CALL element_2d_p2(ww, dd, pp)

    DO m = 1, me

       DO l = 1, l_G

          DO k = 1, k_d
             r = rr(k, jj(:,m))
             DO k1 = 1, k_d
                dr(k, k1) = SUM(r * dd(k1,:,l))
             ENDDO
          ENDDO

          rjac = dr(1,1)*dr(2,2) - dr(1,2)*dr(2,1)

          DO n = 1, n_w
             dw(1, n, l, m)   &
                  = (+ dd(1,n,l)*dr(2,2) - dd(2,n,l)*dr(2,1))/rjac
             dw(2, n, l, m)   &
                  = (- dd(1,n,l)*dr(1,2) + dd(2,n,l)*dr(1,1))/rjac
          ENDDO

          rj(l, m) = rjac * pp(l)

       ENDDO

    ENDDO

    !TEST
    CALL element_1d_p2(wws, dds, pps)
    !TEST

    !  surface elements   
    DO ms = 1, mes
       m = mesh%neighs(ms)
       !Il faut savoir quelle face est la bonne
       DO n = 1, 3
          IF (MINVAL(ABS(js(:,ms)-jj(n,m)))==0) CYCLE
          face = n 
       END DO
       CALL element_2d_p2_boundary (face, dd, ww_s)
       n1 = MODULO(face,3)+1; n2 = MODULO(face+1,3)+1
       IF (js(1,ms)==jj(n1,m) .AND. js(2,ms)==jj(n2,m)) THEN
          nface(1)=n1; nface(2)=n2; nface(3)=face+3
          ngauss(1)=1; ngauss(2)=2; ngauss(3)=3;
          sign = 1
       ELSE IF (js(1,ms)==jj(n2,m) .AND. js(2,ms)==jj(n1,m)) THEN
          nface(1)=n2; nface(2)=n1; nface(3)=face+3
          ngauss(1)=3; ngauss(2)=2; ngauss(3)=1;
          sign = -1
       ELSE
          WRITE(*,*) ' BUG : gauss_points_2d_p2'
          STOP
       END IF

       DO ls = 1, l_Gs

          !write(*,*) wws(:,ls)
          !write(*,*) ww_s(nface,ngauss(ls)), sign

          DO k = 1, k_d
             r = rr(k, jj(:,m))
             DO k1 = 1, k_d
                dr(k, k1) = SUM(r * dd(k1,:,ls))
             ENDDO
          ENDDO

          rjac = (dr(1,1)*dr(2,2) - dr(1,2)*dr(2,1))

          DO n = 1, n_w
             dw_s(1, n, ngauss(ls), ms)   &
                  = (+ dd(1,n,ls)*dr(2,2) - dd(2,n,ls)*dr(2,1))/rjac
             dw_s(2, n, ngauss(ls), ms)   &
                  = (- dd(1,n,ls)*dr(1,2) + dd(2,n,ls)*dr(1,1))/rjac
          ENDDO

       ENDDO

    ENDDO

    CALL element_1d_p2(wws, dds, pps)

    DO ms = 1, mes

       DO ls = 1, l_Gs


          DO k = 1, k_d
             rs = rr(k, js(:,ms))
             drs(1, k) = SUM(rs * dds(1,:,ls))
          ENDDO

          rjacs = SQRT( drs(1,1)**2 + drs(1,2)**2 )

          rnorms(1, ls, ms) = - drs(1,2)/rjacs   ! outward normal
          rnorms(2, ls, ms) = + drs(1,1)/rjacs   ! outward normal
          !erreur de signe corrigee le 10 Oct 2003. La normale etait interieure
          rjs(ls, ms) = rjacs * pps(ls)

          m = mesh%neighs(ms)
          !Il faut savoir quelle face est la bonne
          DO n = 1, n_w
             IF (MINVAL(ABS(js(:,ms)-jj(n,m)))==0) CYCLE
             face = n 
          END DO
          rs(1:k_d) = rr(:,jj(face,m)) - (rr(:,js(1,ms))+rr(:,js(2,ms)))/2
          x = SUM(rnorms(:,ls,ms)*rs(1:k_d))
          IF (x>0) THEN
             rnorms(:,ls,ms) = -rnorms(:,ls,ms)
             !write(*,*) 'retournement de n', ms, ls 
             !erreur de signe corrigee le 8 juin 2004
          END IF

       ENDDO

    ENDDO


    DO ms = 1, mes  ! necessary only for evaluating gradient
       ! tangential to the surface (ex COMMON Gauss_tan)

       DO ls = 1, l_Gs
          !TEST
          !dwps(1, :, ls, ms) = dds(1, :, ls) /  (rjs(ls, ms)/pps(ls))
          !do n = 1, n_ws
          !rjac = SUM(rnorms(:,ls,ms)*dw_s(:,n,ls,ms))
          !write(*,*) n, ' js ', js(:,ms), rjs(ls, ms)/pps(ls)
          !write(*,*) ms, ls, dw_s(1, n, ls, ms) - rjac*rnorms(1,ls,ms), dw_s(2, n, ls, ms) - rjac*rnorms(2,ls,ms)
          !write(*,*) ms, ls, dwps(1, n, ls, ms)*rnorms(2,ls,ms), -dwps(1, n, ls, ms)*rnorms(1,ls,ms), &
          !           SQRT(SUM(rnorms(:,ls,ms)**2)),rnorms(:,ls,ms)
          !write(*,*) ' '
          !end do
          !TEST
          dwps(1, :, ls, ms) = dds(1, :, ls) * pps(ls)
          dws(1, :, ls, ms) = dds(1, :, ls)

       ENDDO

    ENDDO

    PRINT*, 'end of gen_Gauss'

  CONTAINS

    SUBROUTINE element_2d_p2 (w, d, p)

      !     triangular element with quadratic interpolation
      !     and seven Gauss integration points

      !        w(n_w, l_G) : values of shape functions at Gauss points
      !     d(2, n_w, l_G) : derivatives values of shape functions at Gauss points
      !             p(l_G) : weight for Gaussian quadrature at Gauss points

      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(   n_w, l_G), INTENT(OUT) :: w
      REAL(KIND=8), DIMENSION(2, n_w, l_G), INTENT(OUT) :: d
      REAL(KIND=8), DIMENSION(l_G),         INTENT(OUT) :: p

      REAL(KIND=8), DIMENSION(l_G) :: xx, yy
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

      sq = SQRT(three*five)

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


      DO j = 1, l_G

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

    SUBROUTINE element_2d_p2_boundary (face, d, w)

      !     triangular element with quadratic interpolation
      !     and seven Gauss integration points

      !     d(2, n_w, l_G) : derivatives values of shape functions at Gauss points

      IMPLICIT NONE

      INTEGER,                                INTENT(IN)  :: face
      REAL(KIND=8), DIMENSION(2, n_w,  l_Gs), INTENT(OUT) :: d
      REAL(KIND=8), DIMENSION(   n_w, l_G),   INTENT(OUT) :: w

      REAL(KIND=8), DIMENSION(l_G) :: xx, yy
      INTEGER :: j

      REAL(KIND=8) :: zero = 0,  half  = 0.5d0,  one  = 1,  &
           two  = 2,  three = 3,    four = 4,  &
           five = 5

      REAL(KIND=8) ::  f1,   f2,   f3,   f4,   f5,   f6,  &
           df1x, df2x, df3x, df4x, df5x, df6x, &
           df1y, df2y, df3y, df4y, df5y, df6y, &
           x, y


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

      IF (face==1) THEN 
         xx(1) = half + half*SQRT(three/five)
         yy(1) = half - half*SQRT(three/five)

         xx(2) = half
         yy(2) = half

         xx(3) = half - half*SQRT(three/five)
         yy(3) = half + half*SQRT(three/five)
      ELSE IF (face==2) THEN
         xx(1) = zero
         yy(1) = half + half*SQRT(three/five)

         xx(2) = zero
         yy(2) = half

         xx(3) = zero
         yy(3) = half - half*SQRT(three/five)
      ELSE 
         xx(1) = half - half*SQRT(three/five)
         yy(1) = zero

         xx(2) = half
         yy(2) = zero

         xx(3) = half + half*SQRT(three/five)
         yy(3) = zero
      END IF


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

    END SUBROUTINE element_2d_p2_boundary

    !------------------------------------------------------------------------------

    SUBROUTINE element_1d_p2(w, d, p)

      !     one-dimensional element with linear interpolation
      !     and two Gauss integration points

      !        w(n_w, l_G) : values of shape functions at Gauss points
      !     d(1, n_w, l_G) : derivatives values of shape functions at Gauss points
      !             p(l_G) : weight for Gaussian quadrature at Gauss points

      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(   n_ws, l_Gs), INTENT(OUT) :: w
      REAL(KIND=8), DIMENSION(1, n_ws, l_Gs), INTENT(OUT) :: d
      REAL(KIND=8), DIMENSION(l_Gs),          INTENT(OUT) :: p

      REAL(KIND=8), DIMENSION(l_Gs) :: xx
      INTEGER :: j

      REAL(KIND=8) :: zero = 0,  one = 1,  two = 2,  three = 3,  &
           five = 5, eight= 8,   nine = 9

      REAL(KIND=8) :: f1, f2, f3, df1, df2, df3, x

      f1(x) = (x - one)*x/two
      f2(x) = (x + one)*x/two
      f3(x) = (x + one)*(one - x)

      df1(x) = (two*x - one)/two
      df2(x) = (two*x + one)/two
      df3(x) = -two*x

      xx(1) = -SQRT(three/five)
      xx(2) =  zero
      xx(3) =  SQRT(three/five) 

      p(1)  =  five/nine
      p(2)  =  eight/nine 
      p(3)  =  five/nine 

      DO j = 1, l_Gs

         w(1, j) =  f1(xx(j))
         d(1, 1, j) = df1(xx(j))

         w(2, j) =  f2(xx(j))
         d(1, 2, j) = df2(xx(j))

         w(3, j) =  f3(xx(j))
         d(1, 3, j) = df3(xx(j))


      ENDDO

    END SUBROUTINE element_1d_p2

  END SUBROUTINE gauss_points_2d_p2

END MODULE mod_gauss_points_2d_p2

