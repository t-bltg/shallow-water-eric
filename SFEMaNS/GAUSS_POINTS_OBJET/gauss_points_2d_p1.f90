!
!Authors: Jean-Luc Guermond, Lugi Quartapelle, Copyright 1994
!
MODULE mod_gauss_points_2d_p1 
  PRIVATE

  PUBLIC gauss_points_2d_p1

CONTAINS

  SUBROUTINE gauss_points_2d_p1(mesh)

    USE def_type_mesh

    IMPLICIT NONE

    TYPE(mesh_type), TARGET :: mesh 

    INTEGER,                       POINTER :: me, mes
    INTEGER, DIMENSION(:, :),      POINTER :: jj
    INTEGER, DIMENSION(:, :),      POINTER :: js
    REAL(KIND=8), DIMENSION(:, :), POINTER :: rr

    INTEGER, PARAMETER :: k_d = 2,  n_w  = 3,  l_G  = 3,   &
         n_ws = 2,  l_Gs = 2

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
    INTEGER :: m, l, k, k1, n,  ms, ls, face, n1, n2, cote
    REAL(KIND=8), DIMENSION(2) :: rnor, rsd ! JLG April 2009
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

    ! JLG (April 2009)
    IF (mesh%edge_stab) THEN
       ALLOCATE(mesh%gauss%dwni(n_w, l_Gs, 2, mesh%mi))
       ALLOCATE(mesh%gauss%rji(l_Gs, mesh%mi))
       ALLOCATE(mesh%gauss%rnormsi(k_d, l_Gs, 2, mesh%mi)) !JLG June 4 2012
       ALLOCATE(mesh%gauss%wwsi(n_ws, l_gs)) !JLG June 4 2012
    END IF
    ! JLG (April 2009)

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

    CALL element_2d(ww, dd, pp)

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

    CALL element_1d(wws, dds, pps)

    !  surface elements
    DO ms = 1, mes
       m = mesh%neighs(ms)
       !Il faut savoir quelle face est la bonne
       DO n = 1, n_w
          IF (MINVAL(ABS(js(:,ms)-jj(n,m)))==0) CYCLE
          face = n 
       END DO
       CALL element_2d_boundary (face, dd, ww_s)
       n1 = MODULO(face,3)+1; n2 = MODULO(face+1,3)+1
       IF (js(1,ms)==jj(n1,m) .AND. js(2,ms)==jj(n2,m)) THEN
          nface(1)=n1; nface(2)=n2
          ngauss(1)=1; ngauss(2)=2
       ELSE IF (js(1,ms)==jj(n2,m) .AND. js(2,ms)==jj(n1,m)) THEN
          nface(1)=n2; nface(2)=n1
          ngauss(1)=2; ngauss(2)=1
       ELSE
          WRITE(*,*) ' BUG : gauss_points_2d_p1'
          STOP
       END IF

       DO ls = 1, l_Gs
          !write(*,*) wws(:,ls)
          !write(*,*) ww_s(nface,ngauss(ls))
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

             !dw_s(1, n, ls, ms)   &
             !   = (+ dd(1,nface(n),ngauss(ls))*dr(2,2) - dd(2,nface(n),ngauss(ls))*dr(2,1))/rjac
             !dw_s(2, n, ls, ms)   &
             !   = (- dd(1,nface(n),ngauss(ls))*dr(1,2) + dd(2,nface(n),ngauss(ls))*dr(1,1))/rjac
             !write(*,*) ms, n, dw_s(1, n, ls, ms), dw_s(2, n, ls, ms), mesh%sides(ms)
             !write(*,*) ms, n, rr(1,js(n,ms)), rr(2,js(n,ms))

          ENDDO

       ENDDO

    ENDDO

    !CALL element_1d(wws, dds, pps)
    DO ms = 1, mes

       DO ls = 1, l_Gs

          DO k = 1, k_d
             rs = rr(k, js(:,ms))
             drs(1, k) = SUM(rs * dds(1,:,ls))
          ENDDO

          rjacs = SQRT( drs(1,1)**2 + drs(1,2)**2 )

          rnorms(1, ls, ms) = - drs(1,2)/rjacs   ! outward normal
          rnorms(2, ls, ms) = + drs(1,1)/rjacs   ! outward normal
          !erreur de signe corrigee le 23 mai 2001. La normale etait interieure
          rjs(ls, ms) = rjacs * pps(ls)

          m = mesh%neighs(ms)
          !Il faut savoir quelle face est la bonne
          DO n = 1, n_w
             IF (MINVAL(ABS(js(:,ms)-jj(n,m)))==0) CYCLE
             face = n 
          END DO
          rs = rr(:,jj(face,m)) - (rr(:,js(1,ms))+rr(:,js(2,ms)))/2
          x = SUM(rnorms(:,ls,ms)*rs)
          IF (x>0) THEN
             rnorms(:,ls,ms) = -rnorms(:,ls,ms)
             !write(*,*) 'retournement de n', ms, ls
             !erreur de signe corrigee le 7 juin 2004.
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

    !-----------------------------------------------------------------------------
    ! CELL INTERFACES (JLG, April 2009)
    IF (mesh%edge_stab) THEN
       CALL element_1d(mesh%gauss%wwsi, dds, pps) !JLG June 4 2012
       DO ms = 1, mesh%mi
          DO cote = 1, 2
             m = mesh%neighi(cote,ms)
             !Il faut savoir quelle face est la bonne
             DO n = 1, n_w
                IF (MINVAL(ABS(mesh%jjsi(:,ms)-jj(n,m)))==0) CYCLE
                face = n 
             END DO
             CALL element_2d_boundary (face, dd, ww_s)
             n1 = MODULO(face,3)+1; n2 = MODULO(face+1,3)+1
             IF (mesh%jjsi(1,ms)==jj(n1,m) .AND. mesh%jjsi(2,ms)==jj(n2,m)) THEN
                nface(1)=n1; nface(2)=n2
                ngauss(1)=1; ngauss(2)=2
             ELSE IF (mesh%jjsi(1,ms)==jj(n2,m) .AND. mesh%jjsi(2,ms)==jj(n1,m)) THEN
                nface(1)=n2; nface(2)=n1
                ngauss(1)=2; ngauss(2)=1
             ELSE
                WRITE(*,*) ' BUG : gauss_points_2d_p1'
                STOP
             END IF
             DO ls = 1, l_Gs
                ! Compute normal vector
                DO k = 1, k_d
                   rs = rr(k, mesh%jjsi(:,ms))
                   drs(1, k) = SUM(rs * dds(1,:,ls))
                ENDDO
                rjacs = SQRT( drs(1,1)**2 + drs(1,2)**2 )
                mesh%gauss%rji(ls, ms) = rjacs * pps(ls)
                rnor(1) =  drs(1,2)/rjacs
                rnor(2) = - drs(1,1)/rjacs
                mesh%gauss%rnormsi(1, ls, cote, ms) = rnor(1) !JLG, June 4 2012
                mesh%gauss%rnormsi(2, ls, cote, ms) = rnor(2) !JLG, June 4 2012
                rsd = rr(:,jj(face,m)) - (rr(:,mesh%jjsi(1,ms))+rr(:,mesh%jjsi(2,ms)))/2
                x = SUM(rnor(:)*rsd)
                IF (x>0) THEN ! Verify the normal is outward
                   rnor = -rnor
                   !WRITE(*,*) 'retournement de n', ms, ls, cote, x
                END IF
                ! End computation normal vector
                !WRITE(*,*) SUM((rr(:,mesh%jjsi(1,ms))-rr(:,mesh%jjsi(2,ms)))*rnor)
                DO k = 1, k_d
                   r = rr(k, jj(:,m))
                   DO k1 = 1, k_d
                      dr(k, k1) = SUM(r * dd(k1,:,ls))
                   ENDDO
                ENDDO
                rjac = (dr(1,1)*dr(2,2) - dr(1,2)*dr(2,1))
                DO n = 1, n_w
                   !mesh%gauss%dwi(1, n, ngauss(ls), cote, ms)   &
                   !     = (+ dd(1,n,ls)*dr(2,2) - dd(2,n,ls)*dr(2,1))/rjac
                   !mesh%gauss%dwi(2, n, ngauss(ls), cote, ms)   &
                   !     = (- dd(1,n,ls)*dr(1,2) + dd(2,n,ls)*dr(1,1))/rjac
                   mesh%gauss%dwni(n, ngauss(ls), cote, ms)   &
                        = rnor(1)*(+ dd(1,n,ls)*dr(2,2) - dd(2,n,ls)*dr(2,1))/rjac &
                        + rnor(2)*(- dd(1,n,ls)*dr(1,2) + dd(2,n,ls)*dr(1,1))/rjac
                ENDDO
             ENDDO
          END DO
       ENDDO


       DO ms = 1, mesh%mi
          DO cote = 1, 2
             m = mesh%neighi(cote,ms)
             DO ls = 1, l_Gs
                x = ABS(SQRT(SUM(mesh%gauss%rnormsi(:, ls, cote, ms)**2))-1.d0)
                !write(*,*) ' x ',x, ABS(SQRT(SUM(mesh%gauss%rnormsi(:, ls, cote, ms)**2)))
                IF (x.GE.1.d-13) THEN
                   write(*,*) 'Bug1 in normal '
                   STOP
                END IF
                rjac = (rr(1, mesh%jjsi(2,ms)) -  rr(1, mesh%jjsi(1,ms)))*mesh%gauss%rnormsi(1, ls, cote, ms) &
                     + (rr(2, mesh%jjsi(2,ms)) -  rr(2, mesh%jjsi(1,ms)))*mesh%gauss%rnormsi(2, ls, cote, ms)
                x = SQRT((rr(1, mesh%jjsi(2,ms)) -  rr(1, mesh%jjsi(1,ms)))**2 +&
                     (rr(2, mesh%jjsi(2,ms)) -  rr(2, mesh%jjsi(1,ms)))**2)
                !write(*,*) ' rjac/x ', rjac/x, x, rjac
                IF (ABS(rjac/x ).GE. 1.d-13)  THEN
                   write(*,*) 'Bug2 in normal '
                   STOP
                END IF
             END DO
          END DO
       END DO
    END IF
    !-----------------------------------------------------------------------------

    !TEST 
    !WRITE(*,*) MAXVAL(ABS(dw(:,:,1,mesh%neighi(1,:)) - mesh%gauss%dwi(:,:,1,1,:)))
    !WRITE(*,*) MAXVAL(ABS(dw(:,:,2,mesh%neighi(2,:)) - mesh%gauss%dwi(:,:,1,2,:)))
    !WRITE(*,*) MAXVAL(ABS(jj(:,mesh%neighi(1,:))-mesh%jji(:,1,:)))
    !WRITE(*,*) MAXVAL(ABS(jj(:,mesh%neighi(2,:))-mesh%jji(:,2,:)))
    !TEST

    PRINT*, 'end of gen_Gauss'

  CONTAINS

    SUBROUTINE element_2d (w, d, p)

      !     triangular element with linear interpolation
      !     and three Gauss integration points

      !        w(n_w, l_G) : values of shape functions at Gauss points
      !     d(2, n_w, l_G) : derivatives values of shape functions at Gauss points
      !             p(l_G) : weight for Gaussian quadrature at Gauss points

      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(   n_w, l_G), INTENT(OUT) :: w
      REAL(KIND=8), DIMENSION(2, n_w, l_G), INTENT(OUT) :: d
      REAL(KIND=8), DIMENSION(l_G),         INTENT(OUT) :: p

      REAL(KIND=8), DIMENSION(l_G) :: xx, yy
      INTEGER :: j

      REAL(KIND=8) :: zero = 0,  one = 1,  two = 2,  three = 3,  six = 6

      REAL(KIND=8) :: f1, f2, f3, x, y
      f1(x, y) = one - x - y
      f2(x, y) = x
      f3(x, y) = y

      xx(1) = one/six;  xx(2) = two/three;  xx(3) = one/six
      yy(1) = one/six;  yy(2) = one/six;    yy(3) = two/three

      DO j = 1, l_G

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

    SUBROUTINE element_2d_boundary (face, d, w)

      !     triangular element with linear interpolation
      !     and three Gauss integration points

      !        w(n_w, l_G) : values of shape functions at Gauss points
      !     d(2, n_w, l_G) : derivatives values of shape functions at Gauss points
      !             p(l_G) : weight for Gaussian quadrature at Gauss points

      IMPLICIT NONE

      INTEGER,                                INTENT(IN)  :: face
      REAL(KIND=8), DIMENSION(2, n_w,  l_Gs), INTENT(OUT) :: d
      REAL(KIND=8), DIMENSION(   n_w, l_G),   INTENT(OUT) :: w

      REAL(KIND=8), DIMENSION(l_Gs) :: xx, yy
      INTEGER :: j

      REAL(KIND=8) :: zero = 0,  one = 1, three=3, half = 0.5d0

      REAL(KIND=8) :: f1, f2, f3, x, y
      f1(x, y) = 1-x-y
      f2(x, y) = x
      f3(x, y) = y

      IF (face==1) THEN
         xx(1) = half + half/SQRT(three)
         yy(1) = half - half/SQRT(three)

         xx(2) = half - half/SQRT(three)
         yy(2) = half + half/SQRT(three)
      ELSE IF (face==2) THEN
         xx(1) = zero
         yy(1) = half + half/SQRT(three)

         xx(2) = zero
         yy(2) = half - half/SQRT(three)
      ELSE 
         xx(1) = half - half/SQRT(three)
         yy(1) = zero

         xx(2) = half + half/SQRT(three)
         yy(2) = zero
      END IF

      DO j = 1, l_Gs
         w(1,j) = f1(xx(j), yy(j))
         d(1, 1, j) = -one
         d(2, 1, j) = -one

         w(2,j) = f2(xx(j), yy(j))
         d(1, 2, j) = one
         d(2, 2, j) = zero

         w(3,j) = f3(xx(j), yy(j))
         d(1, 3, j) = zero
         d(2, 3, j) = one
      ENDDO

    END SUBROUTINE element_2d_boundary

    !------------------------------------------------------------------------------

    SUBROUTINE element_1d (w, d, p)

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

      REAL(KIND=8) ::  one = 1,  two = 2,  three = 3

      REAL(KIND=8) :: f1, f2, x
      f1(x) = (one - x)/two
      f2(x) = (x + one)/two

      xx(1) = - one/SQRT(three)
      xx(2) = + one/SQRT(three)

      DO j = 1, l_Gs

         w(1, j) = f1(xx(j))
         d(1, 1, j) = - one/two

         w(2, j) = f2(xx(j))
         d(1, 2, j) = + one/two

         p(j) = one

      ENDDO

    END SUBROUTINE element_1d


  END SUBROUTINE gauss_points_2d_p1

END MODULE mod_gauss_points_2d_p1

