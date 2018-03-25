MODULE fem_SGS_M

CONTAINS

SUBROUTINE qs_SGS_p2_2d_M(jj,alpha,ia,ja,aa)
   USE Gauss_points

   IMPLICIT NONE
   INTEGER,      DIMENSION(:,:), INTENT(IN)    :: jj
   REAL(KIND=8),                 INTENT(IN)    :: alpha
   INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia
   INTEGER,      DIMENSION(:),   INTENT(IN)    :: ja
   REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: aa

   REAL(KIND=8), DIMENSION(9)   :: xk, yk  
   REAL(KIND=8), DIMENSION(6,9) :: phi 
   INTEGER,      DIMENSION(6,4) :: ng 
   INTEGER,      DIMENSION(9)   :: ig 
   REAL(KIND=8), DIMENSION(2,9) :: dwt
   INTEGER :: mg, mloc, m, n, k, l, ni, nj, i, j, p
   REAL(KIND=8) :: h, al

   REAL(KIND=8) :: half  = 0.5,  one  = 1,  &
                   two  = 2,     four = 4 

   REAL(KIND=8) ::  f1,   f2,   f3,   f4,   f5,   f6,  &
                    x, y 

!===Define the P2 base functions
      f1(x, y) = (half - x - y) * (one - x - y) * two
      f2(x, y) = x * (x - half) * two
      f3(x, y) = y * (y - half) * two
      f4(x, y) = x * y * four
      f5(x, y) = y * (one - x - y) * four
      f6(x, y) = x * (one - x - y) * four

!=======================Topology of the macro P2 element
!                                  3
!                               13  15
!                              5  12   4 
!                            8   9  11  14
!                          1   7   6  10  2
!
!===Build phi_ij = phi_i(x_{j+6})    1=< i =< 6, 1 =< j =< 9
xk(1) = 0.25; yk(1) =0.     !Point nb 7
xk(2) = 0.;   yk(2) =0.25   !Point nb 8
xk(3) = 0.25; yk(3) =0.25   !Point nb 9
xk(4) = 0.75; yk(4) =0.     !point nb 10
xk(5) = 0.5;  yk(5) =0.25   !point nb 11
xk(6) = 0.25; yk(6) =0.5    !point nb 12
xk(7) = 0.;   yk(7) =0.75   !point nb 13
xk(8) = 0.75; yk(8) =0.25   !point nb 14
xk(9) = 0.25; yk(9) =0.75   !point nb 15

DO nj = 1, 9
   phi(1,nj) = f1(xk(nj),yk(nj)) 
   phi(2,nj) = f2(xk(nj),yk(nj)) 
   phi(3,nj) = f3(xk(nj),yk(nj)) 
   phi(4,nj) = f4(xk(nj),yk(nj)) 
   phi(5,nj) = f5(xk(nj),yk(nj)) 
   phi(6,nj) = f6(xk(nj),yk(nj)) 
END DO

!===Build ng : Gives the global index from the local index 
!   index1 = local node index, index2 = local element number
ng(1,1) = 1; ng(2,1) = 6; ng(3,1) = 5; ng(4,1) = 9; ng(5,1) = 8; ng(6,1) = 7
ng(1,2) = 2; ng(2,2) = 4; ng(3,2) = 6; ng(4,2) =11; ng(5,2) =10; ng(6,2) =14
ng(1,3) = 3; ng(2,3) = 5; ng(3,3) = 4; ng(4,3) =12; ng(5,3) =15; ng(6,3) =13
ng(1,4) = 6; ng(2,4) = 4; ng(3,4) = 5; ng(4,4) =12; ng(5,4) = 9; ng(6,4) =11

IF (MOD(SIZE(jj,2),4).NE.0) THEN
   WRITE(*,*) 'BUG: qs_SGS_p2_2d_M'
   WRITE(*,*) 'Mesh not compatible with hierarchical structure'
   STOP
END IF

DO mg = 1, SIZE(jj,2)/4
!===Build the local connectivity array
   ig(1) = jj(1,4*(mg-1)+1) 
   ig(2) = jj(1,4*(mg-1)+2) 
   ig(3) = jj(1,4*(mg-1)+3) 
   ig(4) = jj(2,4*(mg-1)+4) 
   ig(5) = jj(3,4*(mg-1)+4) 
   ig(6) = jj(1,4*(mg-1)+4) 

   DO mloc = 1, 4
      m = 4*(mg-1) + mloc !  mmg(mloc,mg)  !Global index of element
!===Build the local connectivity array
      DO n = 1, 3
         ig(n+6) = jj(n+3,m)
      END DO

!===Compute the local mesh size
      h = 0
      DO l = 1, l_G
         h = h + rj(l,m)
      END DO 
      h = alpha * SQRT(h)

      DO l = 1, l_G
!===Build the weight
         al = h * rj(l,m)

!===Build the gradient of the fluctuating bases functions: dwt
         dwt = 0
         DO k = 1, k_d
   
            DO n = 1, 3; nj = ng(n+3,mloc)-6
               DO ni = 1, 6 
                  dwt(k,ni) = dwt(k,ni) - phi(ni,nj)*dw(k,n+3,l,m)
               END DO
            END DO

            DO ni = 7, 9 
               dwt(k,ni) = dw(k,ni-3,l,m)
            END DO

         END DO

         DO ni = 1, 9
            i = ig(ni)
            DO nj = 1, 9
               j = ig(nj)
!TEST
!if (ni.le.6) then 
!   if (ni.gt.4) then 
!      dwt(:,ni) =0 
!   else
!      dwt(:,ni) = dw(:,ni,l,m)
!      i = jj(ni,m)
!   end if
!end if
!if (nj.le.6) then 
!   if (nj.gt.4) then 
!      dwt(:,nj) =0 
!   else
!      dwt(:,nj) = dw(:,nj,l,m)
!      j = jj(nj,m)
!   end if
!end if
!TEST
               x = 0
               DO k = 1, k_d
                  x = x + dwt(k,ni)*dwt(k,nj)
               END DO
               x =  al * x
               DO p = ia(i),  ia(i+1) - 1
                  IF (ja(p) == j) THEN;  aa(p) = aa(p) + x;  EXIT;  ENDIF
               ENDDO
            END DO
         END DO

      END DO
   END DO
END DO

END SUBROUTINE qs_SGS_p2_2d_M


SUBROUTINE ph_SGS_p2_2d(jj,uu,ff)
   USE Gauss_points

   IMPLICIT NONE
   INTEGER,      DIMENSION(:,:), INTENT(IN)  :: jj
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: uu 
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: ff 

   REAL(KIND=8), DIMENSION(9)   :: xk, yk  
   REAL(KIND=8), DIMENSION(6,9) :: phi 
   INTEGER,      DIMENSION(6,4) :: ng 
   INTEGER,      DIMENSION(9)   :: ig 
   INTEGER :: mg, mloc, m, n, ni, nj

   REAL(KIND=8) :: half  = 0.5,  one  = 1,  &
                   two  = 2,     four = 4 

   REAL(KIND=8) ::  f1,   f2,   f3,   f4,   f5,   f6,  &
                    x, y 

!===Define the P2 base functions
      f1(x, y) = (half - x - y) * (one - x - y) * two
      f2(x, y) = x * (x - half) * two
      f3(x, y) = y * (y - half) * two
      f4(x, y) = x * y * four
      f5(x, y) = y * (one - x - y) * four
      f6(x, y) = x * (one - x - y) * four

!=======================Topology of the macro P2 element
!                                  3
!                               13  15
!                              5  12   4 
!                            8   9  11  14
!                          1   7   6  10  2
!
!===Build phi_ij = phi_i(x_{j+6})    1=< i =< 6, 1 =< j =< 9
xk(1) = 0.25; yk(1) =0.     !Point nb 7
xk(2) = 0.;   yk(2) =0.25   !Point nb 8
xk(3) = 0.25; yk(3) =0.25   !Point nb 9
xk(4) = 0.75; yk(4) =0.     !point nb 10
xk(5) = 0.5;  yk(5) =0.25   !point nb 11
xk(6) = 0.25; yk(6) =0.5    !point nb 12
xk(7) = 0.;   yk(7) =0.75   !point nb 13
xk(8) = 0.75; yk(8) =0.25   !point nb 14
xk(9) = 0.25; yk(9) =0.75   !point nb 15

DO nj = 1, 9
   phi(1,nj) = f1(xk(nj),yk(nj)) 
   phi(2,nj) = f2(xk(nj),yk(nj)) 
   phi(3,nj) = f3(xk(nj),yk(nj)) 
   phi(4,nj) = f4(xk(nj),yk(nj)) 
   phi(5,nj) = f5(xk(nj),yk(nj)) 
   phi(6,nj) = f6(xk(nj),yk(nj)) 
END DO

!===Build ng : Gives the global index from the local index 
!   index1 = local node index, index2 = local element number
ng(1,1) = 1; ng(2,1) = 6; ng(3,1) = 5; ng(4,1) = 9; ng(5,1) = 8; ng(6,1) = 7
ng(1,2) = 2; ng(2,2) = 4; ng(3,2) = 6; ng(4,2) =11; ng(5,2) =10; ng(6,2) =14
ng(1,3) = 3; ng(2,3) = 5; ng(3,3) = 4; ng(4,3) =12; ng(5,3) =15; ng(6,3) =13
ng(1,4) = 6; ng(2,4) = 4; ng(3,4) = 5; ng(4,4) =12; ng(5,4) = 9; ng(6,4) =11

IF (MOD(SIZE(jj,2),4).NE.0) THEN
   WRITE(*,*) 'BUG: qs_SGS_p2_2d_M'
   WRITE(*,*) 'Mesh not compatible with hierarchical structure'
   STOP
END IF

DO mg = 1, SIZE(jj,2)/4
!===Build the local connectivity array
   ig(1) = jj(1,4*(mg-1)+1) 
   ig(2) = jj(1,4*(mg-1)+2) 
   ig(3) = jj(1,4*(mg-1)+3) 
   ig(4) = jj(2,4*(mg-1)+4) 
   ig(5) = jj(3,4*(mg-1)+4) 
   ig(6) = jj(1,4*(mg-1)+4) 

   DO n = 1, 6
      ff(ig(n)) = uu(ig(n))
   END DO

   DO mloc = 1, 4
      m = 4*(mg-1) + mloc !  mmg(mloc,mg)  !Global index of element
!===Build the local connectivity array
      DO n = 1, 3
         ig(n+6) = jj(n+3,m)
      END DO

      DO n = 1, 3
         ff(ig(n+6)) = 0
         DO ni =1, 6
            ff(ig(n+6)) = ff(ig(n+6)) + uu(ig(ni))*phi(ni,ng(n+3,mloc)-6)
         END DO
      END DO

   END DO
END DO

END SUBROUTINE ph_SGS_p2_2d

END MODULE fem_SGS_M
