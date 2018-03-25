MODULE fem_SGS_M

CONTAINS

SUBROUTINE qs_SGS_p1_2d_M(jj,alpha,ia,ja,aa)
   USE Gauss_points

   IMPLICIT NONE
   INTEGER,      DIMENSION(:,:), INTENT(IN)    :: jj
   REAL(KIND=8),                 INTENT(IN)    :: alpha
   INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia
   INTEGER,      DIMENSION(:),   INTENT(IN)    :: ja
   REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: aa

   REAL(KIND=8), DIMENSION(3)   :: xk, yk  
   REAL(KIND=8), DIMENSION(3,3) :: phi 
   INTEGER,      DIMENSION(3,4) :: ng 
   INTEGER,      DIMENSION(6)   :: ig 
   REAL(KIND=8), DIMENSION(2,6) :: dwt
!TESt
   REAL(KIND=8), DIMENSION(6) :: wt
!TESt
   INTEGER :: mg, mloc, m, n, k, l, ni, nj, i, j, p, n0
   REAL(KIND=8) :: h, al

   REAL(KIND=8) :: half  = 0.5,  one  = 1,  &
                   two  = 2,     four = 4 

   REAL(KIND=8) ::  f1,   f2,   f3,  x,  y 

!===Define the P1 base functions
      f1(x, y) = one - x - y
      f2(x, y) = x 
      f3(x, y) = y

!=======================Topology of the macro P1 element
!                                  3
!                                 
!                              5       4 
!                                   
!                          1       6       2
!
!===Build phi_ij = phi_i(x_{j+3})    1=< i =< 3, 1 =< j =< 3
xk(1) = 0.5;  yk(1) =0.5    !Point nb 4
xk(2) = 0.;   yk(2) =0.5    !Point nb 5
xk(3) = 0.5;  yk(3) =0.     !Point nb 6

DO nj = 1, 3
   phi(1,nj) = f1(xk(nj),yk(nj)) 
   phi(2,nj) = f2(xk(nj),yk(nj)) 
   phi(3,nj) = f3(xk(nj),yk(nj)) 
END DO

!===Build ng : Gives the global index from the local index 
!   index1 = local node index, index2 = local element number
ng(1,1) = 1; ng(2,1) = 6; ng(3,1) = 5
ng(1,2) = 2; ng(2,2) = 4; ng(3,2) = 6
ng(1,3) = 3; ng(2,3) = 5; ng(3,3) = 4
ng(1,4) = 6; ng(2,4) = 4; ng(3,4) = 5

IF (MOD(SIZE(jj,2),4).NE.0) THEN
   WRITE(*,*) 'BUG: qs_SGS_p1_2d_M'
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
      IF (mloc == 4) THEN
         n0 = 1 
      ELSE
         n0 = 2
      END IF


!===Compute the local mesh size
      h = 0
      DO l = 1, l_G
         h = h + rj(l,m)
      END DO 
      h = alpha * SQRT(h)

      DO l = 1, l_G
!===Build the weight
         al = h * rj(l,m)


!TEST
!DO n = 1, 3
!   DO k = 1, k_d
!     if (m==5)  write(*,*) 'l=',l,'k,n',k,n,'dw(k,n,l,m)',dw(k,n,l,m)
!   END DO
!END DO
!if(m==5 .AND. l==l_G) stop
!al = 1
!TEST

!===Build the gradient of the fluctuating base functions: dwt
         dwt = 0
         DO k = 1, k_d
   
            DO n = n0, 3; nj = ng(n,mloc)-3
              DO ni = 1, 3
                 dwt(k,ni) = dwt(k,ni) - phi(ni,nj)*dw(k,n,l,m)
!test
!dwt(k,ng(1,mloc)) = dw(k,1,l,m)
!test
              END DO
           END DO

            DO n = n0, 3; ni = ng(n,mloc) 
               dwt(k,ni) = dw(k,n,l,m)
            END DO
!TEST
!if(k_d==2) goto 100
!if (mloc==1) THEN
!m = 1  !On travaille sur un seul macro triangle : maill.test
!   dwt(k,1) = -0.5*dw(k,2,l,m) -0.5*dw(k,3,l,m)
!   dwt(k,2) = -0.5*dw(k,2,l,m) 
!   dwt(k,3) = -0.5*dw(k,3,l,m)
!   dwt(k,4) = 0.
!   dwt(k,5) = dw(k,3,l,m)
!   dwt(k,6) = dw(k,2,l,m)
!else if (mloc==2) THEN
!m = 2
!   dwt(k,1) = -0.5*dw(k,3,l,m) 
!   dwt(k,2) = -0.5*dw(k,2,l,m) -0.5*dw(k,3,l,m)
!   dwt(k,3) = -0.5*dw(k,2,l,m)
!   dwt(k,4) = dw(k,2,l,m)
!   dwt(k,5) = 0.
!   dwt(k,6) = dw(k,3,l,m)
!else if (mloc==3) THEN
!m = 3
!   dwt(k,1) = -0.5*dw(k,2,l,m) 
!   dwt(k,2) = -0.5*dw(k,3,l,m)
!   dwt(k,3) = -0.5*dw(k,2,l,m) -0.5*dw(k,3,l,m)
!   dwt(k,4) = dw(k,3,l,m)
!   dwt(k,5) = dw(k,2,l,m)
!   dwt(k,6) = 0.
!else 
!m = 4 
!   dwt(k,1) = -0.5*dw(k,1,l,m) -0.5*dw(k,3,l,m)
!   dwt(k,2) = -0.5*dw(k,1,l,m) -0.5*dw(k,2,l,m)
!   dwt(k,3) = -0.5*dw(k,3,l,m) -0.5*dw(k,2,l,m)
!   dwt(k,4) = dw(k,2,l,m)
!   dwt(k,5) = dw(k,3,l,m)
!   dwt(k,6) = dw(k,1,l,m)
!end if
!!dwt =0.
!100 continue
!!TEST

         END DO



!TEST
        wt = 0
        DO n = n0, 3; nj = ng(n,mloc)-3
           DO ni = 1, 3
              wt(ni) = wt(ni) - phi(ni,nj)*ww(n,l)
           END DO
        END DO

        DO n = n0, 3; ni = ng(n,mloc)
           wt(ni) = ww(n,l)
        END DO
!wt =0
!TEST
        

         DO ni = 1, 6
            i = ig(ni)
!i=ni !!!!!!!!!!!!!!!!!!1!            i = ig(ni)
            DO nj = 1, 6
               j = ig(nj)
!j=nj !!!!!!!!!!!!!!!!!111111!               j = ig(nj)

               x = 0
               DO k = 1, k_d
                  x = x + dwt(k,ni)*dwt(k,nj)
               END DO
               x =  al * x
!test
!               x =  al * (x+wt(ni)*wt(nj))
!test
               DO p = ia(i),  ia(i+1) - 1
                  IF (ja(p) == j) THEN;  aa(p) = aa(p) + x;  EXIT;  ENDIF
               ENDDO

            END DO
         END DO

      END DO
   END DO
END DO

END SUBROUTINE qs_SGS_p1_2d_M

SUBROUTINE ph_sgs_p1_2d(jj,uu,ff)
   IMPLICIT NONE
   INTEGER,      DIMENSION(:,:), INTENT(IN)    :: jj
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)    :: uu 
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT)   :: ff 

   INTEGER  :: i1, i2, i3, i4, i5, i6, mg
   
ff = uu

DO mg = 1, SIZE(jj,2)/4
!===Build the local connectivity array
   i1 = jj(1,4*(mg-1)+1) 
   i2 = jj(1,4*(mg-1)+2) 
   i3 = jj(1,4*(mg-1)+3) 
   i4 = jj(2,4*(mg-1)+4) 
   i5 = jj(3,4*(mg-1)+4) 
   i6 = jj(1,4*(mg-1)+4) 

   ff(i4) =  (uu(i2) + uu(i3))/2 
   ff(i5) =  (uu(i3) + uu(i1))/2 
   ff(i6) =  (uu(i1) + uu(i2))/2 

END DO

END SUBROUTINE ph_sgs_p1_2d

END MODULE fem_SGS_M
