MODULE bubble_elim

IMPLICIT NONE 

PRIVATE
PUBLIC :: elim_bubble_rhs_p1, elim_bubble_matrix_p1, solve_bubble_p1, &
          extract_matrix_p1, extract_rhs_p1, reconstruct_p1

CONTAINS

SUBROUTINE elim_bubble_rhs_p1(jj, bcd,  ff)
!=======================================
!  f_H  =  f_H  -  B D^(-1) f_h^H

INTEGER,      DIMENSION(:,:), INTENT(IN)    :: jj
REAL(KIND=8), DIMENSION(:,:), INTENT(IN)    :: bcd
REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: ff

REAL(KIND=8), DIMENSION(SIZE(jj,2)) :: f_b 

INTEGER:: me, np, kd, m

me = SIZE(jj,2)
np = SIZE(ff) - me
kd = (SIZE(bcd,2) - 1)/2 - 1 

f_b = ff(np+1:np+me)/bcd(:,kd+2)

DO m = 1, me
  ff(jj(:,m)) = ff(jj(:,m)) - bcd(m,1:kd+1)*f_b(m) 
END DO

END SUBROUTINE elim_bubble_rhs_p1

SUBROUTINE elim_bubble_matrix_p1(jj, bcd, ia, ja, aa)
!==================================================
!  AA  =  AA  -  B D^(-1) C

INTEGER,      DIMENSION(:,:), INTENT(IN)    :: jj
REAL(KIND=8), DIMENSION(:,:), INTENT(IN)    :: bcd
INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia, ja 
REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: aa 

REAL(KIND=8) :: x
INTEGER:: me, np, kd, i, j, ni, nj, m, p, nw

me = SIZE(jj,2)
np = SIZE(ia) - 1 
kd = (SIZE(bcd,2) - 1)/2 - 1
nw = SIZE(jj,1)

DO m = 1, me
   DO ni = 1, nw;  i = jj(ni,m)
      DO nj = 1, nw;  j = jj(nj,m)
         x = bcd(m,ni)*bcd(m,kd+2+nj)/bcd(m,kd+2)
         DO p = ia(i),  ia(i+1) - 1
            IF (ja(p) == j) THEN;  aa(p) = aa(p) - x;  EXIT;  ENDIF
         ENDDO
      END DO
   END DO
END DO

END SUBROUTINE elim_bubble_matrix_p1

SUBROUTINE solve_bubble_p1(jj, bcd, ff, uu)
!==================================================
!  u_h^H  =  D^(-1)(F_h^H  -  C u_H)

INTEGER,      DIMENSION(:,:), INTENT(IN)    :: jj
REAL(KIND=8), DIMENSION(:,:), INTENT(IN)    :: bcd
REAL(KIND=8), DIMENSION(:),   INTENT(IN)    :: ff 
REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: uu 

REAL(KIND=8) :: x
INTEGER:: me, np, kd, ni, m, nw

me = SIZE(jj,2)
np = SIZE(uu) - me 
kd = (SIZE(bcd,2) - 1)/2 - 1
nw = SIZE(jj,1)

DO m = 1, me
   x = ff(np+m)
   DO ni = 1, nw
         x = x - bcd(m,kd+2+ni)*uu(jj(ni,m)) 
   END DO
   uu(np+m) = x / bcd(m,kd+2)
END DO

END SUBROUTINE solve_bubble_p1

SUBROUTINE reconstruct_p1(jj, uu)
!==================================================
!  u_h  =  u_H + u_h^H

INTEGER,      DIMENSION(:,:), INTENT(IN)    :: jj
REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: uu

REAL(KIND=8) :: x, alpha
INTEGER:: me, np, ni, m, nw

me = SIZE(jj,2)
np = SIZE(uu) - me
nw = SIZE(jj,1)
alpha = 1.d0/nw

DO m = 1, me
   x = 0 
   DO ni = 1, nw
         x = x + uu(jj(ni,m))
   END DO
   uu(np+m) = uu(np+m) + x*alpha
END DO

END SUBROUTINE reconstruct_p1


SUBROUTINE extract_matrix_p1(jj, ia_b, ja_b, ia_test, ja_test, aa_b, aa, bcd)
!============================================================================
!
!  Create aa and bcd from aa_b for P1 interpolation
!

INTEGER,      DIMENSION(:,:), INTENT(IN)    :: jj
INTEGER,      POINTER,        DIMENSION(:)  :: ia_b, ja_b
INTEGER,      POINTER,        DIMENSION(:)  :: ia_test, ja_test
REAL(KIND=8), DIMENSION(:),   INTENT(OUT)   :: aa_b, aa 
REAL(KIND=8), DIMENSION(:,:), INTENT(OUT)   :: bcd

REAL(KIND=8) :: alpha, beta
INTEGER:: me, np, nw, i, j, n, m, p, p_b

me = SIZE(jj,2)
np = SIZE(ia_b) - 1 - me ! For P1 ONLY 
nw = SIZE(jj,1)

! First step: copy aa_b in aa and bcd
p = 0
DO i = 1, np
   DO p_b = ia_b(i), ia_b(i+1) - 1
      j = ja_b(p_b)
      IF (j .LE. np) THEN   ! Create aa
         p = p + 1
         aa(p) = aa_b(p_b) 
      ELSE                  ! Create bb
         m = j - np         ! For P1 interpolation ONLY 
         DO n = 1, nw
            IF (i .EQ. jj(n,m)) THEN
               bcd(m,n) = aa_b(p_b)
               EXIT
            END IF
         END DO
      END IF
   END DO
END DO

DO m = 1, me
   i = np + m
   DO p_b = ia_b(i), ia_b(i+1) - 1
      j = ja_b(p_b)
      IF (j .LE. np) THEN    ! Create cc
         DO n = 1, nw 
            IF (j .EQ. jj(n,m)) THEN
               bcd(m,nw+1+n) = aa_b(p_b)
               EXIT
            END IF
         END DO
      ELSE                   ! Create dd
         bcd(m,nw+1) = aa_b(p_b)
      END IF 
   END DO
END DO

!---Now we correct aa and bcd 

alpha = 1.d0/nw
beta = alpha**2

DO m = 1, me
   DO n = 1, nw
      i = jj(n,m)
      DO np = 1, nw
         j = jj(np,m) 
         DO p = ia_test(i), ia_test(i+1) - 1 
            IF (ja_test(p) .EQ. j) THEN
               aa(p) = aa(p) + alpha*(bcd(m,n)+bcd(m,nw+1+np))+beta*bcd(m,nw+1)
               EXIT
            END IF
         END DO
      END DO
   END DO
END DO

!New loop since do not touch bcd before aa is finished

DO m = 1, me 
   DO n = 1, nw
      bcd(m,n)      = bcd(m,n)      + alpha*bcd(m,nw+1)
      bcd(m,nw+1+n) = bcd(m,nw+1+n) + alpha*bcd(m,nw+1)
   END DO
END DO

END SUBROUTINE extract_matrix_p1 

SUBROUTINE extract_rhs_p1(jj, ff)
!=======================================
! reconstruct fH from fh
!
INTEGER,      DIMENSION(:,:), INTENT(IN)    :: jj
REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: ff

REAL(KIND=8) :: alpha
INTEGER:: me, np, nw, m, n, i

me = SIZE(jj,2)
np = SIZE(ff) - me
nw = SIZE(jj,1)
alpha = 1.d0/nw

DO m = 1, me
   DO n = 1, nw
      i = jj(n,m)
      ff(i) = ff(i) + alpha*ff(np+m)
   END DO
END DO

END SUBROUTINE extract_rhs_p1

END MODULE bubble_elim
