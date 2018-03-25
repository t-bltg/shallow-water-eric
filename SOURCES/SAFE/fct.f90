MODULE fct
CONTAINS
SUBROUTINE FCT_generic(unext,maxn,minn,mat,dg)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: dg, maxn, minn
    TYPE(matrice_bloc),         INTENT(IN)  :: mat
    REAL(KIND=8), DIMENSION(:), INTENT(INOUT):: unext
    REAL(KIND=8), DIMENSION(SIZE(unext))       :: Qplus, Qminus, Pplus, Pminus, Rplus, Rminus
    REAL(KIND=8), PARAMETER :: smallplus = 1.d-15, smallminus = -1.d-15
    REAL(KIND=8) :: x, fij, lij
    INTEGER      :: i, j, p 
    Qplus  = dg*(maxn-unext)
    Qminus = dg*(minn-unext)
    Pplus  = smallplus
    Pminus = smallminus
    DO i = 1, SIZE(unext)
       DO p = mat%ia(i), mat%ia(i+1) - 1
          j = mat%ja(p)
          fij = mat%aa(p)
          IF (fij.GE.0.d0) THEN
             Pplus(i)  = Pplus(i) + fij
          ELSE
             Pminus(i) = Pminus(i) + fij
          END IF
       END DO
       Rplus(i)  =  MIN(Qplus(i)/Pplus(i),1.d0)
       Rminus(i) =  MIN(Qminus(i)/Pminus(i),1.d0)
    END DO
    DO i = 1, SIZE(unext)
       x = 0.d0
       DO p = mat%ia(i), mat%ia(i+1) - 1
          j = mat%ja(p)
          fij = mat%aa(p)
          IF (fij.GE.0.d0) THEN
             lij = MIN(Rplus(i),Rminus(j))
          ELSE
             lij = MIN(Rminus(i),Rplus(j))
          END IF
          x = x + lij*fij
       END DO
       unext(i) = unext(i) + x/dg(i)
    END DO
  END SUBROUTINE FCT_generic

  SUBROUTINE FCT_positivity(unext,mat,dg)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: dg
    TYPE(matrice_bloc),         INTENT(IN)  :: mat
    REAL(KIND=8), DIMENSION(:), INTENT(INOUT):: unext
    REAL(KIND=8), DIMENSION(SIZE(unext))    :: Qminus, Pminus, Rminus
    REAL(KIND=8), PARAMETER :: smallminus = -1.d-15
    REAL(KIND=8) :: x, fij, lij
    INTEGER      :: i, j, p 

    Qminus = dg*(0.d0-unext)
    Pminus = smallminus
    DO i = 1, SIZE(unext)
       DO p = mat%ia(i), mat%ia(i+1) - 1
          j = mat%ja(p)
          fij = mat%aa(p)
          IF (fij.LT.0.d0) THEN
             Pminus(i) = Pminus(i) + fij
          END IF
       END DO
       Rminus(i) =  MIN(Qminus(i)/Pminus(i),1.d0)
    END DO
    DO i = 1, SIZE(unext)
       x = 0.d0
       DO p = mat%ia(i), mat%ia(i+1) - 1
          j = mat%ja(p)
          fij = mat%aa(p)
          IF (fij.GE.0.d0) THEN
             lij = Rminus(j)
          ELSE
             lij = Rminus(i)
          END IF
          x = x + lij*fij
       END DO
       unext(i) = unext(i) + x/dg(i)
    END DO
  END SUBROUTINE FCT_positivity
END MODULE fct
