MODULE fct
  USE matrix_type
  USE input_data
CONTAINS
  SUBROUTINE FCT_positivity(unext,mat,lij)
    IMPLICIT NONE
    TYPE(matrice_bloc),         INTENT(IN)  :: mat
    TYPE(matrice_bloc),         INTENT(OUT) :: lij
    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: unext
    REAL(KIND=8), DIMENSION(SIZE(unext))    :: Qminus, Pminus, Rminus
    REAL(KIND=8), PARAMETER :: smallminus = -1.d-15
    REAL(KIND=8) :: fij
    INTEGER      :: i, j, p 

    Qminus = 0.d0-unext
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
       !===Go back to first-order if water height too small
       IF (unext(i).LE.inputs%htiny) THEN
          Rminus(i) = 0.d0
       END IF
    END DO
    DO i = 1, SIZE(unext)
       DO p = mat%ia(i), mat%ia(i+1) - 1
          j = mat%ja(p)
          fij = mat%aa(p)
          IF (fij.GE.0.d0) THEN
             lij%aa(p) = Rminus(j)
          ELSE
             lij%aa(p) = Rminus(i)
          END IF
       END DO
    END DO
  END SUBROUTINE FCT_positivity

END MODULE fct
