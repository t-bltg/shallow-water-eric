PROGRAM test
   IMPLICIT NONE
   INCLUDE 'inc.h'
   REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: vin 
   ALLOCATE(vin(100))
   
   CALL mpialltoall(vin)
   STOP
END PROGRAM test 
