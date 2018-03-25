SUBROUTINE mpialltoall(u,longueur_tranche, MPID,v)
   IMPLICIT NONE
   REAL(KIND=8), DIMENSION(:,:,:) :: u
   REAL(KIND=8), DIMENSION(:,:,:) :: v
   INTEGER :: longueur_tranche, MPID
END SUBROUTINE mpialltoall
