INTERFACE
   SUBROUTINE mpialltoall(u, longueur_tranche, MPID, v)
      REAL(KIND=8), DIMENSION(:,:,:) :: u
      REAL(KIND=8), DIMENSION(:,:,:) :: v 
      INTEGER :: longueur_tranche, MPID
   END SUBROUTINE mpialltoall
END INTERFACE
