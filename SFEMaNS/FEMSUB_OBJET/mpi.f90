!
!Authors Jean-Luc Guermond, Copyrights 2005
!
SUBROUTINE MPI_INIT(code)
  IMPLICIT NONE
  INTEGER:: code
  code = 0
  RETURN
END SUBROUTINE MPI_INIT

SUBROUTINE MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)
   IMPLICIT NONE
   INTEGER:: MPI_COMM_WORLD, rang, code
   code = 0
   rang = 0
   MPI_COMM_WORLD = 0
   RETURN
END SUBROUTINE MPI_COMM_RANK

SUBROUTINE MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,code)
   IMPLICIT NONE
   INTEGER:: MPI_COMM_WORLD, nb_procs, code
   code = 0
   nb_procs = 1
   MPI_COMM_WORLD = 0
   RETURN
END SUBROUTINE MPI_COMM_SIZE

SUBROUTINE MPI_ALLREDUCE(norm_loc,norm_tot,index,MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, code)
   IMPLICIT NONE
   INTEGER:: MPI_COMM_WORLD, MPI_DOUBLE_PRECISION, MPI_SUM, code, index 
   REAL(KIND=8) :: norm_loc, norm_tot
   code = 0
   MPI_COMM_WORLD = 0
   norm_tot = norm_loc
   RETURN
END SUBROUTINE MPI_ALLREDUCE

SUBROUTINE MPI_SENDRECV(array_in, occur_in, MPI_DOUBLE_PRECISION_in, proc_in, etiquette_in, &
               array_out, occur_out, MPI_DOUBLE_PRECISION_out, proc_out, etiquette_out, &
               MPI_COMM_WORLD, statut, code)
   IMPLICIT NONE
   INTEGER:: MPI_COMM_WORLD, MPI_DOUBLE_PRECISION_in, code, proc_in, etiquette_in, occur_in
   INTEGER:: MPI_DOUBLE_PRECISION_out, proc_out, etiquette_out, occur_out
   INTEGER, DIMENSION(:)      :: statut
   REAL(KIND=8), DIMENSION(:) :: array_in, array_out
   code = 0
   MPI_COMM_WORLD = 0
   array_out = array_in
END SUBROUTINE MPI_SENDRECV

FUNCTION MPI_WTIME() RESULT(time)
   IMPLICIT NONE

   REAL(KIND=8) :: time
   INTEGER :: count, count_rate, count_max

   CALL SYSTEM_CLOCK(COUNT, COUNT_RATE, COUNT_MAX)
   time = (1.d0*count)/count_rate

END FUNCTION MPI_WTIME

SUBROUTINE MPI_FINALIZE(code)
   IMPLICIT NONE
   INTEGER:: code
   STOP
END SUBROUTINE MPI_FINALIZE

SUBROUTINE MPI_BARRIER(MPI_COMM_WORLD,code)
   IMPLICIT NONE
   INTEGER:: MPI_COMM_WORLD, code
   RETURN
END SUBROUTINE MPI_BARRIER

SUBROUTINE MPI_ALLTOALL (dist_field, longueur_tranche, MPID, combined_field, longueur_tranche_2, &
                    MPID_2, MPI_COMM_WORLD, code)
   IMPLICIT NONE
   INTEGER :: MPID, MPID_2, MPI_COMM_WORLD, code, longueur_tranche, longueur_tranche_2
   REAL(KIND=8), DIMENSION(longueur_tranche) :: dist_field, combined_field
   combined_field = dist_field
   code = 0
   MPI_COMM_WORLD = 0
   RETURN
END SUBROUTINE MPI_ALLTOALL
