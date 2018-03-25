INCLUDE 'mkl_pardiso.f90'
MODULE pardiso_solve
CONTAINS
  SUBROUTINE solve_pardiso(a,ia,ja,b,x,mnum_in,maxfct_in)
    USE mkl_pardiso
    IMPLICIT NONE
    TYPE pardiso_pt
       TYPE(MKL_PARDISO_HANDLE), ALLOCATABLE :: pt(:)
    END type pardiso_pt
    type(pardiso_pt), ALLOCATABLE, SAVE :: my_pardiso_pt(:)
    INTEGER, PARAMETER  :: dp = KIND(1.0D0)
    INTEGER,       INTENT(IN) :: ia( : )
    INTEGER,       INTENT(IN) :: ja( : )
    REAL(KIND=DP), INTENT(IN) :: a( : )
    REAL(KIND=DP), INTENT(INOUT) :: b( : )
    REAL(KIND=DP), INTENT(OUT):: x( : )
    INTEGER, INTENT(IN) :: mnum_in
    INTEGER, OPTIONAL   :: maxfct_in

    !.. Internal solver memory pointer 
    TYPE(MKL_PARDISO_HANDLE), ALLOCATABLE, SAVE  :: pt(:)
    INTEGER, ALLOCATABLE,                   SAVE :: iparm(:)
    INTEGER,                                SAVE :: maxfct
    LOGICAL :: once=.TRUE.
    INTEGER :: mnum, mtype, phase, n, nrhs, error, msglvl, nnz
    INTEGER :: i, idum(1)
    REAL(KIND=DP) :: ddum(1)

    IF (once) THEN
       once = .FALSE.
       IF (.NOT.PRESENT(maxfct_in)) THEN
          WRITE(*,*) 'BUG in SUBROUTINE, .NOT.PRESENT(maxfct_in)'
          STOP
       END IF
       IF (maxfct_in.LE.0) THEN
          maxfct = 1
       ELSE
          maxfct = maxfct_in
       END IF
       ALLOCATE(my_pardiso_pt(maxfct))
       ALLOCATE( iparm ( 64 ) )
       DO i = 1, 64
          iparm(i) = 0
       END DO
       iparm(1) = 1 ! no solver default
       iparm(2) = 2 ! fill-in reordering from METIS
       iparm(4) = 0 ! no iterative-direct algorithm
       iparm(5) = 0 ! no user fill-in reducing permutation
       iparm(6) = 0 ! =0 solution on the first n compoments of x
       iparm(8) = 9 ! numbers of iterative refinement steps
       iparm(10) = 13 ! perturbe the pivot elements with 1E-13
       iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
       iparm(13) = 0 ! maximum weighted matching algorithm is switched-off (default for symmetric).
       ! Try iparm(13) = 1 in case of inappropriate accuracy
       iparm(14) = 0 ! Output: number of perturbed pivots
       iparm(18) = -1 ! Output: number of nonzeros in the factor LU
       iparm(19) = -1 ! Output: Mflops for LU factorization
       iparm(20) = 0 ! Output: Numbers of CG Iterations
    END IF

    IF (ABS(mnum_in)> maxfct) THEN
       WRITE(*,*) ' BUG in solve_pardiso, ABS(mnum_in)> maxfct'
       STOP
    END IF

    error  = 0 ! initialize error flag
    msglvl = 0 ! 1 print statistical information
    nrhs = 1 
    mtype= 11            ! unsymmetric matrix

    IF (mnum_in<0) THEN !Factorization
       mnum = ABS(mnum_in)
       IF (.NOT.ALLOCATED(my_pardiso_pt(mnum)%pt)) ALLOCATE (my_pardiso_pt(mnum)%pt(64))
       DO i = 1, 64
          my_pardiso_pt(mnum)%pt(i)%DUMMY =  0 
       END DO

       phase = 11 ! only reordering and symbolic factorization
       n = SIZE(ia)-1   
       CALL pardiso(my_pardiso_pt(mnum)%pt, maxfct, mnum, mtype, phase, n, a, ia, ja, &
            idum, nrhs, iparm, msglvl, ddum, ddum, error)
       !WRITE(*,*) 'Reordering completed ... '
       IF (error /= 0) THEN
          WRITE(*,*) 'The following ERROR was detected: ', error
          STOP
       END IF
       !WRITE(*,*) 'Number of nonzeros in factors = ',iparm(18)
       !WRITE(*,*) 'Number of factorization MFLOPS = ',iparm(19)
       !.. Factorization.
       phase = 22 ! only factorization
       CALL pardiso (my_pardiso_pt(mnum)%pt, maxfct, mnum, mtype, phase, n, a, ia, ja, &
            idum, nrhs, iparm, msglvl, ddum, ddum, error)
       !WRITE(*,*) 'Factorization completed ... '
       IF (error /= 0) THEN
          WRITE(*,*) 'The following ERROR was detected: ', error
          STOP
       ENDIF
    END IF

    mnum = ABS(mnum_in)
    n = SIZE(ia)-1  
    !.. Back substitution and iterative refinement
    iparm(8) = 2 ! max numbers of iterative refinement steps
    phase = 33 ! only factorization
    CALL pardiso (my_pardiso_pt(mnum)%pt, maxfct, mnum, mtype, phase, n, a, ia, ja, &
         idum, nrhs, iparm, msglvl, b, x, error)
    !WRITE(*,*) 'Solve completed ... '
    IF (error /= 0) THEN
       WRITE(*,*) 'The following ERROR was detected: ', error
       STOP
    ENDIF

  END SUBROUTINE solve_pardiso
END MODULE pardiso_solve
