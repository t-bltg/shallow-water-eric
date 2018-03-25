PROGRAM test
    IMPLICIT NONE
    include 'mpif.h'
    REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE  :: V1, V2
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE    :: V1r, V2r
    INTEGER :: m_max_c, np, i, n, code
    REAL(KIND=8) :: t

    CALL MPI_INIT(code)
    !WRITE(*,*) ' M_max_c, np ', m_max_c, np
    !READ(*,*) m_max_c, np
    m_max_c = 24; np=10000

    ALLOCATE(V1r(m_max_c,np), V1(6,np,m_max_c))
    ALLOCATE(V2r(np, m_max_c), V2(np,m_max_c, 6))

    t = MPI_WTIME()
    DO i = 1, m_max_c
       V1r(i,:) = V1(1,:,i)
    END DO
    t = MPI_WTIME() -t
    WRITE(*,*) ' Temps 1', t

    t = MPI_WTIME()
    DO i = 1, m_max_c
       V2r(:,i) = V1(1,:,i)
    END DO
    t = MPI_WTIME() -t
    WRITE(*,*) ' Temps 2', t
    t = MPI_WTIME()
    DO i = 1, m_max_c
       V2r(:,i) = V2(:,i,1)
    END DO
    t = MPI_WTIME() -t
    WRITE(*,*) ' Temps 3', t

    t = MPI_WTIME()
    DO i =1 , m_max_c
       V1(1,:,i) = V1r(2*i-1,:)
    END DO
    t = MPI_WTIME() -t
    WRITE(*,*) ' Temps4 ', t

    t = MPI_WTIME()
    DO i =1 , m_max_c
       V2(:,i, 1) = V1r(2*i-1,:)
    END DO
    t = MPI_WTIME() -t
    WRITE(*,*) ' Temps5 ', t

    STOP
    CALL MPI_FINALIZE(code)
  END PROGRAM test
