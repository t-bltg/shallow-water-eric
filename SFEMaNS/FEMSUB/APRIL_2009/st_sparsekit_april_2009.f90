!
!Authors: Jean-Luc Guermond, Lugi Quartapelle, Copyright 1994
!
MODULE st_matrix

CONTAINS

  !=========================================================================

  SUBROUTINE st_sparsekit(jj,  ja, ia)
    ! VERY VERY BAD

    !  input coefficient structure of the matrix and
    !  perform the ordering of the unknowns

    !  jj(nodes_per_element, number_of_elements)
    !                  --->  node number in the grid

    IMPLICIT NONE

    INTEGER,      DIMENSION(:,:), INTENT(IN)    :: jj
    INTEGER,      DIMENSION(:),   INTENT(OUT)   :: ja
    INTEGER,      DIMENSION(:),   INTENT(OUT)   :: ia


    INTEGER :: me, np, nw
    INTEGER :: i, j, p, p0, l, m, ni, nj, nnz


    nw = SIZE(jj, 1)
    me = SIZE(jj, 2)
    np = SIZE(ia) - 1

    DO i = 1, np
       ja(i) = i
       ia(i) = i
    ENDDO
    nnz = np

    ia(np+1) = np + 1

    DO m = 1, me

       DO ni = 1, nw;  i = jj(ni, m)
          bc_nj:   DO nj = 1, nw;  j = jj(nj, m)

             !            IF (j > i  .OR.  .NOT.symm) THEN   !for non symmetric ordering only

             bc_p:          DO p = ia(i),  ia(i+1) - 1
                IF (ja(p) == j) CYCLE bc_nj
                IF (ja(p) > j ) THEN 
                   p0 = p; EXIT  bc_p
                ELSE 
                   p0 = p + 1
                END IF
             ENDDO bc_p

             DO p = ia(np+1) - 1,  p0,  -1
                ja(p+1) = ja(p)
             ENDDO

             ja(p0) = j

             DO l = i + 1,  np + 1
                ia(l) = ia(l) + 1
             ENDDO
             nnz = nnz + 1 

             !            ENDIF        !Comment to be removed if symmetric ordering chosen 

          ENDDO bc_nj
       ENDDO

    ENDDO

    IF (nnz.GT.SIZE(ja)) THEN
       WRITE(*,*) 'ST_SPARSEKIT: dimension de ja doit etre >= ',nnz
       STOP
    END IF

  END SUBROUTINE st_sparsekit

  !-------------------------------------------------------------------------------

  SUBROUTINE st_csr_p2(nw_c, jj, neigh, ja, ia)

    !  input coefficient structure of the matrix and
    !  perform the ordering of the unknowns

    !  jj(nodes_per_element, number_of_elements)
    !                  --->  node number in the grid

    IMPLICIT NONE

    INTEGER,                      INTENT(IN)    :: nw_c
    INTEGER,      DIMENSION(:,:), INTENT(IN)    :: jj, neigh

    INTEGER,      DIMENSION(:),   INTENT(OUT)   :: ja
    INTEGER,      DIMENSION(:),   INTENT(OUT)   :: ia

    INTEGER :: me, np, nw
    INTEGER :: m, n, jl, i, nja

    INTEGER, DIMENSION(:), ALLOCATABLE :: element_index
    LOGICAL, DIMENSION(:), ALLOCATABLE :: visited

    IF ((nw_c .NE. 3) .AND. (nw_c .NE. 4)) THEN
       WRITE(*,*) ' ST_CSR_P2: problem with nw_c = ', nw_c
       STOP
    END IF

    nw = SIZE(jj, 1)
    me = SIZE(jj, 2)
    np = SIZE(ia) - 1

    ALLOCATE ( element_index(np), visited(me) )

    visited = .FALSE.

    DO m = 1, me
       DO n = 1, nw
          element_index(jj(n,m)) = m
       ENDDO
    ENDDO

    jl = 0

    DO i = 1, np

       ia(i) = jl+1; m = element_index(i)

       nja = 0
       CALL build_bubble_p2 (nw_c, i, m, jl, jj, nja, ja, neigh, visited)
       !      CALL reset_bubble (m, next_index, visited)
       visited = .FALSE.

       jl = jl+nja
       IF (jl .GT. SIZE(ja)) THEN
          WRITE(*,*) 'ST_SPARSEKIT: dimension de ja doit etre >= ',jl
          STOP
       END IF

    ENDDO

    ia(np+1) = jl+1

    DEALLOCATE ( element_index, visited )

  END SUBROUTINE st_csr_p2

  SUBROUTINE st_csr_p1(jj, neigh, ja, ia)

    !  input coefficient structure of the matrix and
    !  perform the ordering of the unknowns

    !  jj(nodes_per_element, number_of_elements)
    !                  --->  node number in the grid

    IMPLICIT NONE

    INTEGER,      DIMENSION(:,:), INTENT(IN)    :: jj, neigh

    INTEGER,      DIMENSION(:),   INTENT(OUT)   :: ja
    INTEGER,      DIMENSION(:),   INTENT(OUT)   :: ia


    INTEGER :: me, np, nw
    INTEGER :: m, n, jl, i, nja

    INTEGER, DIMENSION(:), ALLOCATABLE :: element_index
    LOGICAL, DIMENSION(:), ALLOCATABLE :: visited


    nw = SIZE(jj, 1)
    me = SIZE(jj, 2)
    np = SIZE(ia) - 1

    ALLOCATE ( element_index(np), visited(me) )

    visited = .FALSE.

    DO m = 1, me
       DO n = 1, nw
          element_index(jj(n,m)) = m
       ENDDO
    ENDDO

    jl = 0

    DO i = 1, np

       ia(i) = jl+1; m = element_index(i)

       nja = 0
       CALL build_bubble_p1 (i, m, jl, jj, nja, ja, neigh, visited)
       !      CALL reset_bubble (m, next_index, visited)
       visited = .FALSE.

       jl = jl+nja

       IF (jl .GT. SIZE(ja)) THEN
          WRITE(*,*) 'ST_SPARSEKIT: dimension de ja doit etre >= ',jl
          STOP
       END IF

    ENDDO

    ia(np+1) = jl+1

    DEALLOCATE ( element_index, visited )

  END SUBROUTINE st_csr_p1

  SUBROUTINE zt_sparsekit(jj, neigh, ja, ia)

    !  input coefficient structure of the matrix and
    !  perform the ordering of the unknowns

    !  jj(nodes_per_element, number_of_elements)
    !                  --->  node number in the grid

    IMPLICIT NONE

    INTEGER,      DIMENSION(:,:), INTENT(IN)    :: jj, neigh

    INTEGER,      DIMENSION(:),   INTENT(OUT)   :: ja
    INTEGER,      DIMENSION(:),   INTENT(OUT)   :: ia


    INTEGER :: me, np, nw
    INTEGER :: m, n, jl, i, nja

    INTEGER, DIMENSION(:), ALLOCATABLE :: element_index
    LOGICAL, DIMENSION(:), ALLOCATABLE :: visited


    nw = SIZE(jj, 1)
    me = SIZE(jj, 2)
    np = SIZE(ia) - 1

    ALLOCATE ( element_index(np), visited(me) )

    visited = .FALSE.

    DO m = 1, me
       DO n = 1, nw
          element_index(jj(n,m)) = m
       ENDDO
    ENDDO

    jl = 0

    DO i = 1, np

       ia(i) = jl+1; m = element_index(i)

       nja = 0
       CALL build_bubble_p1 (i, m, jl, jj, nja, ja, neigh, visited)
       !      CALL reset_bubble (m, next_index, visited)
       visited = .FALSE.

       jl = jl+nja

       IF (jl .GT. SIZE(ja)) THEN
          WRITE(*,*) 'ST_SPARSEKIT: dimension de ja doit etre >= ',jl
          STOP
       END IF

    ENDDO

    ia(np+1) = jl+1

    DEALLOCATE ( element_index, visited )

  END SUBROUTINE zt_sparsekit

  !=========================================================================

  RECURSIVE SUBROUTINE build_bubble_p1 (p_index, e_index, jl, jj, nja, ja, neigh, visited)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: p_index, e_index, jl
    INTEGER, INTENT(INOUT) :: nja ! a initialiser avant le premier appel
    INTEGER, DIMENSION(:,:), INTENT(IN) :: jj, neigh
    INTEGER, DIMENSION(:), INTENT(INOUT) :: ja
    LOGICAL, DIMENSION(:), INTENT(INOUT) :: visited

    INTEGER :: nw, i, j, neighbour_index

    nw = SIZE(jj,1)

    DO i=1,nw
       CALL insert_node (jj(i,e_index), jl, nja, ja)
    ENDDO

    visited(e_index) = .TRUE.

    DO i=1,nw
       neighbour_index = neigh(i,e_index)
       IF (neighbour_index .GT. 0) THEN
          IF ( .NOT. visited(neighbour_index) ) THEN
             points_loop:  DO j=1,nw
                IF (jj(j,neighbour_index) .EQ. p_index) THEN
                   CALL build_bubble_p1 (p_index, neighbour_index, jl, jj, nja, ja, neigh, visited)
                   EXIT points_loop
                ENDIF
             ENDDO points_loop
          ENDIF
       ENDIF
    ENDDO

  END SUBROUTINE build_bubble_p1

  !=========================================================================

  RECURSIVE SUBROUTINE build_bubble_p2(nw_c, p_index, e_index, jl, jj, nja, ja, &
       neigh, visited)

    IMPLICIT NONE

    INTEGER, INTENT(IN)  :: nw_c, p_index, e_index, jl
    INTEGER, INTENT(INOUT) :: nja ! a initialiser avant le premier appel
    INTEGER, DIMENSION(:,:), INTENT(IN) :: jj, neigh
    INTEGER, DIMENSION(:), INTENT(INOUT) :: ja
    LOGICAL, DIMENSION(:), INTENT(INOUT) :: visited

    INTEGER :: nw, i, j, neighbour_index

    nw = SIZE(jj,1)

    DO i=1,nw
       CALL insert_node (jj(i,e_index), jl, nja, ja)
    ENDDO

    visited(e_index) = .TRUE.

    DO i=1,nw_c
       neighbour_index = neigh(i,e_index)
       IF (neighbour_index .GT. 0) THEN
          IF ( .NOT. visited(neighbour_index) ) THEN
             points_loop: DO j=1,nw
                IF (jj(j,neighbour_index) .EQ. p_index) THEN
                   CALL build_bubble_p2 (nw_c, p_index, neighbour_index, jl, jj, nja, ja, neigh, visited)
                   EXIT points_loop
                ENDIF
             END DO points_loop
          END IF
       END IF
    ENDDO

  END SUBROUTINE build_bubble_p2

  !=========================================================================

  SUBROUTINE reset_bubble (e_index, next_index, visited)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: e_index
    INTEGER, DIMENSION(:), INTENT(INOUT) :: next_index
    LOGICAL, DIMENSION(:), INTENT(INOUT) :: visited

    INTEGER :: i, j

    i = e_index

    DO
       visited(i) = .FALSE.
       j = next_index(i)
       next_index(i) = 0
       i = j
       IF ( i .EQ. 0 ) EXIT
    ENDDO

  END SUBROUTINE reset_bubble

  !=========================================================================

  SUBROUTINE insert_node (j, jl, nja, ja)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: j, jl

    INTEGER, INTENT(INOUT) :: nja
    INTEGER, DIMENSION(:), INTENT(INOUT) :: ja

    LOGICAL :: addj
    INTEGER :: n, m

    addj = .TRUE.

    DO n = jl+1, jl+nja
       IF ( j .LE. ja(n) ) THEN
          IF ( j .EQ. ja(n) ) THEN
             addj = .FALSE.
             EXIT
          ELSE
             nja = nja+1
             DO m = jl+nja, n+1, -1
                ja(m) = ja(m-1)
             ENDDO
             ja(n) = j
             addj = .FALSE.
             EXIT
          ENDIF
       ENDIF
    ENDDO

    IF ( addj ) THEN
       nja = nja+1
       ja(jl+nja) = j
    ENDIF

  END SUBROUTINE insert_node

  !=========================================================================

  !TESTTTTTTTTT
  SUBROUTINE st_sparsekit_test(jj,  ja, ia)

    !  input coefficient structure of the matrix and
    !  perform the ordering of the unknowns

    !  jj(nodes_per_element, number_of_elements)
    !                  --->  node number in the grid

    IMPLICIT NONE

    INTEGER,      DIMENSION(:,:), INTENT(IN)    :: jj

    INTEGER,      DIMENSION(:),   INTENT(OUT)   :: ja
    INTEGER,      DIMENSION(:),   INTENT(OUT)   :: ia


    INTEGER :: me, np, nw
    INTEGER :: i, j, p, p0, l, m, ni, nj, nnz


    nw = SIZE(jj, 1)
    me = SIZE(jj, 2)
    np = SIZE(ia) - 1

    DO i = 1, np
       ja(i) = i
       ia(i) = i
    ENDDO
    nnz = np

    ia(np+1) = np + 1

    DO m = 1, me
       PRINT*, ' m = ', m

       DO ni = 1, nw;  i = jj(ni, m)
          bc_nj:   DO nj = 1, nw;  j = jj(nj, m)

             !            IF (j > i  .OR.  .NOT.symm) THEN   !for non symmetric ordering only

             bc_p:          DO p = ia(i),  ia(i+1) - 1
                IF (ja(p) == j) CYCLE bc_nj
                IF (ja(p) > j ) THEN 
                   p0 = p; EXIT  bc_p
                ELSE 
                   p0 = p + 1
                END IF
             ENDDO bc_p

             DO p = ia(np+1) - 1,  p0,  -1
                ja(p+1) = ja(p)
             ENDDO

             ja(p0) = j

             DO l = i + 1,  np + 1
                ia(l) = ia(l) + 1
             ENDDO
             nnz = nnz + 1 

             !            ENDIF        !Comment to be removed if symmetric ordering chosen 

          ENDDO bc_nj
       ENDDO

    ENDDO

    IF (nnz.GT.SIZE(ja)) THEN
       WRITE(*,*) 'ST_SPARSEKIT: dimension de ja doit etre >= ',nnz
       STOP
    END IF

  END SUBROUTINE st_sparsekit_test

  !=========================================================================

  SUBROUTINE st_csr_bloc(ia,ja,ia_b,ja_b,n_b)
    !
    ! Builds the CSR structure of a 3D vector matrix
    ! from its scalar counterpart.
    !
    IMPLICIT NONE

    INTEGER, DIMENSION(:), INTENT(IN)  :: ia, ja
    INTEGER, DIMENSION(:), POINTER     :: ia_b, ja_b
    INTEGER,               INTENT(IN)  :: n_b  ! Number of blocs

    INTEGER :: ki, kj, i, ib, jb, p, pb, bloc_size, np, nnz

    np = SIZE(ia) - 1
    nnz = ia(np+1) - ia(1)
    ALLOCATE (ia_b(n_b*np+1), ja_b(n_b*n_b*nnz))

    bloc_size = SIZE(ia) - 1

    ia_b(1) = ia(1) 

    DO ki = 1, n_b
       DO i = 2, bloc_size + 1 
          ib = i + (ki-1)*bloc_size  
          ia_b(ib) = ia_b(ib-1) + n_b*(ia(i) - ia(i-1)) 
       END DO
    END DO

    DO ki = 1, n_b
       DO i = 1, bloc_size
          ib = i + (ki-1)*bloc_size

          DO kj = 1,  n_b
             DO p = ia(i), ia(i+1) - 1          
                jb = ja(p) + (kj-1)*bloc_size 

                pb = ia_b(ib)  +  p - ia(i)  +  (kj-1)*(ia(i+1)-ia(i)) 
                ja_b(pb) = jb 

             END DO
          END DO

       END DO
    END DO

  END SUBROUTINE st_csr_bloc

  !=========================================================================

  SUBROUTINE st_csr_p2_opt(np, jj, neigh, ia, ja)

    !  input coefficient structure of the matrix and
    !  perform the ordering of the unknowns

    !  jj(nodes_per_element, number_of_elements)
    !                  --->  node number in the grid

    IMPLICIT NONE

    INTEGER,                 INTENT(IN)  :: np
    INTEGER, DIMENSION(:,:), INTENT(IN)  :: jj, neigh
    INTEGER, DIMENSION(:),   POINTER     :: ia, ja

    INTEGER :: me, nw, nmax, nw_c
    INTEGER :: m, n, jl, i, nja

    INTEGER, DIMENSION(:), ALLOCATABLE :: ja_work
    INTEGER, DIMENSION(:), ALLOCATABLE :: element_index
    LOGICAL, DIMENSION(:), ALLOCATABLE :: visited

    nw_c = SIZE(neigh,1)

    IF ((nw_c .NE. 3) .AND. (nw_c .NE. 4)) THEN
       WRITE(*,*) ' ST_CSR_P2: problem with nw_c = ', nw_c
       STOP
    END IF

    nw = SIZE(jj, 1)
    me = SIZE(jj, 2)

    ALLOCATE ( element_index(np), visited(me) , ja_work(60*np), ia(np+1))

    visited = .FALSE.

    DO m = 1, me
       DO n = 1, nw
          element_index(jj(n,m)) = m
       ENDDO
    ENDDO

    jl = 0

    DO i = 1, np

       ia(i) = jl+1; m = element_index(i)

       nja = 0
       CALL build_bubble_p2 (nw_c, i, m, jl, jj, nja, ja_work, neigh, visited)
       visited = .FALSE.

       jl = jl+nja
       IF (jl .GT. SIZE(ja_work)) THEN
          WRITE(*,*) 'ST_SPARSEKIT: dimension de ja doit etre >= ',jl
          STOP
       END IF

    ENDDO

    ia(np+1) = jl+1
    nmax = ia(np+1) -ia(1)
    ALLOCATE(ja(nmax))
    ja = ja_work(1:nmax)

    DEALLOCATE ( element_index, visited, ja_work )

  END SUBROUTINE st_csr_p2_opt

  !=========================================================================

  SUBROUTINE st_csr_p1_opt(np, jj, neigh, ia, ja)

    !  input coefficient structure of the matrix and
    !  perform the ordering of the unknowns

    !  jj(nodes_per_element, number_of_elements)
    !                  --->  node number in the grid

    IMPLICIT NONE

    INTEGER, DIMENSION(:,:), INTENT(IN)  :: jj, neigh
    INTEGER,                 INTENT(IN)  :: np
    INTEGER, DIMENSION(:),   POINTER     :: ia, ja

    INTEGER :: me, nw, nmax
    INTEGER :: m, n, jl, i, nja

    INTEGER, DIMENSION(:), ALLOCATABLE :: ja_work, element_index
    LOGICAL, DIMENSION(:), ALLOCATABLE :: visited

    nw = SIZE(jj, 1)
    me = SIZE(jj, 2)

    ALLOCATE ( element_index(np), visited(me), ja_work(60*np), ia(np+1))

    visited = .FALSE.

    DO m = 1, me
       DO n = 1, nw
          element_index(jj(n,m)) = m
       ENDDO
    ENDDO

    jl = 0

    DO i = 1, np

       ia(i) = jl+1; m = element_index(i)

       nja = 0
       CALL build_bubble_p1 (i, m, jl, jj, nja, ja_work, neigh, visited)
       visited = .FALSE.

       jl = jl+nja

       IF (jl .GT. SIZE(ja_work)) THEN
          WRITE(*,*) 'ST_SPARSEKIT: dimension de ja doit etre >= ',jl
          STOP
       END IF

    ENDDO

    ia(np+1) = jl+1
    nmax = ia(np+1)-ia(1)
    ALLOCATE(ja(nmax))
    ja = ja_work(1:nmax)

    DEALLOCATE ( element_index, visited, ja_work )

  END SUBROUTINE st_csr_p1_opt

  !=========================================================================

  SUBROUTINE st_csr_sgs_p1(np, jj, neigh, ia, ja)

    IMPLICIT NONE

    INTEGER,                 INTENT(IN)  :: np
    INTEGER, DIMENSION(:,:), INTENT(IN)  :: jj, neigh
    INTEGER, DIMENSION(:),   POINTER     :: ia, ja

    INTEGER :: index, i, j, k, p, pp, jl, nja, nmax, lmax=100

    INTEGER, DIMENSION(:), ALLOCATABLE :: ia_work, ja_work, liste

    CALL st_csr_p1_opt(np, jj, neigh, ia, ja)

    ALLOCATE (ja_work(3*SIZE(ja)), liste(lmax), ia_work(np+1))

    ia_work(1) = 1

    DO i = 1, np 

       index = 0
       DO p = ia(i), ia(i+1) - 1
          j = ja(p)
          DO pp = ia(j), ia(j+1) - 1
             index = index + 1
             liste(index) = ja(pp)
          END DO
       END DO
       IF (index .GE. lmax) THEN 
          WRITE(*,*) 'lmax pas assez grand dans st_csr_sgs_p1'
          STOP
       END IF

       nja = ia(i+1) - ia(i)
       jl  = ia_work(i) - 1
       ja_work(jl+1:jl+nja) = ja(ia(i):ia(i+1)-1)
       DO k = 1, index
          CALL insert_node(liste(k),jl,nja,ja_work) 
       ENDDO
       ia_work(i+1) = ia_work(i) + nja

    ENDDO

    ia = ia_work
    nmax = ia(np+1)-ia(1)
    NULLIFY(ja)
    ALLOCATE(ja(nmax))
    ja = ja_work(1:nmax)

    DEALLOCATE (liste, ia_work, ja_work)

  END SUBROUTINE st_csr_sgs_p1

  !=========================================================================

  SUBROUTINE st_csr_sgs_p2(np, jj, neigh, ia, ja)

    IMPLICIT NONE

    INTEGER,                 INTENT(IN)  :: np
    INTEGER, DIMENSION(:,:), INTENT(IN)  :: jj, neigh
    INTEGER, DIMENSION(:),   POINTER     :: ia, ja

    INTEGER :: index, i, j, k, p, pp, jl, nja, nmax, lmax=500

    INTEGER, DIMENSION(:), ALLOCATABLE :: ia_work, ja_work, liste

    CALL st_csr_p2_opt(np, jj, neigh, ia, ja)

    ALLOCATE (ja_work(8*SIZE(ja)), liste(lmax), ia_work(np+1))

    ia_work(1) = 1

    DO i = 1, np

       index = 0
       DO p = ia(i), ia(i+1) - 1
          j = ja(p)
          DO pp = ia(j), ia(j+1) - 1
             index = index + 1
             liste(index) = ja(pp)
          END DO
       END DO
       IF (index .GE. lmax) THEN
          WRITE(*,*) 'lmax pas assez grand dans st_csr_sgs_p1', index
          STOP
       END IF

       nja = ia(i+1) - ia(i)
       jl  = ia_work(i) - 1
       ja_work(jl+1:jl+nja) = ja(ia(i):ia(i+1)-1)
       DO k = 1, index
          CALL insert_node(liste(k),jl,nja,ja_work)
       ENDDO
       ia_work(i+1) = ia_work(i) + nja

    ENDDO

    ia = ia_work
    nmax = ia(np+1)-ia(1)
    NULLIFY(ja)
    ALLOCATE(ja(nmax))
    ja = ja_work(1:nmax)

    DEALLOCATE (liste, ia_work, ja_work)

  END SUBROUTINE st_csr_sgs_p2


  SUBROUTINE st_csr(jj, ia, ja)

    !  input coefficient structure of the matrix and
    !  perform the ordering of the unknowns

    !  jj(nodes_per_element, number_of_elements)
    !                  --->  node number in the grid

    IMPLICIT NONE

    INTEGER, DIMENSION(:,:), INTENT(IN)  :: jj
    INTEGER, DIMENSION(:),   POINTER     :: ia, ja

    INTEGER :: nparm=90
    INTEGER :: me, nw, nmax, np
    INTEGER :: m, ni, nj, i, j, n_a_d

    INTEGER, DIMENSION(:,:), ALLOCATABLE :: ja_work
    INTEGER, DIMENSION(:),   ALLOCATABLE :: nja, a_d

    nw = SIZE(jj, 1)
    me = SIZE(jj, 2)
    np = MAXVAL(jj)

    ALLOCATE (ja_work(np,nparm), ia(np+1), a_d(nparm), nja(np))

    ja_work = 0
    nja = 1
    DO i = 1, np
       ja_work(i,1) = i
    END DO

    DO m = 1, me
       DO ni = 1, nw
          i = jj(ni,m)

          DO nj = 1, nw
             j = jj(nj,m)
             IF (MINVAL(ABS(ja_work(i,1:nja(i))-j)) /= 0) THEN 
                nja(i) = nja(i) + 1
                ja_work(i,nja(i)) = j 
             END IF
          END DO


       END DO
    END DO

    IF (MAXVAL(nja)>nparm) THEN 
       WRITE(*,*) 'ST_SPARSEKIT: dimension de ja doit etre >= ',nparm
       STOP
    END IF


    nmax = 0
    DO i = 1, np
       nmax = nmax + nja(i)
    END DO
    ALLOCATE(ja(nmax))

    ia(1) = 1
    DO i = 1, np
       CALL tri_jlg (ja_work(i,1:nja(i)), a_d, n_a_d)
       IF (n_a_d /= nja(i)) THEN
          WRITE(*,*) ' BUG : st_p1_CSR'
          WRITE(*,*) 'n_a_d ', n_a_d, 'nja(i)', nja(i) 
          STOP
       END IF
       ia(i+1) = ia(i) + nja(i)
       ja(ia(i):ia(i+1)-1) = a_d(1:nja(i))
    END DO


    DEALLOCATE ( ja_work, nja, a_d )

  END SUBROUTINE st_csr

  SUBROUTINE st_p1_csr(jj, ia, ja)

    !  input coefficient structure of the matrix and
    !  perform the ordering of the unknowns

    !  jj(nodes_per_element, number_of_elements)
    !                  --->  node number in the grid

    IMPLICIT NONE

    INTEGER, DIMENSION(:,:), INTENT(IN)  :: jj
    INTEGER, DIMENSION(:),   POINTER     :: ia, ja

    INTEGER :: nparm=90
    INTEGER :: me, nw, nmax, np
    INTEGER :: m, ni, nj, i, j, n_a_d

    INTEGER, DIMENSION(:,:), ALLOCATABLE :: ja_work
    INTEGER, DIMENSION(:),   ALLOCATABLE :: nja, a_d

    nw = SIZE(jj, 1)
    me = SIZE(jj, 2)
    np = MAXVAL(jj)

    ALLOCATE (ja_work(np,nparm), ia(np+1), a_d(nparm), nja(np))

    ja_work = 0
    nja = 1
    DO i = 1, np
       ja_work(i,1) = i
    END DO

    DO m = 1, me
       DO ni = 1, nw
          i = jj(ni,m)

          DO nj = 1, nw
             j = jj(nj,m)
             IF (MINVAL(ABS(ja_work(i,1:nja(i))-j)) /= 0) THEN 
                nja(i) = nja(i) + 1
                ja_work(i,nja(i)) = j 
             END IF
          END DO


       END DO
    END DO

    IF (MAXVAL(nja)>nparm) THEN 
       WRITE(*,*) 'ST_SPARSEKIT: dimension de ja doit etre >= ',nparm
       STOP
    END IF


    nmax = 0
    DO i = 1, np
       nmax = nmax + nja(i)
    END DO
    ALLOCATE(ja(nmax))

    ia(1) = 1
    DO i = 1, np
       CALL tri_jlg (ja_work(i,1:nja(i)), a_d, n_a_d)
       IF (n_a_d /= nja(i)) THEN
          WRITE(*,*) ' BUG : st_p1_CSR'
          WRITE(*,*) 'n_a_d ', n_a_d, 'nja(i)', nja(i) 
          STOP
       END IF
       ia(i+1) = ia(i) + nja(i)
       ja(ia(i):ia(i+1)-1) = a_d(1:nja(i))
    END DO


    DEALLOCATE ( ja_work, nja, a_d )

  END SUBROUTINE st_p1_csr

  SUBROUTINE tri_jlg (a,  a_d, n_a_d)

    !  sort in ascending order of the integer array  a  and generation
    !  of the integer array  a_d  whose first  n_a_d  leading entries
    !  contain different values in ascending order, while all the
    !  remaining entries are set to zero

    !  sorting by Shell's method.

    IMPLICIT NONE

    INTEGER, DIMENSION(:), INTENT(INOUT) :: a
    INTEGER, DIMENSION(:), INTENT(OUT)   :: a_d
    INTEGER,               INTENT(OUT)   :: n_a_d

    INTEGER :: n, na, inc, i, j, k, ia

    na = SIZE(a)

    !  sort phase

    IF (na == 0) THEN
       n_a_d = 0
       RETURN
    ENDIF

    inc = 1
    DO WHILE (inc <= na)
       inc = inc * 3
       inc = inc + 1
    ENDDO

    DO WHILE (inc > 1)
       inc = inc/3
       DO i = inc + 1, na
          ia = a(i)
          j = i
          DO WHILE (a(j-inc) > ia)
             a(j) = a(j-inc)
             j = j - inc
             IF (j <= inc) EXIT
          ENDDO
          a(j) = ia
       ENDDO
    ENDDO

    !  compression phase

    n = 1
    a_d(n) = a(1)
    DO k = 2, na
       IF (a(k) > a(k-1)) THEN
          n = n + 1
          a_d(n) = a(k)
          !TEST JUIN 13 2008
       ELSE
          WRITE(*,*) 'We have a problem in the compression phase of tri_jlg', k, k-1
          !TEST JUIN 13 2008
       ENDIF
    ENDDO

    n_a_d = n

    a_d(n_a_d + 1 : na) = 0

  END SUBROUTINE tri_jlg


  SUBROUTINE st_mhd_csr(jj, jj_sub_c, i_d, list_dom, ia, ja)

    !  input coefficient structure of the matrix and
    !  perform the ordering of the unknowns

    !  jj(nodes_per_element, number_of_elements)
    !                  --->  node number in the grid

    IMPLICIT NONE

    INTEGER, DIMENSION(:,:), INTENT(IN)  :: jj, jj_sub_c
    INTEGER, DIMENSION(:),   INTENT(IN)  :: i_d, list_dom
    INTEGER, DIMENSION(:),   POINTER     :: ia, ja

    INTEGER :: nparm=200
    INTEGER :: me, nw, nmax, np, me_sub_c, nw_sub_c, np_sub_c, np_tot
    INTEGER :: m, ni, nj, i, j, n_a_d, k, h

    INTEGER, DIMENSION(:,:), ALLOCATABLE :: ja_work
    INTEGER, DIMENSION(:),   ALLOCATABLE :: nja, a_d
    LOGICAL :: test


    nw = SIZE(jj, 1)
    me = SIZE(jj, 2)
    np = MAXVAL(jj)

    nw_sub_c = SIZE(jj_sub_c, 1)
    me_sub_c = SIZE(jj_sub_c, 2)
    IF (me_sub_c/=0) THEN 
       np_sub_c = MAXVAL(jj_sub_c)
    ELSE
       np_sub_c = 0
    END IF

    np_tot = 2*np + np_sub_c

    ALLOCATE (ja_work(np_tot,nparm), ia(np_tot+1), a_d(nparm), nja(np_tot))

    ja_work = 0
    nja = 1
    DO i = 1, np_tot
       ja_work(i,1) = i
    END DO

    ! Boucle sur le deuxieme bloc
    DO m = 1, me
       test = .TRUE. ! In conductor
       IF (MINVAL(ABS(i_d(m)-list_dom)) ==0) test = .FALSE. ! In void

       DO ni = 1, nw
          DO k = 1, 2
             i = jj(ni,m) + (k-1)*np

             DO nj = 1, nw
                DO h = 1, 2
                   j = jj(nj,m) + (h-1)*np
                   IF (MINVAL(ABS(ja_work(i,1:nja(i))-j)) /= 0) THEN 
                      nja(i)    = nja(i) + 1
                      ja_work(i,nja(i)) = j
                   END IF
                END DO
             END DO

             IF (test) CYCLE

             DO nj = 1, nw_sub_c
                j = jj_sub_c(nj,m) + 2*np
                IF (MINVAL(ABS(ja_work(i,1:nja(i))-j)) /= 0) THEN 
                   nja(i) = nja(i) + 1
                   ja_work(i,nja(i)) = j
                END IF
             END DO

          END DO
       END DO

       ! Boucle sur le deuxieme bloc
       IF (test) CYCLE

       DO ni = 1, nw_sub_c
          i = jj_sub_c(ni,m) + 2*np

          DO nj = 1, nw
             j = jj(nj,m)
             IF (MINVAL(ABS(ja_work(i,1:nja(i))-j)) /= 0) THEN 
                nja(i) = nja(i) + 2
                ja_work(i,nja(i)-1) = j
                ja_work(i,nja(i))   = j + np
             END IF
          END DO

          DO nj = 1, nw_sub_c
             j = jj_sub_c(nj,m) + 2*np
             IF (MINVAL(ABS(ja_work(i,1:nja(i))-j)) /= 0) THEN 
                nja(i) = nja(i) + 1
                ja_work(i,nja(i)) = j
             END IF
          END DO

       END DO
    END DO


    IF (MAXVAL(nja)>nparm) THEN 
       WRITE(*,*) 'ST_SPARSEKIT: dimension de ja doit etre >= ',nparm
       STOP
    END IF

    nmax = 0
    DO i = 1, np_tot
       nmax = nmax + nja(i)
    END DO
    ALLOCATE(ja(nmax))
    ia(1) = 1
    DO i = 1, np_tot
       CALL tri_jlg (ja_work(i,1:nja(i)), a_d, n_a_d)
       IF (n_a_d /= nja(i)) THEN
          WRITE(*,*) ' BUG : st_p1_CSR'
          WRITE(*,*) 'n_a_d ', n_a_d, 'nja(i)', nja(i) 
          STOP
       END IF
       ia(i+1) = ia(i) + nja(i)
       ja(ia(i):ia(i+1)-1) = a_d(1:nja(i))
    END DO


    DEALLOCATE ( ja_work, nja, a_d )

    DO m = 1, me
       DO ni = 1, nw
          DO k= 1, 2; i = jj(ni,m) + (k-1)*np
             DO nj = 1, nw
                DO h= 1, 2; j = jj(nj,m) + (h-1)*np
                   IF(MINVAL(ABS(ja(ia(i):ia(i+1)-1)-j)) /=0) THEN
                      WRITE(*,*) ' BUG in st_mhd_csr'
                      WRITE(*,*) 'i ', i, ' j ', j
                      WRITE(*,*) ja(ia(i):  ia(i+1) - 1)
                      STOP
                   END IF
                END DO
             END DO
          END DO
       END DO
    END DO

  END SUBROUTINE st_mhd_csr

  SUBROUTINE st_mhd_csr_te(jj, jj_sub_c, i_d, list_dom, ia, ja)

    !  input coefficient structure of the matrix and
    !  perform the ordering of the unknowns

    !  jj(nodes_per_element, number_of_elements)
    !                  --->  node number in the grid

    IMPLICIT NONE

    INTEGER, DIMENSION(:,:), INTENT(IN)  :: jj, jj_sub_c
    INTEGER, DIMENSION(:),   INTENT(IN)  :: i_d, list_dom
    INTEGER, DIMENSION(:),   POINTER     :: ia, ja

    INTEGER :: nparm=100
    INTEGER :: me, nw, nmax, np, me_sub_c, nw_sub_c, np_sub_c, np_tot
    INTEGER :: m, ni, nj, i, j, k, h, n_a_d, i_b, j_b, kd, bloc_size
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: ja_work
    INTEGER, DIMENSION(:),   ALLOCATABLE :: nja, a_d
    LOGICAL :: test

    nw = SIZE(jj, 1)
    me = SIZE(jj, 2)
    np = MAXVAL(jj)
    kd = 2

    nw_sub_c = SIZE(jj_sub_c, 1)
    me_sub_c = SIZE(jj_sub_c, 2)
    np_sub_c = MAXVAL(jj_sub_c)
    bloc_size = np_sub_c

    np_tot = np + 2*np_sub_c

    ALLOCATE (ja_work(np_tot,nparm), ia(np_tot+1), a_d(nparm), nja(np_tot))

    ja_work = 0
    nja = 1
    DO i = 1, np_tot
       ja_work(i,1) = i
    END DO


    DO m = 1, me

       test = .TRUE. ! In conductor
       IF (MINVAL(ABS(i_d(m)-list_dom))==0) test = .FALSE. ! In void

       DO ni = 1, nw;  i_b = jj(ni,m)
          i = i_b

          DO nj = 1, nw;  j_b = jj(nj,m)
             j = j_b 
             !  bloc 1
             IF (MINVAL(ABS(ja_work(i,1:nja(i))-j)) /= 0) THEN 
                nja(i)    = nja(i) + 1
                ja_work(i,nja(i)) = j
             END IF
          ENDDO

          IF (test) CYCLE
          DO nj = 1, nw_sub_c;  j_b = jj_sub_c(nj,m)   
             DO h = 1, kd
                j = j_b + (h-1)*bloc_size + np
                !  bloc 2
                IF (MINVAL(ABS(ja_work(i,1:nja(i))-j)) /= 0) THEN 
                   nja(i)    = nja(i) + 1
                   ja_work(i,nja(i)) = j
                END IF
             ENDDO
          ENDDO

       ENDDO

       IF (test) CYCLE

       DO ni = 1, nw_sub_c;  i_b = jj_sub_c(ni, m)
          DO k = 1, kd
             i = i_b + (k-1)*bloc_size + np

             DO nj = 1, nw;  j_b = jj(nj,m)
                j = j_b
                !  bloc 3
                IF (MINVAL(ABS(ja_work(i,1:nja(i))-j)) /= 0) THEN 
                   nja(i)    = nja(i) + 1
                   ja_work(i,nja(i)) = j
                END IF
             ENDDO

             DO nj = 1, nw_sub_c;  j_b = jj_sub_c(nj,m)
                DO h = 1, kd; 
                   j = j_b + (h-1)*bloc_size + np
                   !  bloc 4
                   IF (MINVAL(ABS(ja_work(i,1:nja(i))-j)) /= 0) THEN 
                      nja(i)    = nja(i) + 1
                      ja_work(i,nja(i)) = j
                   END IF
                ENDDO
             ENDDO


          ENDDO
       ENDDO

    ENDDO !fin boucle m


    IF (MAXVAL(nja)>nparm) THEN 
       WRITE(*,*) 'ST_SPARSEKIT: dimension de ja doit etre >= ',nparm
       STOP
    END IF

    nmax = 0
    DO i = 1, np_tot
       nmax = nmax + nja(i)
    END DO
    ALLOCATE(ja(nmax))
    ia(1) = 1
    DO i = 1, np_tot
       CALL tri_jlg (ja_work(i,1:nja(i)), a_d, n_a_d)
       IF (n_a_d /= nja(i)) THEN
          WRITE(*,*) ' BUG : st_CSR_mhd_te'
          WRITE(*,*) 'n_a_d ', n_a_d, 'nja(i)', nja(i) 
          STOP
       END IF
       ia(i+1) = ia(i) + nja(i)
       ja(ia(i):ia(i+1)-1) = a_d(1:nja(i))
    END DO

    DEALLOCATE ( ja_work, nja, a_d )

  END SUBROUTINE st_mhd_csr_te

  !-----------------------------------------------------------------------

  SUBROUTINE st_scr_maxwell_tm(H_mesh, phi_mesh, INTERFACE, ia_glob, ja_glob)

    USE def_type_mesh

    IMPLICIT NONE

    TYPE(mesh_type),           INTENT(IN) :: H_mesh, phi_mesh
    TYPE(mesh_type_interface), INTENT(IN) :: INTERFACE
    INTEGER, DIMENSION(:),     POINTER    :: ia_glob, ja_glob

    INTEGER, DIMENSION(:),     POINTER    :: ia_m, ja_m, ia_s, ja_s, ia_b, &
         ja_b, nja_glob, nja_b, a_d
    INTEGER, DIMENSION(:,:),   POINTER    :: ja_work
    INTEGER, PARAMETER :: param=90
    INTEGER :: np_m, np_s, np_b, np_glob, i_glob, j_glob, i, i_b, i_s, &
         nnz, p_i, p_f, nw, m, n, nn, n_a_d

    CALL st_csr(H_mesh%jj, ia_m, ja_m)
    np_m = SIZE(ia_m)-1
    CALL st_csr(phi_mesh%jj, ia_s, ja_s)
    CALL st_csr_bloc(ia_m, ja_m, ia_b, ja_b, 2)
    NULLIFY(ia_m,ja_m)

    np_b = SIZE(ia_b)-1
    np_s = SIZE(ia_s)-1
    np_glob = np_b + np_s

    ALLOCATE(ja_work(np_glob,param),a_d(param),nja_b(np_b),nja_glob(np_glob))

    nja_glob = 0
    DO i_b = 1, np_b
       nnz = ia_b(i_b+1) - ia_b(i_b)
       p_i = ia_b(i_b)
       p_f = ia_b(i_b+1) - 1
       ja_work(i_b,1:nnz) = ja_b(p_i:p_f)
       nja_b(i_b) = nnz
    END DO
    nja_glob(1:np_b) = nja_b

    nw = SIZE(interface%master_node,1)
    DO m = 1, interface%me
       DO n = 1, nw
          IF (interface%master_node(n,m) < 0) CYCLE
          i = interface%master_node(n,m)
          DO nn = 1, nw
             j_glob = phi_mesh%jj(nn,interface%slave_elem(m)) + np_b
             IF (MINVAL(ABS(ja_work(i,nja_b(i):nja_glob(i))-j_glob)) /= 0) THEN 
                nja_glob(i) = nja_glob(i) + 1
                nja_glob(i+np_m) = nja_glob(i+np_m) + 1
                ja_work(i,nja_glob(i)) = j_glob
                ja_work(i+np_m,nja_glob(i+np_m)) = j_glob
                nja_glob(j_glob) = nja_glob(j_glob) + 2
                ja_work(j_glob,nja_glob(j_glob)-1) = i
                ja_work(j_glob,nja_glob(j_glob)) = i + np_m
             END IF
          END DO
       END DO
    END DO

    DO i_s = 1, np_s
       i_glob =  i_s + np_b
       nnz = ia_s(i_s+1) - ia_s(i_s)
       p_i = ia_s(i_s)
       p_f = ia_s(i_s+1)-1
       ja_work(i_glob,nja_glob(i_glob)+1:nja_glob(i_glob)+nnz) = ja_s(p_i:p_f) + np_b
       nja_glob(i_glob) =  nja_glob(i_glob) + nnz
    END DO

    NULLIFY(ia_s,ja_s,nja_b)

    nnz = SUM(nja_glob)
    ALLOCATE(ia_glob(np_glob+1),ja_glob(nnz))
    ia_glob(1) = 1
    p_f = ia_glob(1) - 1
    DO i_glob = 1, np_glob 
       CALL tri_jlg (ja_work(i_glob,1:nja_glob(i_glob)), a_d, n_a_d)
       IF (n_a_d /= nja_glob(i_glob)) THEN
          WRITE(*,*) ' BUG in st_scr_maxwell_tm '
          WRITE(*,*) 'n_a_d ', n_a_d, 'nja_glob(i_glob)', nja_glob(i_glob)
          STOP
       END IF
       p_i = p_f + 1
       p_f =  p_f + nja_glob(i_glob)
       ia_glob(i_glob+1) = ia_glob(i_glob) + nja_glob(i_glob)
       ja_glob(p_i:p_f) = a_d(1:nja_glob(i_glob))
    END DO

    IF (p_f /= nnz) THEN
       WRITE(*,*) ' BUG in st_scr_maxwell_tm '
       STOP
    END IF

    NULLIFY(ja_work,nja_glob)  

  END SUBROUTINE st_scr_maxwell_tm
  !-----------------------------------------------------------------------

  SUBROUTINE st_scr_maxwell_tmag(H_mesh, phi_mesh, INTERFACE, ia_glob, ja_glob)

    USE def_type_mesh

    IMPLICIT NONE

    TYPE(mesh_type),           INTENT(IN) :: H_mesh, phi_mesh
    TYPE(interface_type),      INTENT(IN) :: INTERFACE
    INTEGER, DIMENSION(:),     POINTER    :: ia_glob, ja_glob

    INTEGER, DIMENSION(:),     POINTER    :: ia_m, ja_m, ia_s, ja_s, ia_b, &
         ja_b, nja_glob, nja_b, a_d
    INTEGER, DIMENSION(:,:),   POINTER    :: ja_work
    INTEGER, PARAMETER :: param=90
    INTEGER :: np_m, np_s, np_b, np_glob, i_glob, j_glob, i, i_b, i_s, &
         nnz, p_i, p_f, ms, n_a_d, m1, m2, n1, n2

    CALL st_csr(H_mesh%jj, ia_m, ja_m)
    np_m = SIZE(ia_m)-1
    CALL st_csr(phi_mesh%jj, ia_s, ja_s)
    CALL st_csr_bloc(ia_m, ja_m, ia_b, ja_b, 2)
    NULLIFY(ia_m,ja_m)

    np_b = SIZE(ia_b)-1
    np_s = SIZE(ia_s)-1
    np_glob = np_b + np_s

    ALLOCATE(ja_work(np_glob,param),a_d(param),nja_b(np_b),nja_glob(np_glob))

    nja_glob = 0
    DO i_b = 1, np_b
       nnz = ia_b(i_b+1) - ia_b(i_b)
       p_i = ia_b(i_b)
       p_f = ia_b(i_b+1) - 1
       ja_work(i_b,1:nnz) = ja_b(p_i:p_f)
       nja_b(i_b) = nnz
    END DO
    nja_glob(1:np_b) = nja_b

    !   DO ms = 1, interface%mes
    !      DO ns1 = 1, H_mesh%gauss%n_ws
    !         i = interface%jjs1(ns1,ms)
    !         DO ns2 = 1, phi_mesh%gauss%n_ws
    !            j_glob = interface%jjs2(ns2,ms) + np_b
    !            IF (MINVAL(ABS(ja_work(i,nja_b(i):nja_glob(i))-j_glob)) == 0) CYCLE 
    !            nja_glob(i) = nja_glob(i) + 1
    !            nja_glob(i+np_m) = nja_glob(i+np_m) + 1
    !            ja_work(i,nja_glob(i)) = j_glob
    !            ja_work(i+np_m,nja_glob(i+np_m)) = j_glob
    !            nja_glob(j_glob) = nja_glob(j_glob) + 2
    !            ja_work(j_glob,nja_glob(j_glob)-1) = i
    !            ja_work(j_glob,nja_glob(j_glob)) = i + np_m
    !         END DO
    !      END DO
    !   END DO

    DO ms = 1, interface%mes
       m1 = H_mesh%neighs(interface%mesh1(ms))
       m2 = phi_mesh%neighs(interface%mesh2(ms))
       DO n1 = 1, H_mesh%gauss%n_w
          i = H_mesh%jj(n1,m1)
          DO n2 = 1, phi_mesh%gauss%n_w
             j_glob =  phi_mesh%jj(n2,m2) + np_b
             IF (MINVAL(ABS(ja_work(i,nja_b(i):nja_glob(i))-j_glob)) == 0) CYCLE 
             nja_glob(i) = nja_glob(i) + 1
             nja_glob(i+np_m) = nja_glob(i+np_m) + 1
             ja_work(i,nja_glob(i)) = j_glob
             ja_work(i+np_m,nja_glob(i+np_m)) = j_glob
             nja_glob(j_glob) = nja_glob(j_glob) + 2
             ja_work(j_glob,nja_glob(j_glob)-1) = i
             ja_work(j_glob,nja_glob(j_glob)) = i + np_m
          END DO
       END DO
    END DO




    DO i_s = 1, np_s
       i_glob =  i_s + np_b
       nnz = ia_s(i_s+1) - ia_s(i_s)
       p_i = ia_s(i_s)
       p_f = ia_s(i_s+1)-1
       ja_work(i_glob,nja_glob(i_glob)+1:nja_glob(i_glob)+nnz) = ja_s(p_i:p_f) + np_b
       nja_glob(i_glob) =  nja_glob(i_glob) + nnz
    END DO

    NULLIFY(ia_s,ja_s,nja_b)

    nnz = SUM(nja_glob)
    ALLOCATE(ia_glob(np_glob+1),ja_glob(nnz))
    ia_glob(1) = 1
    p_f = ia_glob(1) - 1
    DO i_glob = 1, np_glob 
       CALL tri_jlg (ja_work(i_glob,1:nja_glob(i_glob)), a_d, n_a_d)
       IF (n_a_d /= nja_glob(i_glob)) THEN
          WRITE(*,*) ' BUG in st_scr_maxwell_tmag '
          WRITE(*,*) 'n_a_d ', n_a_d, 'nja_glob(i_glob)', nja_glob(i_glob)
          STOP
       END IF
       p_i = p_f + 1
       p_f =  p_f + nja_glob(i_glob)
       ia_glob(i_glob+1) = ia_glob(i_glob) + nja_glob(i_glob)
       ja_glob(p_i:p_f) = a_d(1:nja_glob(i_glob))
    END DO

    IF (p_f /= nnz) THEN
       WRITE(*,*) ' BUG in st_scr_maxwell_tm '
       STOP
    END IF

    NULLIFY(ja_work,nja_glob)  

  END SUBROUTINE st_scr_maxwell_tmag

  SUBROUTINE st_scr_maxwell_tmag3d(H_mesh, phi_mesh, INTERFACE, ia_glob, ja_glob)

    USE def_type_mesh

    IMPLICIT NONE

    TYPE(mesh_type),           INTENT(IN) :: H_mesh, phi_mesh
    TYPE(interface_type),      INTENT(IN) :: INTERFACE
    INTEGER, DIMENSION(:),     POINTER    :: ia_glob, ja_glob

    INTEGER, DIMENSION(:),     POINTER    :: ia_m, ja_m, ia_s, ja_s, ia_b, &
         ja_b, nja_glob, nja_b, a_d, &
         ia_s1, ja_s1
    INTEGER, DIMENSION(:,:),   POINTER    :: ja_work
    INTEGER, PARAMETER :: param=90
    INTEGER :: np_m, np_s, np_b, np_glob, i_glob, j_glob, i, i_b, i_s, &
         nnz, p_i, p_f, ms, n_a_d, m1, m2, n1, n2

    CALL st_csr(H_mesh%jj, ia_m, ja_m)
    np_m = SIZE(ia_m)-1
    CALL st_csr(phi_mesh%jj, ia_s1, ja_s1)
    CALL st_csr_bloc(ia_s1, ja_s1, ia_s, ja_s, 2)
    CALL st_csr_bloc(ia_m, ja_m, ia_b, ja_b, 6)
    NULLIFY(ia_m,ja_m)

    np_b = SIZE(ia_b)-1
    np_s = SIZE(ia_s)-1
    np_glob = np_b + np_s

    ALLOCATE(ja_work(np_glob,param),a_d(param),nja_b(np_b),nja_glob(np_glob))

    nja_glob = 0
    DO i_b = 1, np_b
       nnz = ia_b(i_b+1) - ia_b(i_b)
       p_i = ia_b(i_b)
       p_f = ia_b(i_b+1) - 1
       ja_work(i_b,1:nnz) = ja_b(p_i:p_f)
       nja_b(i_b) = nnz
    END DO
    nja_glob(1:np_b) = nja_b

    !   DO ms = 1, interface%mes
    !      DO ns1 = 1, H_mesh%gauss%n_ws
    !         i = interface%jjs1(ns1,ms)
    !         DO ns2 = 1, phi_mesh%gauss%n_ws
    !            j_glob = interface%jjs2(ns2,ms) + np_b
    !            IF (MINVAL(ABS(ja_work(i,nja_b(i):nja_glob(i))-j_glob)) == 0) CYCLE 
    !            nja_glob(i) = nja_glob(i) + 1
    !            nja_glob(i+np_m) = nja_glob(i+np_m) + 1
    !            ja_work(i,nja_glob(i)) = j_glob
    !            ja_work(i+np_m,nja_glob(i+np_m)) = j_glob
    !            nja_glob(j_glob) = nja_glob(j_glob) + 2
    !            ja_work(j_glob,nja_glob(j_glob)-1) = i
    !            ja_work(j_glob,nja_glob(j_glob)) = i + np_m
    !         END DO
    !      END DO
    !   END DO

    DO ms = 1, interface%mes
       m1 = H_mesh%neighs(interface%mesh1(ms))
       m2 = phi_mesh%neighs(interface%mesh2(ms))
       DO n1 = 1, H_mesh%gauss%n_w
          i = H_mesh%jj(n1,m1)
          DO n2 = 1, phi_mesh%gauss%n_w
             j_glob =  phi_mesh%jj(n2,m2) + np_b
             IF (MINVAL(ABS(ja_work(i,nja_b(i):nja_glob(i))-j_glob)) == 0) CYCLE 
             nja_glob(i) = nja_glob(i) + 1
             nja_glob(i+np_m) = nja_glob(i+np_m) + 1
             ja_work(i,nja_glob(i)) = j_glob
             ja_work(i+np_m,nja_glob(i+np_m)) = j_glob
             nja_glob(j_glob) = nja_glob(j_glob) + 2
             ja_work(j_glob,nja_glob(j_glob)-1) = i
             ja_work(j_glob,nja_glob(j_glob)) = i + np_m
          END DO
       END DO
    END DO




    DO i_s = 1, np_s
       i_glob =  i_s + np_b
       nnz = ia_s(i_s+1) - ia_s(i_s)
       p_i = ia_s(i_s)
       p_f = ia_s(i_s+1)-1
       ja_work(i_glob,nja_glob(i_glob)+1:nja_glob(i_glob)+nnz) = ja_s(p_i:p_f) + np_b
       nja_glob(i_glob) =  nja_glob(i_glob) + nnz
    END DO

    NULLIFY(ia_s,ja_s,nja_b)

    nnz = SUM(nja_glob)
    ALLOCATE(ia_glob(np_glob+1),ja_glob(nnz))
    ia_glob(1) = 1
    p_f = ia_glob(1) - 1
    DO i_glob = 1, np_glob 
       CALL tri_jlg (ja_work(i_glob,1:nja_glob(i_glob)), a_d, n_a_d)
       IF (n_a_d /= nja_glob(i_glob)) THEN
          WRITE(*,*) ' BUG in st_scr_maxwell_tmag '
          WRITE(*,*) 'n_a_d ', n_a_d, 'nja_glob(i_glob)', nja_glob(i_glob)
          STOP
       END IF
       p_i = p_f + 1
       p_f =  p_f + nja_glob(i_glob)
       ia_glob(i_glob+1) = ia_glob(i_glob) + nja_glob(i_glob)
       ja_glob(p_i:p_f) = a_d(1:nja_glob(i_glob))
    END DO

    IF (p_f /= nnz) THEN
       WRITE(*,*) ' BUG in st_scr_maxwell_tm '
       STOP
    END IF

    NULLIFY(ja_work,nja_glob)  

  END SUBROUTINE st_scr_maxwell_tmag3d

  SUBROUTINE st_scr_maxwell_tmag3d_decouple(H_mesh, phi_mesh, INTERFACE, ia_glob, ja_glob)

    USE def_type_mesh

    IMPLICIT NONE

    TYPE(mesh_type),           INTENT(IN) :: H_mesh, phi_mesh
    TYPE(interface_type),      INTENT(IN) :: INTERFACE
    INTEGER, DIMENSION(:),     POINTER    :: ia_glob, ja_glob

    INTEGER, DIMENSION(:),     POINTER    :: ia_m, ja_m, ia_s, ja_s, ia_b, &
         ja_b, nja_glob, nja_b, a_d, &
         ia_s1, ja_s1
    INTEGER, DIMENSION(:,:),   POINTER    :: ja_work
    INTEGER, PARAMETER :: param=90
    INTEGER :: np_m, np_s, np_b, np_glob, i_glob, j_glob, i, i_b, i_s, &
         nnz, p_i, p_f, ms, n_a_d, m1, m2, n1, n2


    CALL st_csr(H_mesh%jj, ia_m, ja_m)
    np_m = SIZE(ia_m)-1
    CALL st_csr(phi_mesh%jj, ia_s, ja_s)
    CALL st_csr_bloc(ia_m, ja_m, ia_b, ja_b, 3)
    NULLIFY(ia_m,ja_m)           

    np_b = SIZE(ia_b)-1
    np_s = SIZE(ia_s)-1
    np_glob = np_b + np_s

    ALLOCATE(ja_work(np_glob,param),a_d(param),nja_b(np_b),nja_glob(np_glob))

    nja_glob = 0
    DO i_b = 1, np_b
       nnz = ia_b(i_b+1) - ia_b(i_b)
       p_i = ia_b(i_b)
       p_f = ia_b(i_b+1) - 1
       ja_work(i_b,1:nnz) = ja_b(p_i:p_f)
       nja_b(i_b) = nnz
    END DO
    nja_glob(1:np_b) = nja_b

    !   DO ms = 1, interface%mes
    !      DO ns1 = 1, H_mesh%gauss%n_ws
    !         i = interface%jjs1(ns1,ms)
    !         DO ns2 = 1, phi_mesh%gauss%n_ws
    !            j_glob = interface%jjs2(ns2,ms) + np_b
    !            IF (MINVAL(ABS(ja_work(i,nja_b(i):nja_glob(i))-j_glob)) == 0) CYCLE 
    !            nja_glob(i) = nja_glob(i) + 1
    !            nja_glob(i+np_m) = nja_glob(i+np_m) + 1
    !            ja_work(i,nja_glob(i)) = j_glob
    !            ja_work(i+np_m,nja_glob(i+np_m)) = j_glob
    !            nja_glob(j_glob) = nja_glob(j_glob) + 2
    !            ja_work(j_glob,nja_glob(j_glob)-1) = i
    !            ja_work(j_glob,nja_glob(j_glob)) = i + np_m
    !         END DO
    !      END DO
    !   END DO

    DO ms = 1, interface%mes
       m1 = H_mesh%neighs(interface%mesh1(ms))
       m2 = phi_mesh%neighs(interface%mesh2(ms))
       DO n1 = 1, H_mesh%gauss%n_w
          i = H_mesh%jj(n1,m1)
          DO n2 = 1, phi_mesh%gauss%n_w
             j_glob =  phi_mesh%jj(n2,m2) + np_b
             IF (MINVAL(ABS(ja_work(i,nja_b(i):nja_glob(i))-j_glob)) == 0) CYCLE 
             nja_glob(i) = nja_glob(i) + 1
             nja_glob(i+np_m) = nja_glob(i+np_m) + 1
             nja_glob(i+2*np_m) = nja_glob(i+2*np_m) + 1
             ja_work(i,nja_glob(i)) = j_glob
             ja_work(i+np_m,nja_glob(i+np_m)) = j_glob
             ja_work(i+2*np_m,nja_glob(i+2*np_m)) = j_glob
             nja_glob(j_glob) = nja_glob(j_glob) + 3
             ja_work(j_glob,nja_glob(j_glob)-2) = i
             ja_work(j_glob,nja_glob(j_glob)-1) = i + np_m
             ja_work(j_glob,nja_glob(j_glob)) = i + 2*np_m
          END DO
       END DO
    END DO


    DO i_s = 1, np_s
       i_glob =  i_s + np_b
       nnz = ia_s(i_s+1) - ia_s(i_s)
       p_i = ia_s(i_s)
       p_f = ia_s(i_s+1)-1
       ja_work(i_glob,nja_glob(i_glob)+1:nja_glob(i_glob)+nnz) = ja_s(p_i:p_f) + np_b
       nja_glob(i_glob) =  nja_glob(i_glob) + nnz
    END DO

    NULLIFY(ia_s,ja_s,nja_b)

    nnz = SUM(nja_glob)
    ALLOCATE(ia_glob(np_glob+1),ja_glob(nnz))
    ia_glob(1) = 1
    p_f = ia_glob(1) - 1
    DO i_glob = 1, np_glob 
       CALL tri_jlg (ja_work(i_glob,1:nja_glob(i_glob)), a_d, n_a_d)
       IF (n_a_d /= nja_glob(i_glob)) THEN
          WRITE(*,*) ' BUG in st_scr_maxwell_tmag '
          WRITE(*,*) 'n_a_d ', n_a_d, 'nja_glob(i_glob)', nja_glob(i_glob)
          STOP
       END IF
       p_i = p_f + 1
       p_f =  p_f + nja_glob(i_glob)
       ia_glob(i_glob+1) = ia_glob(i_glob) + nja_glob(i_glob)
       ja_glob(p_i:p_f) = a_d(1:nja_glob(i_glob))
    END DO

    IF (p_f /= nnz) THEN
       WRITE(*,*) ' BUG in st_scr_maxwell_tm '
       STOP
    END IF

    NULLIFY(ja_work,nja_glob)  

  END SUBROUTINE st_scr_maxwell_tmag3d_decouple

  SUBROUTINE st_scr_maxwell_tmag3d_decouple_mu(H_mesh, phi_mesh, INTERFACE, interface_H_mu, ia_glob, ja_glob)

    USE def_type_mesh

    IMPLICIT NONE

    TYPE(mesh_type),           INTENT(IN) :: H_mesh, phi_mesh
    TYPE(interface_type),      INTENT(IN) :: INTERFACE, interface_H_mu
    INTEGER, DIMENSION(:),     POINTER    :: ia_glob, ja_glob

    INTEGER, DIMENSION(:),     POINTER    :: ia_m, ja_m, ia_s, ja_s, ia_b, &
         ja_b, nja_glob, nja_b, a_d, &
         ia_s1, ja_s1
    INTEGER, DIMENSION(:,:),   POINTER    :: ja_work
    INTEGER, PARAMETER :: param=120
    INTEGER :: np_m, np_s, np_b, np_glob, i_glob, j_glob, i, i_b, i_s, &
         nnz, p_i, p_f, ms, n_a_d, m1, m2, ni, nj, ki, kj, ib, jb, ci, cj, mi, mj, n1, n2


    CALL st_csr(H_mesh%jj, ia_m, ja_m)
    np_m = SIZE(ia_m)-1
    CALL st_csr(phi_mesh%jj, ia_s, ja_s)
    CALL st_csr_bloc(ia_m, ja_m, ia_b, ja_b, 3)
    NULLIFY(ia_m,ja_m)

    np_b = SIZE(ia_b)-1
    np_s = SIZE(ia_s)-1
    np_glob = np_b + np_s

    ALLOCATE(ja_work(np_glob,param),a_d(param),nja_b(np_b),nja_glob(np_glob))

    nja_glob = 0
    DO i_b = 1, np_b
       nnz = ia_b(i_b+1) - ia_b(i_b)
       p_i = ia_b(i_b)
       p_f = ia_b(i_b+1) - 1
       ja_work(i_b,1:nnz) = ja_b(p_i:p_f)
       nja_b(i_b) = nnz
    END DO
    nja_glob(1:np_b) = nja_b


    !  Interface_H_mu
    DO ms = 1, interface_H_mu%mes
       m1 = H_mesh%neighs(interface_H_mu%mesh1(ms))
       m2 = H_mesh%neighs(interface_H_mu%mesh2(ms))
       DO ci = 1, 2
          DO ki = 1, 3
             DO ni = 1, H_mesh%gauss%n_w
                IF (ci==1) THEN
                   mi = m1
                ELSE
                   mi = m2
                END IF
                ib = H_mesh%jj(ni,mi) + (ki-1)*np_m
                DO cj = 1, 2
                   DO kj = 1, 3
                      DO nj = 1, H_mesh%gauss%n_w
                         IF (cj==1) THEN
                            mj = m1
                         ELSE
                            mj = m2
                         END IF
                         jb =  H_mesh%jj(nj,mj) + (kj-1)*np_m 
                         IF (MINVAL(ABS(ja_work(ib,1:nja_glob(ib))-jb)) == 0) CYCLE 
                         nja_glob(ib) = nja_glob(ib) + 1
                         ja_work(ib,nja_glob(ib)) = jb
                      END DO
                   END DO
                END DO
             END DO
          END DO
       END DO
    END DO

    ! Interface_H_phi
    DO ms = 1, interface%mes
       m1 = H_mesh%neighs(interface%mesh1(ms))
       m2 = phi_mesh%neighs(interface%mesh2(ms))
       DO n1 = 1, H_mesh%gauss%n_w
          i = H_mesh%jj(n1,m1)
          DO n2 = 1, phi_mesh%gauss%n_w
             j_glob =  phi_mesh%jj(n2,m2) + np_b
             IF (MINVAL(ABS(ja_work(i,1:nja_glob(i))-j_glob)) == 0) CYCLE
             nja_glob(i) = nja_glob(i) + 1
             nja_glob(i+np_m) = nja_glob(i+np_m) + 1
             nja_glob(i+2*np_m) = nja_glob(i+2*np_m) + 1
             ja_work(i,nja_glob(i)) = j_glob
             ja_work(i+np_m,nja_glob(i+np_m)) = j_glob
             ja_work(i+2*np_m,nja_glob(i+2*np_m)) = j_glob
             nja_glob(j_glob) = nja_glob(j_glob) + 3
             ja_work(j_glob,nja_glob(j_glob)-2) = i
             ja_work(j_glob,nja_glob(j_glob)-1) = i + np_m
             ja_work(j_glob,nja_glob(j_glob)) = i + 2*np_m
          END DO
       END DO
    END DO

    DO i_s = 1, np_s
       i_glob =  i_s + np_b
       nnz = ia_s(i_s+1) - ia_s(i_s)
       p_i = ia_s(i_s)
       p_f = ia_s(i_s+1)-1
       ja_work(i_glob,nja_glob(i_glob)+1:nja_glob(i_glob)+nnz) = ja_s(p_i:p_f) + np_b
       nja_glob(i_glob) =  nja_glob(i_glob) + nnz
    END DO

    NULLIFY(ia_s,ja_s,nja_b)

    nnz = SUM(nja_glob)
    ALLOCATE(ia_glob(np_glob+1),ja_glob(nnz))
    ia_glob(1) = 1
    p_f = ia_glob(1) - 1
    DO i_glob = 1, np_glob 
       CALL tri_jlg (ja_work(i_glob,1:nja_glob(i_glob)), a_d, n_a_d)
       IF (n_a_d /= nja_glob(i_glob)) THEN
          WRITE(*,*) ' BUG in st_scr_maxwell_tmag_mu, i_glob = ', i_glob
          WRITE(*,*) 'n_a_d ', n_a_d, 'nja_glob(i_glob)', nja_glob(i_glob)
          WRITE(*,*) ja_work(i_glob,:nja_glob(i_glob))
          STOP
       END IF
       p_i = p_f + 1
       p_f =  p_f + nja_glob(i_glob)
       ia_glob(i_glob+1) = ia_glob(i_glob) + nja_glob(i_glob)
       ja_glob(p_i:p_f) = a_d(1:nja_glob(i_glob))
    END DO

    IF (p_f /= nnz) THEN
       WRITE(*,*) ' BUG in st_scr_maxwell_tmag_mu '
       STOP
    END IF

    NULLIFY(ja_work,nja_glob)  

  END SUBROUTINE st_scr_maxwell_tmag3d_decouple_mu

END MODULE st_matrix










