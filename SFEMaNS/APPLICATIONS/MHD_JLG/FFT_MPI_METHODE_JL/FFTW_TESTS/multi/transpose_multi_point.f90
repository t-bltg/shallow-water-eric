program alltoall
 implicit none

 INCLUDE 'mpif.h'
 INTEGER, PARAMETER :: m_max=4, np=7 
 INTEGER :: np_tot, nb_valeurs, m_max_c, bloc_size
 INTEGER :: nb_procs, rang, longueur_tranche, i, n, code, MPID
 REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: u1_dist, u1_glued 

 call MPI_INIT (code)
 call MPI_COMM_SIZE ( MPI_COMM_WORLD ,nb_procs,code)
 call MPI_COMM_RANK ( MPI_COMM_WORLD ,rang,code)

 m_max_c=m_max/nb_procs
 IF (MOD(m_max,nb_procs)/=0 .OR. m_max_c==0) THEN
    WRITE(*,*) ' BUG '
    STOP
 END IF

 IF (MODULO(np,nb_procs)==0) THEN
    bloc_size = np/nb_procs
 ELSE
    bloc_size = np/nb_procs + 1
 END IF
 np_tot = nb_procs*bloc_size

 ALLOCATE(u1_dist(m_max_c,np_tot))
 ALLOCATE(u1_glued(m_max_c,np_tot))

 DO n =  1, np
    DO i = 1, m_max_c
       u1_dist(i,n) = (n-1)*m_max_c + i-1 
    END DO
 END DO
 DO n =  np+1, np_tot
    u1_dist(:,n) = 0
 END DO

 longueur_tranche=bloc_size*m_max_c

 WRITE(*,*) ' AVANT'
 WRITE(*,'(I2,12(F3.0,1x)') rang, u1_dist 
 MPID=MPI_DOUBLE_PRECISION
 call MPI_ALLTOALL (u1_dist,longueur_tranche, MPID , u1_glued,longueur_tranche, &
                    MPID , MPI_COMM_WORLD ,code)

 WRITE(*,*) ' APRES'
 WRITE(*,'(I2,12(F3.0,1x)') rang, u1_glued 
 !print *,'Moi, processus ',rang, ', j''ai recu' ,donnees(1),' ... ', &
 !donnees(longueur_tranche+1),' ... ',donnees(nb_valeurs)

 call MPI_FINALIZE (code)

 end program alltoall
