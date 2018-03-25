PROGRAM fft_parallele
  
  USE sft_parallele
  IMPLICIT NONE
  include 'mpif.h'

  INTEGER, PARAMETER :: m_max=24, np = 10000
  INTEGER            :: code, rang, nb_procs, m_max_c

  REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: V1r_p, V1t_p, V1z_p, V2r_p, V2t_p, V2z_p
  REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: V1r,   V1t,   V1z,   V2r,   V2t,   V2z
  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: V1_in, V2_in, V_out
  REAL(KIND=8), DIMENSION(3) :: temps  
  REAL(KIND=8) :: t
  INTEGER :: n, i

  CALL MPI_INIT(code)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,code)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)

  temps = 0
  m_max_c=m_max/nb_procs
  IF (MOD(m_max,nb_procs)/=0 .OR. m_max_c==0) THEN
     WRITE(*,*) ' BUG '
     STOP
  END IF
 
  ALLOCATE(V1r_p(2*m_max_c,np),V1t_p(2*m_max_c,np),V1z_p(2*m_max_c,np))
  ALLOCATE(V2r_p(2*m_max_c,np),V2t_p(2*m_max_c,np),V2z_p(2*m_max_c,np))
  ALLOCATE(V1r(2*m_max_c,np),V1t(2*m_max_c,np),V1z(2*m_max_c,np))
  ALLOCATE(V2r(2*m_max_c,np),V2t(2*m_max_c,np),V2z(2*m_max_c,np))

  DO n = 1, np
     DO i = 1, 2*m_max_c
        V1r_p(i,n) = SIN(COS(1234.d0*i*n))
        V1t_p(i,n) = COS(SIN(1234.d0*i*n))
        V1z_p(i,n) = SIN(SIN(1234.d0*i*n))
        V2r_p(i,n) = COS(COS(1234.d0*i*n))
        V2t_p(i,n) = COS(SIN(4321.d0*i*n))
        V2z_p(i,n) = SIN(COS(4321.d0*i*n))
     END DO
  END DO

  CALL FFT_PAR_MP_b(V1r_p, V1r, 2)
  CALL FFT_PAR_MP_b(V1t_p, V1t, 2)
  CALL FFT_PAR_MP_b(V1z_p, V1z, 2)
  CALL FFT_PAR_MP_b(V2r_p, V2r, 2)
  CALL FFT_PAR_MP_b(V2t_p, V2t, 2)
  CALL FFT_PAR_MP_b(V2z_p, V2z, 2)

  ALLOCATE(V1_in(6,np,m_max_c),V2_in(6,np,m_max_c),V_out(6,np,m_max_c))

  DO i =1 , m_max_c
     V1_in(1,:,i) = V1r(2*i-1,:)
     V1_in(2,:,i) = V1r(2*i,:)
     V1_in(3,:,i) = V1t(2*i-1,:)
     V1_in(4,:,i) = V1t(2*i,:)
     V1_in(5,:,i) = V1z(2*i-1,:)
     V1_in(6,:,i) = V1z(2*i,:)
     V2_in(1,:,i) = V2r(2*i-1,:)
     V2_in(2,:,i) = V2r(2*i,:)
     V2_in(3,:,i) = V2t(2*i-1,:)
     V2_in(4,:,i) = V2t(2*i,:)
     V2_in(5,:,i) = V2z(2*i-1,:)
     V2_in(6,:,i) = V2z(2*i,:)
  END DO

  V1r= V1t_p*V2z_p - V1z_p*V2t_p
  V1t= V1z_p*V2r_p - V1r_p*V2z_p
  V1z= V1r_p*V2t_p - V1t_p*V2r_p

  CALL FFT_PAR_MP_b(V1r, V2r, 2)
  CALL FFT_PAR_MP_b(V1t, V2t, 2)
  CALL FFT_PAR_MP_b(V1z, V2z, 2)

  t = MPI_WTIME()
  CALL FFT_PAR_PROD_VECT(V1_in, V2_in, V_out, temps)
  CALL FFT_PAR_PROD_VECT(V1_in, V2_in, V_out, temps)
  t = MPI_WTIME() - t

  IF(rang==0) WRITE(*,*) ' Temps de communication', temps(1), temps(1)/(2*m_max*np)
  IF(rang==0) WRITE(*,*) ' Temps de calcul       ', temps(2), temps(2)/(2*m_max*np)
  IF(rang==0) WRITE(*,*) ' Temps de transfert    ', temps(3), temps(3)/(2*m_max*np)
  IF(rang==0) WRITE(*,*) ' Temps total ', t, SUM(temps)

  DO i =1 , m_max_c
     V1_in(1,:,i) = V2r(2*i-1,:)
     V1_in(2,:,i) = V2r(2*i,:)
     V1_in(3,:,i) = V2t(2*i-1,:)
     V1_in(4,:,i) = V2t(2*i,:)
     V1_in(5,:,i) = V2z(2*i-1,:)
     V1_in(6,:,i) = V2z(2*i,:)
  END DO 

  WRITE(*,*) ' rang, Erreur', rang, MAXVAL(ABS(V1_in-V_out))
  DEALLOCATE(V1_in, V2_in)
  DEALLOCATE(V1r, V1t, V1z, V2r, V2t, V2z)

  CALL MPI_FINALIZE(code)

END PROGRAM fft_parallele
   
