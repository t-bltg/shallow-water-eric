PROGRAM test

  USE sft_parallele 
  IMPLICIT NONE

  INCLUDE 'mpif.h'
  INCLUDE 'fftw3.f'

  INTEGER, PARAMETER :: m_max=32, np=100000, nb_field=6
  INTEGER :: m_max_c, N_r
  INTEGER :: nb_procs, rang, i, n, nf, p, code, MPID, it=1
  REAL(KIND=8) :: error, t
  REAL(KIND=8), DIMENSION(3) :: temps, cumul_t 
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: u1, u2, u3 
  COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: cu, prod_cu 
  COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: uuc, uc, vc
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:)    :: uur, ur, vvvr 

  ! FFTW parameters
  ! It is extremly important that the fftw plan are declared with  INTEGER(KIND=8)
  INTEGER(KIND=8) :: fftw_plan_test_c2r, fftw_plan_test_r2c
  INTEGER :: fft_dim, howmany, istride, ostride, idist, odist
  INTEGER, DIMENSION(1) :: dim, inembed, onembed
  ! Recall complexes must be rescaled
  ! End FFTW parameters

  CALL MPI_INIT (code)
  CALL MPI_COMM_SIZE (MPI_COMM_WORLD,nb_procs,code)
  CALL MPI_COMM_RANK (MPI_COMM_WORLD,rang,code)

  N_r=2*m_max-1
  m_max_c=m_max/nb_procs
  IF (MOD(m_max,nb_procs)/=0 .OR. m_max_c==0) THEN
     WRITE(*,*) ' BUG '
     CALL MPI_FINALIZE (code)
     STOP
  END IF

  ! Definition of u1 and u2
  ALLOCATE(u1(np,nb_field,m_max_c))
  ALLOCATE(u2(np,nb_field,m_max_c))
  ALLOCATE(u3(np,nb_field,m_max_c))

  DO n =  1, np
     DO i = 1, m_max_c
        DO nf = 1, nb_field/2 
           u1(n,2*nf-1,i) = nf*(n+i + rang*m_max_c)*1.d0/i**2
           u1(n,2*nf,i)   = (nf-i)*(n+i + rang*m_max_c)*1.d0/i**2
           u2(n,2*nf-1,i) = COS(nf*(n+i + rang*m_max_c)*1.d0/i**2)
           u2(n,2*nf,i)   = SIN((nf-i)*(n+i + rang*m_max_c)*1.d0/i**2)**2
        END DO
     END DO
  END DO

  ! Definition of u1 and u2 with all the modes on the processor
  ALLOCATE(uuc(m_max,np,nb_field/2),uc(m_max,np,nb_field/2), vc(m_max,np,nb_field/2))
  ALLOCATE(uur(N_r,np,nb_field/2))
  ALLOCATE(ur(N_r,np,nb_field/2))
  ALLOCATE(vvvr(N_r,np,nb_field/2))

  DO n =  1, np
     DO p = 1, nb_procs
        DO i = 1, m_max_c
           DO nf = 1, nb_field/2 
              IF (i + (p-1)*m_max_c==1) THEN 
                 uc(i+(p-1)*m_max_c,n,nf) = CMPLX(nf*(n+i + (p-1)*m_max_c)*1.d0/i**2,0.d0,KIND=8)
                 vc(i+(p-1)*m_max_c,n,nf) = CMPLX(COS(nf*(n+i + (p-1)*m_max_c)*1.d0/i**2),0.d0,KIND=8)
              ELSE
                 uc(i+(p-1)*m_max_c,n,nf) = 0.5d0*CMPLX(nf*(n+i + (p-1)*m_max_c)*1.d0/i**2, &
                      -(nf-i)*(n+i + (p-1)*m_max_c)*1.d0/i**2,KIND=8)
                 vc(i+(p-1)*m_max_c,n,nf) = 0.5d0*CMPLX(COS(nf*(n+i + (p-1)*m_max_c)*1.d0/i**2),&
                      -SIN((nf-i)*(n+i + (p-1)*m_max_c)*1.d0/i**2)**2,KIND=8)
              END IF
          END DO
       END DO
    END DO
 END DO
 fft_dim=1; istride=1; ostride=1; 
 idist=N_r;   inembed(1)=N_r; DIM(1)=N_r
 odist=m_max; onembed(1)=m_max
 howmany = np*nb_field/2
 uuc = uc
 CALL dfftw_plan_many_dft_c2r(fftw_plan_test_c2r, fft_dim, dim, howmany, uuc, &
                onembed, istride, odist, uur, inembed, istride, idist, FFTW_ESTIMATE)
 CALL dfftw_execute(fftw_plan_test_c2r)
 ur(:,:,:) = uur(:,:,:)

 ! There is a bug with the optimization -O here
 !ur = uur
 uuc = vc
 CALL dfftw_execute(fftw_plan_test_c2r)
 vvvr = uur

 !CROSS PRODUCT
 uur(:,:,1) = ur(:,:,2)*vvvr(:,:,3) - ur(:,:,3)*vvvr(:,:,2) 
 uur(:,:,2) = ur(:,:,3)*vvvr(:,:,1) - ur(:,:,1)*vvvr(:,:,3) 
 uur(:,:,3) = ur(:,:,1)*vvvr(:,:,2) - ur(:,:,2)*vvvr(:,:,1) 
 ur = uur
 !END CROSS PRODUCT

 CALL dfftw_plan_many_dft_r2c(fftw_plan_test_r2c, fft_dim, dim, howmany, ur, &
                inembed, istride, idist, uc, onembed, ostride, odist, FFTW_ESTIMATE);
 CALL dfftw_execute(fftw_plan_test_r2c)
 uc = uc/N_r !========Scaling
 !BACK to COSINE AND SINE FORMAT
 uc(2:,:,:)=2*CONJG(uc(2:,:,:))
 !WRITE(in_unit,*) ' Direct results'
 ! End definition of u1 and u2 with all the modes on the processor

 
 cumul_t = 0.d0
 t = MPI_WTIME()
 DO i = 1, it
    CALL ref(u1, u2, u3, temps)
    cumul_t = cumul_t + temps
 END DO
 IF(rang==0) WRITE(*,*) ' Temps de communication', cumul_t(1), cumul_t(1)/(6*m_max*np*it)
 IF(rang==0) WRITE(*,*) ' Temps de calcul       ', cumul_t(2), cumul_t(2)/(6*m_max*np*it)
 IF(rang==0) WRITE(*,*) ' Temps de transfert    ', cumul_t(3), cumul_t(3)/(6*m_max*np*it)
 t = MPI_WTIME() - t
 IF(rang==0) WRITE(*,*) ' Temps total ', t, SUM(temps)

 DO n = 1, np
    DO i = 1, m_max_c
       DO nf = 1, nb_field/2
          error = MAX(error,ABS(REAL(uc(rang*m_max_c+i,n,nf),KIND=8)-u3(n,2*nf-1,i)) &
                          + ABS(AIMAG(uc(rang*m_max_c+i,n,nf))-u3(n,2*nf,i)))
       END DO
    END DO
 END DO
 WRITE(*,*) ' ERROR ', error 

 cumul_t = 0.d0
 t = MPI_WTIME()
 DO i = 1, it
    CALL fft_par_cross_prod(u1, u2, u3, temps)
    cumul_t = cumul_t + temps
 END DO
 IF(rang==0) WRITE(*,*) ' Temps de communication', cumul_t(1), cumul_t(1)/(2*m_max*np)
 IF(rang==0) WRITE(*,*) ' Temps de calcul       ', cumul_t(2), cumul_t(2)/(2*m_max*np)
 IF(rang==0) WRITE(*,*) ' Temps de transfert    ', cumul_t(3), cumul_t(3)/(2*m_max*np)
 t = MPI_WTIME() - t
 IF(rang==0) WRITE(*,*) ' Temps total ', t, SUM(temps)

 DO n = 1, np
    DO i = 1, m_max_c
       DO nf = 1, nb_field/2
          error = MAX(error,ABS(REAL(uc(rang*m_max_c+i,n,nf),KIND=8)-u3(n,2*nf-1,i)) &
                          + ABS(AIMAG(uc(rang*m_max_c+i,n,nf))-u3(n,2*nf,i)))
       END DO
    END DO
 END DO


 WRITE(*,*) ' ERROR ', error 
 CALL MPI_FINALIZE (code)

 END PROGRAM test
