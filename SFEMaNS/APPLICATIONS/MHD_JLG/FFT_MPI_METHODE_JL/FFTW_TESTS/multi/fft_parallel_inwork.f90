MODULE sft_parallele

IMPLICIT NONE
  PUBLIC :: FFT_PAR_CROSS_PROD
  PRIVATE

CONTAINS
  SUBROUTINE fft_par_cross_prod(V1_in, V2_in, V_out, temps)
 IMPLICIT NONE

 INCLUDE 'fftw3.f'
 INCLUDE 'mpif.h'
 ! Format: V_1in(1:np,1:6,1:m_max_c) 
 ! INPUT ARE COSINE AND SINE COEFFICIENTS
 ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
 REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN)  :: V1_in, V2_in
 REAL(KIND=8), DIMENSION(:,:,:),  INTENT(OUT) :: V_out
 REAL(KIND=8), DIMENSION(:), OPTIONAL, INTENT(INOUT) :: temps

 ! Saved variables
 LOGICAL, SAVE   :: once=.true.
 INTEGER, SAVE   :: np, np_tot, bloc_size, nb_field, m_max, m_max_c, N_r, rang, nb_procs, MPID
 INTEGER(KIND=8), SAVE :: fftw_plan_multi_c2r, fftw_plan_multi_r2c
 COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: cu, prod_cu, cu_c 
 REAL(KIND=8),    ALLOCATABLE, DIMENSION(:,:,:), SAVE :: ru, prod_ru
 ! End saved variables

 INTEGER :: i_field
 INTEGER :: nb_valeurs, nb, nf, shiftc, shiftl, index, jindex, longueur_tranche, longueur_tranche_c, i, n, p, code
 REAL(KIND=8) :: t
 REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: dist_field, combined_field 
 COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: dist_field_c, combined_field_c
 COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: combined_prod_cu, dist_prod_cu, out_prod_cu 

 ! FFTW parameters
 INTEGER   :: fft_dim, howmany, istride, ostride, idist, odist
 INTEGER, DIMENSION(1) :: dim, inembed, onembed
 ! Recall complexes must be rescaled
 ! End FFTW parameters

 !Temps(1) = Temps de communication
 !Temps(2) = Temps de calcul
 !Temps(3) = Temps de changement de format
 temps = 0.d0

 IF (.NOT.once) THEN
    IF ((SIZE(V1_in,1).NE.np) .OR. (SIZE(V1_in,2).NE.nb_field) .OR. (SIZE(V1_in,3).NE.m_max_c)) THEN
       once = .TRUE. !Something wrong happened, I'am reinitializing
    END IF
 END IF

 IF (once) THEN
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,code)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)

    np      = SIZE(V1_in,1) ! Number of points in the meridian plane
    nb_field= SIZE(V1_in,2) ! Number of fields
    m_max_c = SIZE(V1_in,3) ! Number of complex (cosines + sines) coefficients per point
    m_max = m_max_c*nb_procs! Number of comlex coefficients per point per processor
    N_r=2*m_max-1           ! Number of Real coefficients per point
    IF (MOD(nb_field,2)/=0 .OR. m_max_c==0) THEN
       WRITE(*,*) ' BUG '
       STOP
    END IF

    ! Bloc_size is the number of points that are handled by one processor 
    ! once the Fourier modes are all collected
    ! Computation of bloc_size and np_tot
    IF (MODULO(np,nb_procs)==0) THEN
       bloc_size = np/nb_procs
    ELSE
       bloc_size = np/nb_procs + 1
    END IF
    np_tot = nb_procs*bloc_size
    IF (ALLOCATED(ru)) DEALLOCATE(ru,cu,prod_ru,prod_cu)
    ALLOCATE(cu  (m_max,nb_field,bloc_size))
    ALLOCATE(cu_c(m_max,nb_field,bloc_size))
    ALLOCATE(ru(N_r,  nb_field,bloc_size))
    ALLOCATE(prod_cu(m_max,nb_field/2,bloc_size))
    ALLOCATE(prod_ru(N_r,  nb_field/2,bloc_size))
 END IF

 !t = MPI_WTIME()
 ALLOCATE(      dist_field(m_max_c,2*nb_field,np_tot))
 ALLOCATE(  combined_field(m_max_c,2*nb_field,np_tot))
 ALLOCATE(    dist_field_c(m_max_c,nb_field,np_tot))
 ALLOCATE(combined_field_c(m_max_c,nb_field,np_tot))
 ALLOCATE(dist_prod_cu(nb_field/2,bloc_size,m_max), combined_prod_cu(nb_field/2,bloc_size,m_max))
 ALLOCATE(out_prod_cu(m_max_c,np_tot,nb_field/2))
 !WRITE(*,*) ' time to allocate', MPI_WTIME() -t

 ! Packing all 3 complex components of both v1 and v2 input fields
 ! into dist_field, where the dimension indexing the nodal points varies the least rapidly,
 ! so that after distributing the data to the processes, each one will obtain a part 
 ! on nodal points
 t = MPI_WTIME()
 DO nf = 1, nb_field/2 
    DO i = 1, m_max_c
       dist_field(i,2*nf-1,1:np) = V1_in(:,2*nf-1,i) 
       dist_field(i,2*nf  ,1:np) = V1_in(:,2*nf  ,i) 
       dist_field(i,nb_field+2*nf-1,1:np) = V2_in(:,2*nf-1,i)
       dist_field(i,nb_field+2*nf  ,1:np) = V2_in(:,2*nf,  i)
    END DO
 END DO
 DO nf = 1, nb_field/2 
    ! Put cosine and sine format in complex format for mode 0
    dist_field_c(1,nf,           1:np) = cmplx(V1_in(:,2*nf-1,1),0.d0)
    dist_field_c(1,nb_field/2+nf,1:np) = cmplx(V2_in(:,2*nf-1,1),0.d0)
    DO i = 2, m_max_c
       ! Put cosine and sine format in complex format for nonzero modes
       dist_field_c(i,nf,           1:np) = 0.5d0*cmplx(V1_in(:,2*nf-1,i),-V1_in(:,2*nf,i)) 
       dist_field_c(i,nb_field/2+nf,1:np) = 0.5d0*cmplx(V2_in(:,2*nf-1,i),-V2_in(:,2*nf,i)) 
    END DO
 END DO
 IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

 IF (np/=np_tot) dist_field(:,:,np+1:np_tot) = 1.d100

 longueur_tranche=bloc_size*m_max_c*nb_field*2
 longueur_tranche_c=bloc_size*m_max_c*nb_field

 t = MPI_WTIME()
 MPID=MPI_DOUBLE_PRECISION
 call MPI_ALLTOALL (dist_field,longueur_tranche, MPID , combined_field,longueur_tranche, &
                    MPID , MPI_COMM_WORLD ,code)
 call MPI_ALLTOALL (dist_field_c,longueur_tranche_c, MPID , combined_field_c,longueur_tranche, &
                    MPID , MPI_COMM_WORLD ,code)
 IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t

 t = MPI_WTIME()
 DO n = 1, bloc_size 
    DO nb = 1, nb_procs
       shiftc = (nb-1)*bloc_size  
       shiftl = (nb-1)*m_max_c
       jindex = n + shiftc
       DO nf = 1, nb_field
          ! Put real and imaginary parts in a complex 
          ! for each field, nf=1,2,3,4,5,6
          ! nf=1,2,3 => V1_in
          ! nf=4,5,6 => V2_in
          !cu(shiftl+1:shiftl+m_max_c,nf,n) = cmplx(combined_field(:,2*nf-1,jindex),combined_field(:,2*nf,jindex))
          ! INPUT ARE COSINE AND SINE COEFFICIENTS
          ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
          cu(shiftl+1:shiftl+m_max_c,nf,n) = 0.5d0*cmplx(combined_field(:,2*nf-1,jindex),-combined_field(:,2*nf,jindex))
          cu_c(shiftl+1:shiftl+m_max_c,nf,n) = combined_field_c(:,nf,jindex)
       END DO
    END DO
 END DO 
 cu(1,:,:) = 2*cmplx(REAL(cu(1,:,:)),0.d0)
 write(*,*) ' ERROR ', maxval(ABS(REAL(cu-cu_c)) + ABS(AIMAG(cu-cu_c)))
 STOP
 IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

 ! Set the parameters for dfftw
 fft_dim=1; istride=1; ostride=1; 
 idist=N_r;   inembed(1)=N_r; dim(1)=N_r
 odist=m_max; onembed(1)=m_max
 IF (rang==(nb_procs-1)) THEN
    howmany= (np - bloc_size*(nb_procs-1))*nb_field 
 ELSE 
    howmany=bloc_size*nb_field
 END IF

 t = MPI_WTIME()
 IF (once) CALL dfftw_plan_many_dft_c2r(fftw_plan_multi_c2r, fft_dim, dim, howmany, cu, &
           onembed, ostride, odist, ru, inembed, istride, idist, FFTW_ESTIMATE);
 CALL dfftw_execute(fftw_plan_multi_c2r)

 ! SIMPLE PRODUCT
 !DO nf = 1, nb_field/2
 !   ! ru(1:3) contains V1_in and ru(4:6) contains V2_in
 !   prod_ru(:,nf,:)  = ru(:,nf,:)*ru(:,nb_field/2+nf,:)
 !END DO
 ! END SIMPLE PRODUCT

 ! CROSS PRODDUCT
 IF (nb_field==6) THEN
    prod_ru(:,1,:) = ru(:,2,:)*ru(:,6,:) - ru(:,3,:)*ru(:,5,:)
    prod_ru(:,2,:) = ru(:,3,:)*ru(:,4,:) - ru(:,1,:)*ru(:,6,:)
    prod_ru(:,3,:) = ru(:,1,:)*ru(:,5,:) - ru(:,2,:)*ru(:,4,:)
 END IF
 ! CROSS PRODUCT

 howmany = howmany/2
 IF (once) call dfftw_plan_many_dft_r2c(fftw_plan_multi_r2c, fft_dim, dim, howmany, prod_ru, &
               inembed, istride, idist, prod_cu, onembed, ostride, odist, FFTW_ESTIMATE);
 CALL dfftw_execute(fftw_plan_multi_r2c)

 prod_cu = prod_cu/N_r !Scaling 
 IF (PRESENT(temps)) temps(2) = temps(2) + MPI_WTIME() -t
 
 !Now we need to redistribute the Fourier coefficients on each processor

 t = MPI_WTIME()
 combined_prod_cu(:,:,1)=prod_cu(1,:,:)
 DO n=2, m_max
    ! We want to get back to cosines and sines (a+i*b = 2*conjug(c) for nonzero modes)
    combined_prod_cu(:,:,n)=2*CONJG(prod_cu(n,:,:))
 END DO
 IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

 t = MPI_WTIME()
 longueur_tranche=bloc_size*m_max_c*nb_field
 MPID=MPI_DOUBLE_PRECISION
 call MPI_ALLTOALL (combined_prod_cu,longueur_tranche,MPID, dist_prod_cu,longueur_tranche, &
                    MPID,MPI_COMM_WORLD,code)
 IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t

 ! dimensions: 
 ! out_prod_cu(m_max_c, np_tot, 3)
 ! v_out(np, nb_field, m_max_c)
 ! out_prod_cu is complex 
 t = MPI_WTIME()
 DO i_field = 1, nb_field/2
    DO n = 1, bloc_size 
       DO nb = 1, nb_procs
          shiftc = (nb-1)*bloc_size  
          shiftl = (nb-1)*m_max_c
          out_prod_cu(:,n+shiftc,i_field) = dist_prod_cu(i_field,n,shiftl+1:shiftl+m_max_c)
       END DO
    END DO 
 END DO

 DO i_field = 1, nb_field/2
    DO i = 1, m_max_c
       v_out(:, i_field*2-1, i) = real(out_prod_cu(i, 1:np, i_field))
       v_out(:, i_field*2,   i) = aimag(out_prod_cu(i, 1:np, i_field))
    END DO
 END DO      
 IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

 DEALLOCATE(dist_field, combined_field, dist_prod_cu, combined_prod_cu, out_prod_cu)
 IF (once) once = .FALSE.
 END subroutine FFT_PAR_CROSS_PROD

END MODULE sft_parallele
