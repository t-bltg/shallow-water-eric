program alltoall
 implicit none

 INCLUDE 'mpif.h'
 INCLUDE 'fftw3.f'

 INTEGER, PARAMETER :: m_max=6, np=7 
 INTEGER :: np_tot, nb_valeurs, m_max_c, bloc_size, N_r
 INTEGER :: nb, shiftc, shiftl, index, jindex, nb_field
 INTEGER :: nb_procs, rang, longueur_tranche, i, n, p, code, MPID
 REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: dist_field, combined_field 
 REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: u1, u2, v1, v2 
 COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: cu, prod_cu, combined_prod_cu, dist_prod_cu, out_prod_cu 
 REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:)  :: ru, prod_ru
 COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: uuc, uc, vc
 REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:)    :: uur, ur, vr 

 ! FFTW parameters
 INTEGER*8 :: fftw_plan_multi_c2r, fftw_plan_multi_r2c
 INTEGER*8 :: fftw_plan_test_c2r, fftw_plan_test_r2c
 INTEGER   :: fft_dim, howmany, istride, ostride, idist, odist
 INTEGER, DIMENSION(1) :: dim, inembed, onembed
 ! Recal complexes must be rescaled
 ! End FFTW parameters

 character(len=200) :: ifilename, ofilename, procrankstr
 INTEGER ::  out_unit=9, in_unit=10

 call MPI_INIT (code)
 call MPI_COMM_SIZE ( MPI_COMM_WORLD ,nb_procs,code)
 call MPI_COMM_RANK ( MPI_COMM_WORLD ,rang,code)

 N_r=2*m_max-1
 m_max_c=m_max/nb_procs
 IF (MOD(m_max,nb_procs)/=0 .OR. m_max_c==0) THEN
    WRITE(*,*) ' BUG '
    STOP
 END IF

 ! Definition of u1 and u2
 write(procrankstr,*) rang
 procrankstr=adjustl(procrankstr)
 ofilename='dat/out_'//trim(procrankstr)//'.dat'
 ifilename='dat/in_'//trim(procrankstr)//'.dat'
 open(unit=in_unit, file=ifilename, action='write')
 open(unit=out_unit, file=ofilename, action='write')

 ALLOCATE(u1(m_max_c,np))
 ALLOCATE(u2(m_max_c,np))
 ALLOCATE(v1(m_max_c,np))
 ALLOCATE(v2(m_max_c,np))
 DO n =  1, np
    DO i = 1, m_max_c
       !u1(i,n) = (n-1)*m_max_c + i-1 + 10*rang 
       !u2(i,n) = (n-1)*m_max_c + i-1 + 10*rang+100*(rang+1) 
       u1(i,n) = (n+i + rang*m_max_c)*1.d0/i**2
       u2(i,n) = (n-i - rang*m_max_c)*1.d0/i**2
    END DO
 END DO
 v1 = u2
 v2 = u1
 !WRITE(in_unit,*) ' Entrees'
 !WRITE(in_unit,*) ' Real'
 !DO i =1, m_max_c
 !  WRITE(in_unit,'(100(f4.0,2x))') u1(i,:)
 !END DO
 !WRITE(in_unit,*) ' Im'
 !DO i =1, m_max_c
 !  WRITE(in_unit,'(100(f4.0,2x))') u2(i,:)
 !END DO
 ! End definition of u1 and u2

 ! Definition of u1 and u2 with all the modes on the processor
 ALLOCATE(uuc(m_max,np),uc(m_max,np), vc(m_max,np))
 ALLOCATE(uur(N_r,np),ur(N_r,np), vr(N_r,np))

 DO n =  1, np
    DO p = 1, nb_procs
       DO i = 1, m_max_c
          uc(i+(p-1)*m_max_c,n) = cmplx((n+i + (p-1)*m_max_c)*1.d0/i**2,(n-i - (p-1)*m_max_c)*1.d0/i**2)
          vc(i+(p-1)*m_max_c,n) = cmplx((n-i - (p-1)*m_max_c)*1.d0/i**2,(n+i + (p-1)*m_max_c)*1.d0/i**2)
       END DO
    END DO
 END DO
 fft_dim=1; istride=1; ostride=1; 
 idist=N_r;   inembed(1)=N_r; dim(1)=N_r
 odist=m_max; onembed(1)=m_max
 howmany = np
 uuc = uc
 CALL dfftw_plan_many_dft_c2r(fftw_plan_test_c2r, fft_dim, dim, howmany, uuc, &
                onembed, istride, odist, uur, inembed, istride, idist, FFTW_ESTIMATE);
 CALL dfftw_execute(fftw_plan_test_c2r)
 ur = uur
 uuc = vc
 CALL dfftw_execute(fftw_plan_test_c2r)
 vr = uur
 ur = ur*vr
 !WRITE(in_unit,*) ' Input'
 !DO i =1, N_r 
 !  WRITE(in_unit,'(100(f10.5,2x))') ur(i,:)
 !END DO
 CALL dfftw_plan_many_dft_r2c(fftw_plan_test_r2c, fft_dim, dim, howmany, ur, &
                inembed, istride, idist, uc, onembed, ostride, odist, FFTW_ESTIMATE);
 CALL dfftw_execute(fftw_plan_test_r2c)
 uc = uc/N_r !========Scaling
 WRITE(in_unit,*) ' Input'
 DO i =1, m_max
   WRITE(in_unit,'(100(e10.3,2x))') uc(i,:)
 END DO
 ! End definition of u1 and u2 with all the modes on the processor

 IF (MODULO(np,nb_procs)==0) THEN
    bloc_size = np/nb_procs
 ELSE
    bloc_size = np/nb_procs + 1
 END IF
 np_tot = nb_procs*bloc_size

 nb_field = 2

 ALLOCATE(dist_field(m_max_c,2*nb_field,np_tot),combined_field(m_max_c,2*nb_field,np_tot))
 ALLOCATE(cu(m_max,nb_field,bloc_size))
 ALLOCATE(ru(N_r,nb_field,bloc_size))
 ALLOCATE(prod_cu(m_max,nb_field/2,bloc_size))
 ALLOCATE(prod_ru(N_r,nb_field/2,bloc_size))

 dist_field(:,1,1:np) = u1
 dist_field(:,2,1:np) = u2
 dist_field(:,3,1:np) = v1
 dist_field(:,4,1:np) = v2
 IF (np/=np_tot) dist_field(:,:,np+1:np_tot) = 1.d100

 longueur_tranche=bloc_size*m_max_c*nb_field*2

 MPID=MPI_DOUBLE_PRECISION
 call MPI_ALLTOALL (dist_field,longueur_tranche, MPID , combined_field,longueur_tranche, &
                    MPID , MPI_COMM_WORLD ,code)

 DO n = 1, bloc_size 
    DO nb = 1, nb_procs
       shiftc = (nb-1)*bloc_size  
       shiftl = (nb-1)*m_max_c
       jindex = n + shiftc
       cu(shiftl+1:shiftl+m_max_c,1,n) = cmplx(combined_field(:,1,jindex),combined_field(:,2,jindex))
       cu(shiftl+1:shiftl+m_max_c,2,n) = cmplx(combined_field(:,3,jindex),combined_field(:,4,jindex))
    END DO
 END DO 

 !WRITE(out_unit,*) ' sortie'
 !DO i =1, m_max
 !  WRITE(out_unit,'(100(f4.0,2x))') cu(i,1,:)
 !END DO

 fft_dim=1; istride=1; ostride=1; 
 idist=N_r;   inembed(1)=N_r; dim(1)=N_r
 odist=m_max; onembed(1)=m_max
 IF (rang==(nb_procs-1)) THEN
    howmany= (np - bloc_size*(nb_procs-1))*nb_field 
 ELSE 
    howmany=bloc_size*nb_field
 END IF

 CALL dfftw_plan_many_dft_c2r(fftw_plan_multi_c2r, fft_dim, dim, howmany, cu, &
           onembed, ostride, odist, ru, inembed, istride, idist, FFTW_ESTIMATE);
 CALL dfftw_execute(fftw_plan_multi_c2r)

 DO n = 1, nb_field/2
    prod_ru(:,n,:)  = ru(:,2*n-1,:)*ru(:,2*n,:)
 END DO
 !WRITE(out_unit,*) ' sortie'
 !DO i =1, N_r 
 !  WRITE(out_unit,'(100(f10.5,2x))') prod_ru(i,1,:)
 !END DO

 howmany = howmany/2
 call dfftw_plan_many_dft_r2c(fftw_plan_multi_r2c, fft_dim, dim, howmany, prod_ru, &
               inembed, istride, idist, prod_cu, onembed, ostride, odist, FFTW_ESTIMATE);

 CALL dfftw_execute(fftw_plan_multi_r2c)

 prod_cu = prod_cu/N_r !Scaling 
 WRITE(out_unit,*) ' sortie'
 DO i =1, m_max
   WRITE(out_unit,'(100(e10.3,2x))') prod_cu(i,1,:)
 END DO
 
 !Now we need to redistribute the Fourier coefficients of each processor

 ALLOCATE(dist_prod_cu(nb_field/2,bloc_size,m_max), combined_prod_cu(nb_field/2,bloc_size,m_max))
 ALLOCATE(out_prod_cu(m_max_c,np_tot,nb_field/2))

 DO n=1, m_max
    combined_prod_cu(:,:,n)=prod_cu(n,:,:)
 END DO
 longueur_tranche=bloc_size*m_max_c*nb_field
 MPID=MPI_DOUBLE_PRECISION
 call MPI_ALLTOALL (combined_prod_cu,longueur_tranche,MPID, dist_prod_cu,longueur_tranche, &
                    MPID,MPI_COMM_WORLD,code)

 !out_prod_cu(mode,:,field) = dist_prod_cu(field,:,mode,:)

 DO n = 1, bloc_size 
    DO nb = 1, nb_procs
       shiftc = (nb-1)*bloc_size  
       shiftl = (nb-1)*m_max_c
       out_prod_cu(:,n+shiftc,1) = dist_prod_cu(1,n,shiftl+1:shiftl+m_max_c)
    END DO
 END DO 

 WRITE(out_unit,*) ' sortie'
 DO i =1, m_max_c
   WRITE(out_unit,'(100(e10.3,2x))') out_prod_cu(i,:,1)
 END DO


 CLOSE (out_unit)
 CALL MPI_FINALIZE (code)

 end program alltoall
