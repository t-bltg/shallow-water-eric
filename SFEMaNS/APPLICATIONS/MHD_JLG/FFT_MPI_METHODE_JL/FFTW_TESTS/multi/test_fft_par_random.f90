program test_fft_par
!
! by Kasia Boronska
!
!  Test of subroutine fftw_par_cross_product,
! with generated input vectors V=(a1,a2,a3) and U=(b1,b2,b3)
! so the expected product is E=VxU=(a2*b3-a3*b2, a3*b1-a1*b3, a1*b2-a2*b1)
! random seed is put to one, so all the processors will generate the same values of a1...b3.
! fdummy() is a dummy subroutine, fixing aggressive optimisation problem on turing.
!
	USE sft_parallele
	IMPLICIT NONE
	INCLUDE 'mpif.h'
	INCLUDE 'fftw3.f'
	INTEGER :: N_phy , np=29, m_max_c, m_max, code, ip, nb_procs, i, rank
	INTEGER(kind=8) :: fftw_plan_r2c 
	INTEGER, ALLOCATABLE, DIMENSION(:) :: int_real_vect
	REAL(kind=8), ALLOCATABLE, DIMENSION(:,:,:) :: V, U, P, E, V_p, U_p, P_p
	REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: real_vect, a1, a2, a3, b1, b2, b3, e1, e2, e3
	REAL(kind=8), DIMENSION(1:3) :: temps
	COMPLEX(kind=8), ALLOCATABLE, DIMENSION(:) :: complex_vect
	

	! Format: V(1:np,1:6,1:m_max_c) 
	m_max=12
	N_phy=m_max*2-1
	
	call MPI_INIT (code)
	call MPI_COMM_RANK ( MPI_COMM_WORLD ,rank,code)
	call MPI_COMM_SIZE ( MPI_COMM_WORLD ,nb_procs,code)
	! check if m_max/nb_procs is an integer
	if(mod(m_max,nb_procs) /=0) then	
		write(*,*) 'Error! m_max should be a multiple of nb_proc'
		stop 1
	end if
	m_max_c=m_max/nb_procs

	!  allocate
	ALLOCATE(V(1:np, 1:6, 1:m_max))
	ALLOCATE(U(1:np, 1:6, 1:m_max))
	ALLOCATE(P(1:np, 1:6, 1:m_max))
	ALLOCATE(E(1:np, 1:6, 1:m_max))
	ALLOCATE(V_p(1:np,1:6,m_max_c), U_p(1:np,1:6,m_max_c), P_p(1:np,1:6,m_max_c))
	ALLOCATE(real_vect(1:N_phy), a1(1:N_phy),  a2(1:N_phy),  a3(1:N_phy),  b1(1:N_phy),  b2(1:N_phy),  b3(1:N_phy))
	ALLOCATE( e1(1:N_phy),  e2(1:N_phy),  e3(1:N_phy))
	ALLOCATE(complex_vect(1:m_max))
	ALLOCATE(int_real_vect(1:N_phy))
 	CALL init_random_seed
	CALL  dfftw_plan_dft_r2c_1d(fftw_plan_r2c, N_phy, real_vect, complex_vect, FFTW_ESTIMATE)

	! fill in input vectors V and U
	! each input vector has three components of length N_phy: a1,...,b3 
	DO ip=1, np
		! fill in randomly
		call random_number(a1)
		call random_number(a2)
		call random_number(a3)
		call random_number(b1)
		call random_number(b2)
		call random_number(b3)
		!a1
		real_vect=a1
		call fdummy(real_vect,complex_vect)
		call dfftw_execute(fftw_plan_r2c)
		V(ip,1,1:m_max)=dble(complex_vect)
		V(ip,2,1:m_max)=imag(complex_vect)
		real_vect=a2
		call fdummy(real_vect,complex_vect);
		call dfftw_execute(fftw_plan_r2c)
		V(ip,3,1:m_max)=dble(complex_vect)
		V(ip,4,1:m_max)=imag(complex_vect)
		! a3
		real_vect=a3
		call fdummy(real_vect,complex_vect); 
		call dfftw_execute(fftw_plan_r2c)
		V(ip,5,1:m_max)=dble(complex_vect)
		V(ip,6,1:m_max)=imag(complex_vect)
		!b1
		real_vect=b1
		call fdummy(real_vect,complex_vect); 
		call dfftw_execute(fftw_plan_r2c)
		U(ip,1,1:m_max)=dble(complex_vect)
		U(ip,2,1:m_max)=imag(complex_vect)
		real_vect=b2
		call fdummy(real_vect,complex_vect);
		call dfftw_execute(fftw_plan_r2c)
		U(ip,3,1:m_max)=dble(complex_vect)
		U(ip,4,1:m_max)=imag(complex_vect)
		! b3
		real_vect=b3
		call fdummy(real_vect,complex_vect); 
		call dfftw_execute(fftw_plan_r2c)
		U(ip,5,1:m_max)=dble(complex_vect)
		U(ip,6,1:m_max)=imag(complex_vect)
		! expected result : 
		DO i=1, N_phy
			e1(i)=a2(i)*b3(i)-a3(i)*b2(i)
			e2(i)=a3(i)*b1(i)-a1(i)*b3(i)
			e3(i)=a1(i)*b2(i)-a2(i)*b1(i)
		END DO
		real_vect=e1
		call fdummy(real_vect,complex_vect); 
		call dfftw_execute(fftw_plan_r2c)
		E(ip,1,1:m_max)=dble(complex_vect)
		E(ip,2,1:m_max)=imag(complex_vect)
		real_vect=e2
		call fdummy(real_vect,complex_vect); 
		call dfftw_execute(fftw_plan_r2c)
		E(ip,3,1:m_max)=dble(complex_vect)
		E(ip,4,1:m_max)=imag(complex_vect)
		real_vect=e3
		call fdummy(real_vect,complex_vect); 
		call dfftw_execute(fftw_plan_r2c)
		E(ip,5,1:m_max)=dble(complex_vect)
		E(ip,6,1:m_max)=imag(complex_vect)
	END DO

	! normalize
	V=V/dble(N_phy)
	U=U/dble(N_phy)
	E=E/dble(N_phy)

	! fill V_p etc. with parts of the total mode set
	V_p(:,:,:)=V(:,:,1+rank*m_max_c:(rank+1)*m_max_c)
	U_p(:,:,:)=U(:,:,1+rank*m_max_c:(rank+1)*m_max_c)
	
	CALL fft_par_cross_prod(V_p, U_p, P_p , temps)
		
	write(*, *) 'rank / max(abs(P-E)):',  rank, '/', maxval(abs( P_p(:,:,:)-E(:,:,1+rank*m_max_c:(rank+1)*m_max_c) ))
	write(*, *) 'rank / maxloc(abs(P-E)):' , rank, '/', maxloc(abs( P_p(:,:,:)-E(:,:,1+rank*m_max_c:(rank+1)*m_max_c) ))

	! deallocate
	DEALLOCATE(U,V,P,E, V_p, U_p, P_p)
	DEALLOCATE(real_vect,complex_vect, int_real_vect , a1, a2, a3, b1, b2, b3, e1, e2, e3)

	CALL MPI_FINALIZE (code)

end program test_fft_par

SUBROUTINE init_random_seed()
			implicit none
			include 'mpif.h'
            INTEGER :: i, n=11, clock
            INTEGER, DIMENSION(:), ALLOCATABLE :: seed
          
            CALL RANDOM_SEED(size = n)
            ALLOCATE(seed(n))
          
			clock=MPI_WTIME() 
			! it is necessary to use a constant if we use more than one processor
			! otherwise the data will be different in each process
			clock=1
            seed = clock + 37 * (/ (i - 1, i = 1, n) /)
            CALL RANDOM_SEED(PUT = seed)
          
            DEALLOCATE(seed)
END SUBROUTINE
