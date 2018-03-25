program test_fft_par
! Simple test of subroutine fftw_par_cross_product,
! with input vectors V=(a1,a2,a3) and U=(0,0,1)
! so output should be P=VxU=(a2,-a1,0)
! only the zeroth frequencies of each component are non-zero
! so it corresponds to constant values in physical space.
	USE sft_parallele
	IMPLICIT NONE
	INCLUDE 'mpif.h'
	INCLUDE 'fftw3.f'
	INTEGER :: N_phy = 10, np=5, m_max_c, m_max, code, ip
	REAL(kind=8), ALLOCATABLE, DIMENSION(:,:,:) :: V, U, P, E
	REAL(kind=8) :: a1, a2, a3

	! Format: V(1:np,1:6,1:m_max_c) 
	! right now, we have only one processor, so m_max=m_max_c
	m_max=N_phy/2+1
	m_max_c=m_max
	
	call MPI_INIT (code)
	
	!  allocate
	ALLOCATE(V(1:np, 1:6, 1:m_max_c))
	ALLOCATE(U(1:np, 1:6, 1:m_max_c))
	ALLOCATE(P(1:np, 1:6, 1:m_max_c))
	ALLOCATE(E(1:np, 1:6, 1:m_max_c))

	! fill in
	V=0d0
	U=0d0
	E=0d0
	do ip=1, np
		a1=(ip-1)*100+1d0
		a2=(ip-1)*100+2d0
		a3=(ip-1)*100+3d0
		! vector V
		V(ip,1,1)=a1 ! first element of the real part of the first component
		V(ip,3,1)=a2 ! first element of the real part of the second component
		V(ip,5,1)=a3 ! first element of the real part of the third component
		! vector U
		U(ip,1,1)=0d0
		U(ip,3,1)=0d0
		U(ip,5,1)=1d0
		! vector E=VxU
		! (a1,a2,a3)x(0,0,1)=(a2,-a1,0)
		E(ip,1,1)=a2
		E(ip,3,1)=-a1
		E(ip,5,1)=0d0
	end do

	CALL fftw_par_cross_product(V, U, P )
	write(*,*) 'U:' 
	write(*,*) U
	write(*,*) 'V:'
	write(*,*) V
	write(*,*) 'E (expected result):'
	write(*,*) E
	write(*, *) 'P (the result):' 
	write(*,*) P
	write(*, *) 'max(abs(P-E)):' 
	write(*,*) maxval(abs(P-E))

	! deallocate
	DEALLOCATE(U)
	DEALLOCATE(V)
	DEALLOCATE(P)
	DEALLOCATE(E)

	CALL MPI_FINALIZE (code)

end program test_fft_par
