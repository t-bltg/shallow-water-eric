
program transform1D
	! This program creates two fields: U and V and puts on the output
	! FFT(U), FFT(V and FFT(U*V)
	implicit none
	integer :: N_phy, N_cpl, out_unit, i
	integer*8 :: fftw_plan_u_phy, fftw_plan_v_phy, fftw_plan_uv_phy
	real(kind=8), dimension(:), allocatable :: U_phy, V_phy, UV_phy
	complex(kind=8), dimension(:), allocatable :: U_cpl, V_cpl, UV_cpl
 include 'fftw3.f'
	N_phy=7	
	N_cpl=N_phy/2+1 ! Nxout means number of complex values
	out_unit=10
	allocate(U_phy(N_phy), V_phy(N_phy), UV_phy(N_phy))
	allocate(U_cpl(N_cpl), V_cpl(N_cpl), UV_cpl(N_cpl))

	CALL RANDOM_SEED 
	call random_number (U_phy)
	call random_number (V_phy)
	UV_phy=U_phy*V_phy
	
	call dfftw_plan_dft_r2c_1d(fftw_plan_u_phy, N_phy, U_phy, U_cpl, FFTW_ESTIMATE)
	call dfftw_execute(fftw_plan_u_phy)
	call dfftw_plan_dft_r2c_1d(fftw_plan_v_phy, N_phy, V_phy, V_cpl, FFTW_ESTIMATE)
	call dfftw_execute(fftw_plan_v_phy)
	call dfftw_plan_dft_r2c_1d(fftw_plan_uv_phy, N_phy, UV_phy, UV_cpl, FFTW_ESTIMATE)
	call dfftw_execute(fftw_plan_uv_phy)

	print "(100f9.5)", U_cpl	
	print "(100f9.5)", V_cpl
	print "(100f9.5)", UV_cpl
	print *, 'Data printed above: U_cpl (line 1), V_cpl (line 2), UV_cpl (line 3).'
end program transform1D
