program transform1D
	! this is a test for calculating 1D Fourier transform with FFTW
	! simplest data set: one sine
	use test_data_1D
	implicit none
	integer :: N=Nx
	integer, parameter :: in_unit=9, out_unit=10
 include 'fftw3.f'

	call fill_sine
	open(unit=in_unit, file='in.dat', action='write')
	call write_arr(in_unit)
	call dfftw_plan_dft_r2c_1d(fftw_plan, N, xdata, xdataout, FFTW_ESTIMATE)
	call dfftw_execute(fftw_plan)
	call print_arr_out
	

end program transform1D
