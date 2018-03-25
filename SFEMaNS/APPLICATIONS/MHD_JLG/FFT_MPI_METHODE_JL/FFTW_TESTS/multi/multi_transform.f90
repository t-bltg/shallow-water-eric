program transform1D
	! this is a test for calculating 1D Fourier transform with FFTW
	! simplest data set: one sine
	use test_data_multi
	integer, parameter :: in_unit=9, out_unit=10
	integer :: rank=1, howmany=n_sets, istride=1, ostride=1, idist=Nx, odist=Nxout
	integer, dimension(1) :: n, inembed, onembed
	
	implicit none
 include 'fftw3.f'
	rank=1
	howmany=n_sets
	istride=1
	ostride=1
	idist=Nx
	odist=Nxout
	n(1)=Nx
	inembed(1)=Nx
	onembed(1)=Nxout
	call fill_data_m
	open(unit=in_unit, file='dat/in.dat', action='write')
	open(unit=out_unit, file='dat/out.dat', action='write')
	call dfftw_plan_many_dft_r2c(fftw_plan_multi, rank, n, howmany, mdata, & 
		inembed, istride, idist, mdataout, onembed, ostride, odist, FFTW_ESTIMATE);

!	call dfftw_plan_dft_r2c_1d(fftw_plan, N, xdata, xdataout, FFTW_ESTIMATE)
	call dfftw_execute(fftw_plan_multi)
	call print_data_m(in_unit)
	call print_out_m(out_unit)
	close(in_unit)
	close(out_unit)
end program transform1D
