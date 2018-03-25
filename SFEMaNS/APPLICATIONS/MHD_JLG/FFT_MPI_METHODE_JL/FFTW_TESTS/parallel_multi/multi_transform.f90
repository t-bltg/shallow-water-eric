program transform1D
	! this is a test for calculating multiple 1D Fourier transforms with FFTW and MPI
	use test_data_multi
    use mpi
	implicit none
	integer, parameter :: in_unit=9, out_unit=10
	integer :: rank=1, howmany=n_sets, istride=1, ostride=1, idist=Nx, odist=Nxout
	integer, dimension(1) :: n, inembed, onembed
    integer :: code, np=-1, procrank, tag
    integer, dimension(MPI_STATUS_SIZE) :: status
	character(len=200) :: ifilename, ofilename, procrankstr
 include 'fftw3.f'
    call MPI_INIT(code)
    call MPI_COMM_SIZE (MPI_COMM_WORLD, np, code)
    call MPI_COMM_RANK (MPI_COMM_WORLD, procrank, code)
	rank=1
	print *, "hello from proc number", procrank
	howmany=n_sets
	istride=1
	ostride=1
	idist=Nx
	odist=Nxout
	n(1)=Nx
	inembed(1)=Nx
	onembed(1)=Nxout
	m=procrank
	call fill_data_m
	write(procrankstr,*) procrank
	procrankstr=adjustl(procrankstr)
	ifilename='dat/in_'//trim(procrankstr)//'.dat'
	ofilename='dat/out_'//trim(procrankstr)//'.dat'
	open(unit=in_unit, file=ifilename, action='write')
	open(unit=out_unit, file=ofilename, action='write')
	call dfftw_plan_many_dft_r2c(fftw_plan_multi, rank, n, howmany, mdata, & 
		inembed, istride, idist, mdataout, onembed, ostride, odist, FFTW_ESTIMATE);

	call dfftw_execute(fftw_plan_multi)
	call print_data_m(in_unit)
	call print_out_m(out_unit)
	close(in_unit)
	close(out_unit)
    call MPI_FINALIZE(code)
end program transform1D
