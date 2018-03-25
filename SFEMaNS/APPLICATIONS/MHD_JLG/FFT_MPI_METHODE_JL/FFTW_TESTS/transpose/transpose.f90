! test of the simplest transpose operation
! for data distributed on a number of processors
! see README for details.
! bash command to see the results aligned :
! for i in dat/out* ; do c=$(cat $i) ; echo ${i#dat/}: $c ; done

program transpose
	use mpi
	implicit none
	integer :: data_unit, result_unit, sendportionsize
	integer :: veclen
    integer :: code, procrank, np
	real(kind=8), allocatable :: xdata(:), ydata(:)
	character(len=200) :: data_filename, result_filename, procrankstr

    call MPI_INIT(code)
    call MPI_COMM_SIZE (MPI_COMM_WORLD, np, code)
    call MPI_COMM_RANK (MPI_COMM_WORLD, procrank, code)

	veclen=3*np
	allocate(xdata(veclen), ydata(veclen))
	write(procrankstr,*) procrank
	procrankstr=adjustl(procrankstr)
	data_filename='dat/in_'//trim(procrankstr)//'.dat'
	result_filename='dat/out_'//trim(procrankstr)//'.dat'
	data_unit=9
	result_unit=11
	open(data_unit, file=data_filename, action="write")
	open(result_unit, file=result_filename, action="write")
	call fill_data(xdata, veclen, procrank) 
	call print_data(xdata, veclen, data_unit)
	ydata=xdata
	sendportionsize=veclen/np
	! exchange the data 
	call MPI_ALLTOALL (xdata, sendportionsize, MPI_DOUBLE_PRECISION, &
		ydata, sendportionsize, MPI_DOUBLE_PRECISION,  MPI_COMM_WORLD ,code)
	call print_data(ydata, veclen, result_unit)
	deallocate(xdata)
	close(data_unit)
	close(result_unit)
    call MPI_FINALIZE(code)
end program transpose

subroutine print_data(vector, veclen, fileunit)
	implicit none
	integer, intent(in) :: fileunit
	integer, intent(in) :: veclen
	real(kind=8), intent(in), dimension(veclen) :: vector
	integer :: i, funit
	funit=fileunit
	do i=1, veclen
		write(funit,'(100f6.0)') vector(i)
	end do
end subroutine print_data

subroutine fill_data(vector, veclen, procrank) 
	implicit none
	integer, intent(in) :: veclen, procrank
    real(kind=8), intent(out), dimension(*) :: vector
	integer :: i
	do i=1, veclen
		vector(i)=10**(procrank+2)+i
	end do
end subroutine fill_data
 
