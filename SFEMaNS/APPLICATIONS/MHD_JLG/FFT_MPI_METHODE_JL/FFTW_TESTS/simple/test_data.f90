module dims
	implicit none
	integer, parameter :: Nx =17, Nxout=Nx/2+1
	real(kind=8) :: scale=1
	real(kind=8) :: pi=3.14159265358979323846
end module dims

module test_data_1D
	use dims
	implicit none
	real(kind=8),dimension(Nx) :: xdata
	complex(kind=8), dimension(Nxout) :: xdataout
	integer*8 :: fftw_plan
	contains
	subroutine fill_sine 
		implicit none
		integer :: ix
		do ix=1,Nx
			xdata(ix)=scale*dsin((pi*2.*(ix-1))/Nx)
		end do
	end subroutine fill_sine
	subroutine print_arr
		implicit none
		integer :: ix
		do ix=1, Nx
			print *, xdata(ix)
		end do
	end subroutine print_arr

	subroutine write_arr(unit)
		implicit none
		integer, intent(in) :: unit
		integer :: ix
		do ix=1, Nx
			write(unit, *) xdata(ix)
		end do
	end subroutine write_arr

	subroutine print_arr_out
		implicit none
		integer :: ix
		do ix=1, Nxout
			print *, xdataout(ix)
		end do
	end subroutine print_arr_out
end module test_data_1D
