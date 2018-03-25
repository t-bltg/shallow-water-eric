module dims
	implicit none
	integer, parameter :: Nx =17, Nxout=Nx/2+1
	real(kind=8) :: scale=1
	real(kind=8) :: pi=3.14159265358979323846
end module dims

module test_data_multi
    use dims
    implicit none
    public
    integer, parameter :: n_sets=3
    real(kind=8),dimension(Nx,n_sets) :: mdata
    complex(kind=8), dimension(Nxout,n_sets) :: mdataout
	integer :: m=-10000
    integer*8 :: fftw_plan_multi
    contains
    subroutine fill_data_m
        implicit none
        integer :: ix
        ! fill first set with sine
        do ix=1,Nx
            mdata(ix,1)=scale*sin((pi*2.*(ix-1)*m)/Nx)
        end do
        ! fill second set with a constant
        mdata(:,2)=m*10
        ! fill third set with cosine
        do ix=1,Nx
            mdata(ix,3)=scale*cos((pi*2.*(ix-1)*m)/Nx)
        end do
    end subroutine fill_data_m
    subroutine print_data_m(unit)
        implicit none
        integer :: ix, iset
        integer, intent(in), optional :: unit
        if (present(unit)) then
            do iset=1, n_sets
                do ix=1, Nx
                        write(unit,*) ix, mdata(ix,iset)
                end do
					write(unit,*)
            end do
        else
            do iset=1, n_sets
                do ix=1, Nx
                    print *, ix, mdata(ix,iset)
                end do
                print *
            end do
        end if
    end subroutine print_data_m
   
    subroutine print_out_m(unit)
        implicit none
        integer :: ix, iset
        integer, intent(in), optional :: unit
        if (present(unit)) then
            do iset=1, n_sets
                do ix=1, Nxout
                    write(unit,*) ix, mdataout(ix,iset)
                end do
					write(unit,*)
            end do
        else
            do iset=1, n_sets
                do ix=1, Nxout
                    print *, ix, mdataout(ix,iset)
                end do
                print *
			end do
		end if
    end subroutine print_out_m
end module test_data_multi
