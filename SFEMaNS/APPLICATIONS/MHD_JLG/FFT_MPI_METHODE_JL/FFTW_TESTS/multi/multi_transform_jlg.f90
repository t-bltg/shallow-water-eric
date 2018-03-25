program transform1D
        ! this is a test for calculating 1D Fourier transform with FFTW
        ! simplest data set: one sine

        implicit none
        include 'mpif.h'
        include 'fftw3.f'

	integer, parameter :: N_r =16, N_c=N_r/2+1
	real(kind=8) :: scale=1
	real(kind=8) :: pi=3.14159265358979323846d0
        integer, parameter :: np = 3
        real(kind=8),dimension(N_r,np) :: u1, u2 
        complex(kind=8), dimension(N_c,np) :: u1c, u2c 
        complex(kind=8), dimension(N_c,np) :: mdataout
        integer :: ix, n
        integer*8 :: fftw_plan_multi_r2c, fftw_plan_multi_c2r

        integer, parameter :: in_unit=9, out_unit=10
        integer :: rank, howmany, istride, ostride, idist, odist
        integer, dimension(1) :: dim, inembed, onembed
        integer :: code

! Complexe to real do not need scaling.
! Real to complex needs to be scaled. The complex result needs be divided by N_r
        CALL MPI_INIT(code)

        DO n = 1, np
           DO ix = 1, N_c 
              IF (ix==n) THEN
                 IF (ix==1) THEN
                     u2c(ix,n)=cmplx(1.d0,0,8)
                 ELSE
                     u2c(ix,n)=cmplx(1.d0,0,8)/2
                 END IF

              ELSE 
                 u1c(ix,n)=0
              END IF
           end do
	   do ix=1,N_r
              u1(ix,n)=COS((n-1)*(pi*2.d0*(ix-1))/N_r) + SIN((n-1)*(pi*2.d0*(ix-1))/N_r)
           end do
        END DO

        rank=1
        howmany=np
        istride=1
        ostride=1
        idist=N_r
        odist=N_c
        dim(1)=N_r
        inembed(1)=N_r
        onembed(1)=N_c

        call dfftw_plan_many_dft_r2c(fftw_plan_multi_r2c, rank, dim, howmany, u1, &
                inembed, istride, idist, u1c, onembed, ostride, odist, FFTW_ESTIMATE);
        call dfftw_plan_many_dft_c2r(fftw_plan_multi_c2r, rank, dim, howmany, u2c, &
                onembed, istride, odist, u2, inembed, istride, idist, FFTW_ESTIMATE);

        !CALL dfftw_plan_dft_c2r_1d(fftw_plan_multi_c2r, N_r, u1c(:,2), u1(:,2), FFTW_ESTIMATE)
        call dfftw_execute(fftw_plan_multi_r2c)
        call dfftw_execute(fftw_plan_multi_c2r)

        !WRITE(*,*) MAXVAL(ABS(u1-u2/N_r))

        call print_data_m(u2,in_unit)
        u1c=u1c/N_r
        call print_out_m(u1c,out_unit)
        close(out_unit)

        CALL MPI_FINALIZE(code)

CONTAINS
    subroutine print_data_m(mdata,unit)
        implicit none
        integer :: ix, iset
        real(kind=8),dimension(N_r,np) :: mdata
        integer, intent(in), optional :: unit
        if (present(unit)) then
            do iset=1, np
                do ix=1, N_r
                        write(unit,*) ix, mdata(ix,iset)
                end do
					write(unit,*)
            end do
        else
            do iset=1, np
                do ix=1, N_r
                    print *, ix, mdata(ix,iset)
                end do
                print *
            end do
        end if
    end subroutine print_data_m
   
    subroutine print_out_m(mdata,unit)
        implicit none
        integer :: ix, iset
        complex(kind=8), dimension(N_c,np) :: mdata
        integer, intent(in), optional :: unit
        if (present(unit)) then
            do iset=1, np
                do ix=1, N_c
                    write(unit,*) ix, mdata(ix,iset)
                end do
					write(unit,*)
            end do
        else
            do iset=1, np
                do ix=1, N_c
                    print *, ix, mdata(ix,iset)
                end do
                print *
			end do
		end if
    end subroutine print_out_m

end program transform1D
