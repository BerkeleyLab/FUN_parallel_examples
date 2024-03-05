! Credit to: https://fortran-lang.discourse.group/t/fortran-code-snippets/2150/72
program mandelbrot_area
    ! compile and run with:
    !   $ ftn coarrays.f90 -o coarrays.exe
    !   $ srun -N 1 -n 4 -t 00:00:30 --account=your_project --qos=debug --constraint=cpu ./coarrays.exe

    implicit none

    integer, parameter :: rk = selected_real_kind(15, 300)
    integer, parameter :: n = 2**12
    real(rk), parameter :: expected = 1.5065849_rk
    integer :: nproc, me, i, j
    real(rk) :: area, re, im

    nproc = num_images()
    me = this_image()
 
    area = 0

    do i = me, n, nproc
        re = -2._rk + 2.49_rk * real(i-1,rk) / real(n,rk)
        do j = 1, n/2
            im = 1.15_rk * real(j-1,rk) / real(n/2,rk)
            if (re**2 + im**2 > 4._rk) cycle
            if (julia(cmplx(re,im,rk))) area = area + 1._rk
        end do
    end do

    call co_sum(area, result_image=1)

    if (me == 1) then
        area = 4._rk * 2.49_rk * 1.15_rk * area / real(n,rk)**2
        print *, "Area: ", area, ", expected: ", expected, ", difference: ", expected-area
    end if
contains
    pure function julia(c)
        complex(rk), intent(in) :: c
        logical :: julia

        integer :: i
        complex(rk) :: z

        z=c
        do i = 1, n/8
            z= (((((((z**2 + c)**2 + c)**2 + c)**2 + c)**2 + c)**2 + c)**2 + c)**2 + c
            if (z%re**2 + z%im**2 > 4._rk) then
                julia = .false.
                return
            end if
        end do
        julia=.true.
    end function
end program