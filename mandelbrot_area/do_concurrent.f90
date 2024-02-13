! Credit to: https://fortran-lang.discourse.group/t/fortran-code-snippets/2150/72
program mandelbrot_area
    ! compile and run with:
    !   ftn -h thread_do_concurrent do_concurrent.f90 -o do_concurrent.exe
    !   $ ./do_concurrent.exe
    implicit none

    integer, parameter :: rk = selected_real_kind(15, 300)
    integer, parameter :: n = 2**12
    real(rk), parameter :: expected = 1.5065849_rk
    integer :: i, j
    real(rk) :: area, re, im

    area = 0

    do concurrent (i = 1 : n, j = 1 : n/2) local(re, im) reduce(+:area)
        re = -2._rk + 2.49_rk * real(i-1,rk) / real(n,rk)
        im = 1.15_rk * real(j-1,rk) / real(n/2,rk)
        if (re**2 + im**2 <= 4._rk) then
            if (julia(cmplx(re,im,rk))) area = area + 1._rk
        end if
    end do

    area = 4._rk * 2.49_rk * 1.15_rk * area / real(n,rk)**2
    print *, "Area: ", area, ", expected: ", expected, ", difference: ", expected-area
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