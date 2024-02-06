! Credit to: https://fortran-lang.discourse.group/t/fortran-code-snippets/2150/72
program mandelbrot_area
    ! compile and run with:
    !   $ ftn mpi.f90 -o mpi.exe
    !   $ srun -N 1 -n 4 -t 00:00:30 --account=your_project --qos=debug --constraint=cpu ./mpi.exe
    use mpi_f08

    implicit none

    integer, parameter :: rk = selected_real_kind(15, 300)
    integer, parameter :: n = 2**12
    real(rk), parameter :: expected = 1.5065849_rk
    integer :: nproc, me, stat, i, j
    real(rk) :: partial_area, area, re, im

    call mpi_init(stat)
    call mpi_comm_size(mpi_comm_world, nproc, stat)
    call mpi_comm_rank(mpi_comm_world, me, stat)
 
    partial_area = 0

    do i = me+1, n, nproc
        re = -2._rk + 2.49_rk * real(i-1,rk) / real(n,rk)
        do j = 1, n/2
            im = 1.15_rk * real(j-1,rk) / real(n/2,rk)
            if (re**2 + im**2 > 4._rk) cycle
            if (julia(cmplx(re,im,rk))) partial_area = partial_area + 1._rk
        end do
    end do

    call mpi_reduce(partial_area, area, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, stat)

    if (me == 0) then
        area = 4._rk * 2.49_rk * 1.15_rk * area / real(n,rk)**2
        print *, "Area: ", area, ", expected: ", expected, ", difference: ", expected-area
    end if

    call mpi_finalize(stat)
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