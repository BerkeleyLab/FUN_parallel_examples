program serial
    ! compile and run with:
    !   $ ftn coarrays.f90 -o coarrays.exe
    !   $ srun -N 1 -n 4 -t 00:00:30 --account=your_project --qos=debug --constraint=cpu ./coarrays.exe
    use iso_fortran_env, only: int64
    implicit none

    integer, parameter :: n = 2, m = 2, k = 2, num_cycles = 1
    real :: a(n, k), b(k, m), c(n, m)
    integer :: i, seed_size, me, nproc
    integer(int64) :: begin, finish, rate, count_max, elapsed

    me = this_image()
    nproc = num_images()

    if (me == 1) then
        call random_seed(size=seed_size)
        call random_seed(put=[(i, i = 1, seed_size)])
        call random_number(a)
        call random_number(b)
    end if
    call co_broadcast(a, 1)
    call co_broadcast(b, 1)

    if (me == 1) call system_clock(count=begin, count_rate=rate, count_max=count_max)
    do i = 1, num_cycles
        call matrix_multiply(a, b, c, me, nproc)
    end do
    if (me == 1) then
        call system_clock(count=finish)
        elapsed = finish - begin
        if (elapsed < 0) elapsed = elapsed + count_max

        print *, c ! sanity check and prevent optimizing away calculations
        print *, "Took: ", real(elapsed) / rate, " seconds"
    end if
contains
    subroutine matrix_multiply(a, b, c, me, nproc)
        real, intent(in) :: a(:,:), b(:,:)
        integer, intent(in) :: me, nproc
        real, intent(out) :: c(:,:)

        integer :: i, j, k

        if (size(a,1) /= size(c,1)) error stop "a and c must have same number of rows"
        if (size(b,2) /= size(c,2)) error stop "b and c must have same number of columns"
        if (size(a,2) /= size(b,1)) error stop "a must have same number of columns as b has number of rows"

        c = 0.0

        do i = me, size(a,1), nproc
            do j = 1, size(b, 2)
                do k = 1, size(a,2)
                    c(i, j) = c(i, j) + a(i,k)*b(k,j)
                end do
            end do
        end do
        call co_sum(c)
    end subroutine
end program