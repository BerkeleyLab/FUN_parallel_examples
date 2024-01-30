program openmp
    ! compile and run with:
    !   $ ftn -h omp openmp.f90 -o openmp.exe
    !   $ ./openmp.exe
    use iso_fortran_env, only: int64
    implicit none

    integer, parameter :: n = 200, m = 300, k = 400, num_cycles = 100
    real :: a(n, k), b(k, m), c(n, m)
    integer :: i, seed_size
    integer(int64) :: begin, finish, rate, count_max, elapsed

    call random_seed(size=seed_size)
    call random_seed(put=[(i, i = 1, seed_size)])
    call random_number(a)
    call random_number(b)

    call system_clock(count=begin, count_rate=rate, count_max=count_max)
    do i = 1, num_cycles
        call matrix_multiply(a, b, c)
    end do
    call system_clock(count=finish)
    elapsed = finish - begin
    if (elapsed < 0) elapsed = elapsed + count_max

    print *, c(1,1) ! sanity check and prevent optimizing away calculations
    print *, "Took: ", real(elapsed) / rate, " seconds"
contains
    subroutine matrix_multiply(a, b, c)
        real, intent(in) :: a(:,:), b(:,:)
        real, intent(out) :: c(:,:)

        integer :: i, j, k

        if (size(a,1) /= size(c,1)) error stop "a and c must have same number of rows"
        if (size(b,2) /= size(c,2)) error stop "b and c must have same number of columns"
        if (size(a,2) /= size(b,1)) error stop "a must have same number of columns as b has number of rows"

        c = 0.0

        !$omp parallel do collapse(3) reduction(+:c)
        do i = 1, size(a,1)
            do j = 1, size(b, 2)
                do k = 1, size(a,2)
                    c(i, j) = c(i, j) + a(i,k)*b(k,j)
                end do
            end do
        end do
    end subroutine
end program