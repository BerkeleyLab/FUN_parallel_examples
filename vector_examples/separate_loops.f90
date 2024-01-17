program separate_loops
    ! Compile with:
    !   $ ftn -hS -O3 separate_loops.f90
    ! or compile and run with:
    !   $ ftn -O3 separate_loops.f90 -o separate_loops.exe
    !   $ ./separate_loops.exe
    use iso_fortran_env, only: int64

    implicit none

    integer, parameter :: num_entries = 100, num_cycles = 1000
    real, dimension(num_entries) :: xs, squares, cubes, fourths
    integer :: i, j
    integer(int64) :: begin, finish, rate, count_max, elapsed

    xs = [(real(i), i = 1, num_entries)]
    call system_clock(count=begin, count_rate=rate, count_max=count_max)
    do j = 1, num_cycles
        do i = 1, num_entries
            squares(i) = xs(i)**2
        end do

        do i = 1, num_entries
            cubes(i) = xs(i)**3
        end do

        do i = 1, num_entries
            fourths(i) = xs(i)**4
        end do
    end do
    print *, squares, cubes, fourths ! prevents optimizations from removing all code
    call system_clock(count=finish)
    elapsed = finish - begin
    if (elapsed < 0) elapsed = elapsed + count_max
    print *, "Took: ", real(elapsed)/rate, " seconds"
end program