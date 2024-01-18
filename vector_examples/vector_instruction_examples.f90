program vector_instruction_examples
    ! Compile with:
    !   $ ftn -S -O3 vector_instruction_examples.f90
    ! or compile and run with:
    !   $ ftn -O3 vector_instruction_examples.f90 -o vector_instruction_examples.exe
    !   $ ./vector_instruction_examples.exe
    use iso_fortran_env, only: int64

    implicit none

    integer, parameter :: num_entries = 10000, num_cycles = 100000
    real, dimension(num_entries) :: xs, squares, cubes, fourths
    integer :: i, j
    integer(int64) :: begin, finish, rate, count_max, elapsed
    real :: separate_loop_time, single_loop_time, implied_loop_time, exercise_time

    xs = [(real(i), i = 1, num_entries)]

    call system_clock(count=begin, count_rate=rate, count_max=count_max)
    do j = 1, num_cycles
        call separate_loops
    end do
    call system_clock(count=finish)
    elapsed = finish - begin
    if (elapsed < 0) elapsed = elapsed + count_max
    separate_loop_time = real(elapsed)/rate

    call system_clock(count=begin)
    do j = 1, num_cycles
        call single_loop
    end do
    call system_clock(count=finish)
    elapsed = finish - begin
    if (elapsed < 0) elapsed = elapsed + count_max
    single_loop_time = real(elapsed)/rate

    call system_clock(count=begin)
    do j = 1, num_cycles
        call implied_loop
    end do
    call system_clock(count=finish)
    elapsed = finish - begin
    if (elapsed < 0) elapsed = elapsed + count_max
    implied_loop_time = real(elapsed)/rate

    call system_clock(count=begin)
    do j = 1, num_cycles
        call exercise
    end do
    call system_clock(count=finish)
    elapsed = finish - begin
    if (elapsed < 0) elapsed = elapsed + count_max
    exercise_time = real(elapsed)/rate

    print *, "Timings [seconds]:"
    print *, "  separate loops: ", separate_loop_time
    print *, "  single loop:    ", single_loop_time
    print *, "  implied loops:  ", implied_loop_time
    print *, "  exercise:       ", exercise_time
contains
    subroutine separate_loops
        do i = 1, num_entries
            squares(i) = xs(i)**2
        end do

        do i = 1, num_entries
            cubes(i) = xs(i)**3
        end do

        do i = 1, num_entries
            fourths(i) = xs(i)**4
        end do
    end subroutine

    subroutine single_loop
        do i = 1, num_entries
            squares(i) = xs(i)**2
            cubes(i) = xs(i)**3
            fourths(i) = xs(i)**4
        end do
    end subroutine

    subroutine implied_loop
        squares = [(xs(i)**2, i = 1, num_entries)]
        cubes = [(xs(i)**3, i = 1, num_entries)]
        fourths = [(xs(i)**4, i = 1, num_entries)]
    end subroutine

    subroutine exercise
        do i = 1, num_entries
            squares(i) = xs(i)**2
        end do

        do i = 1, num_entries
            cubes(i) = xs(i)**3
        end do

        do i = 1, num_entries
            fourths(i) = xs(i)**4
        end do
    end subroutine
end program