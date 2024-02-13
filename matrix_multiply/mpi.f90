program mpi
    ! compile and run with:
    !   $ ftn mpi.f90 -o mpi.exe
    !   $ srun -N 1 -n 4 -t 00:00:30 --account=your_project --qos=debug --constraint=cpu ./mpi.exe
    use iso_fortran_env, only: int64
    use mpi_f08
    implicit none

    integer, parameter :: n = 200, m = 300, k = 400, num_cycles = 100
    real :: a(n, k), b(k, m), c(n, m)
    integer :: i, seed_size, nproc, me, stat
    integer(int64) :: begin, finish, rate, count_max, elapsed

    call mpi_init(stat)
    call mpi_comm_size(mpi_comm_world, nproc, stat)
    call mpi_comm_rank(mpi_comm_world, me, stat)

    if (me == 0) then
        call random_seed(size=seed_size)
        call random_seed(put=[(i, i = 1, seed_size)])
        call random_number(a)
        call random_number(b)
    end if
    call mpi_bcast(a, size(a), mpi_real, 0, mpi_comm_world, stat)
    call mpi_bcast(b, size(b), mpi_real, 0, mpi_comm_world, stat)

    if (me == 0) call system_clock(count=begin, count_rate=rate, count_max=count_max)
    do i = 1, num_cycles
        call matrix_multiply(a, b, c, me, nproc)
    end do
    if (me == 0) then
        call system_clock(count=finish)
        elapsed = finish - begin
        if (elapsed < 0) elapsed = elapsed + count_max

        print *, c(1,1) ! sanity check and prevent optimizing away calculations
        print *, "Took: ", real(elapsed) / rate, " seconds"
    end if
contains
    subroutine matrix_multiply(a, b, c, me, nproc)
        real, intent(in) :: a(:,:), b(:,:)
        integer, intent(in) :: me, nproc
        real, intent(out) :: c(:,:)

        integer :: i, j, k, stat
        real, allocatable :: local_c(:,:)

        if (size(a,1) /= size(c,1)) error stop "a and c must have same number of rows"
        if (size(b,2) /= size(c,2)) error stop "b and c must have same number of columns"
        if (size(a,2) /= size(b,1)) error stop "a must have same number of columns as b has number of rows"

        allocate(local_c, mold=c)
        local_c = 0.0

        do i = me+1, size(a,1), nproc
            do j = 1, size(b, 2)
                do k = 1, size(a,2)
                    local_c(i, j) = local_c(i, j) + a(i,k)*b(k,j)
                end do
            end do
        end do
        call mpi_reduce(local_c, c, size(c), mpi_real, mpi_sum, 0, mpi_comm_world, stat)
    end subroutine
end program