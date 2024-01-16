program hello_mpi
    ! compile and run with:
    !   $ ftn hello_mpi.f90 -o hello_mpi.exe
    !   $ srun -N 1 -n 4 -t 00:00:10 --account=your_project --qos=debug --constraint=cpu ./hello_mpi.exe
    use mpi_f08

    implicit none

    integer :: stat, num_proc, rank

    call mpi_init(stat)
    call check_stat(stat)

    call mpi_comm_size(mpi_comm_world, num_proc, stat)
    call check_stat(stat)

    call mpi_comm_rank(mpi_comm_world, rank, stat)
    call check_stat(stat)

    print *, "Hello from rank ", rank, " of ", num_proc

    call mpi_finalize(stat)
    call check_stat(stat)
contains
    subroutine check_stat(stat_)
        integer, intent(in) :: stat_

        if (stat_ /= 0) stop stat_ ! -\_("/)_/- something went wrong I guess
    end subroutine
end program