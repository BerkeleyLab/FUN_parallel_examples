program mpi_burgers_solver
    ! compile and run with:
    !   $ ftn mpi.f90 -o mpi.exe
    !   $ srun -N 1 -n 4 -t 00:00:10 --account=your_project --qos=debug --constraint=cpu ./mpi.exe
    use mpi_f08

    implicit none

    real, parameter :: nu=1., final_time=0.1, tolerance=1.e-3, safety_factor=0.1
    integer, parameter :: nodes=240
    real, allocatable :: u(:), u_half(:), half_uu(:), tmp_d_dx(:), tmp_d_dx2(:)
    real :: dx, dt, time
    integer :: me, nproc, stat

    time = 0
    call init(global_field=u, initial_function=ten_sin, num_points=nodes, dx=dx, me=me, nproc=nproc)
    allocate(u_half(size(u)), half_uu(size(u)), tmp_d_dx(size(u)), tmp_d_dx2(size(u)))
    dt = safety_factor * diffusion_stability_limit(nu, dx, order_of_accuracy=2)
    do while (time < final_time)
        half_uu(:) = 0.5*u**2
        tmp_d_dx(:) = d_dx(half_uu,dx,me,nproc)
        tmp_d_dx2(:) = d_dx2(u,dx,me,nproc)
        u_half(:) = u + (dt/2)*(nu*tmp_d_dx2 - tmp_d_dx)
        half_uu(:) = 0.5*u_half**2
        tmp_d_dx(:) = d_dx(half_uu,dx,me,nproc)
        tmp_d_dx2(:) = d_dx2(u_half,dx,me,nproc)
        u(:) = u + dt*(nu*tmp_d_dx2 - tmp_d_dx)
        time = time + dt
    end do
    print *, "On rank ", me, " u = ", u
    call mpi_finalize(stat)
contains
    subroutine init(global_field, initial_function, num_points, dx, me, nproc)
        real, allocatable, intent(out) :: global_field(:)
        interface
            pure function initial_function(x)
                implicit none
                real, intent(in) :: x
                real :: initial_function
            end function
        end interface
        integer, intent(in) :: num_points
        real, intent(out) :: dx
        integer, intent(out) :: me, nproc

        real, parameter :: two_pi = 2 * 3.14159265
        integer :: i, stat

        call mpi_init(stat)
        call mpi_comm_size(mpi_comm_world, nproc, stat)
        call mpi_comm_rank(mpi_comm_world, me, stat)
        if (mod(num_points, nproc) /= 0) error stop "num_points not evenly divisible by nproc"
        associate(num_local_points => num_points / nproc)
            allocate(global_field(num_local_points))
            dx = two_pi / real(num_points)
            associate(local_grid => [(me*num_local_points + i - 1, i = 1, num_local_points)]*dx)
                do i = 1, num_local_points
                    global_field(i) = initial_function(local_grid(i))
                end do
            end associate
        end associate
    end subroutine

    function d_dx(global_field, dx, me, nproc)
        real, intent(in) :: global_field(:)
        real, intent(in) :: dx
        integer, intent(in) :: me, nproc
        real :: d_dx(size(global_field))

        integer :: i
        real :: left_border, right_border
        type(mpi_status) :: stat

        associate(num_local_points => size(global_field))
            associate( &
                    left_neighbor => merge(nproc-1, me-1, me == 0), &
                    right_neighbor => merge(0, me+1, me==nproc-1))
                if (me < nproc-1) call mpi_send(global_field(num_local_points), 1, MPI_REAL, right_neighbor, 0, mpi_comm_world)
                call mpi_recv(left_border, 1, MPI_REAL, left_neighbor, 0, mpi_comm_world, stat)
                if (me == nproc-1) call mpi_send(global_field(num_local_points), 1, MPI_REAL, right_neighbor, 0, mpi_comm_world)
                if (me > 0) call mpi_send(global_field(1), 1, MPI_REAL, left_neighbor, 0, mpi_comm_world)
                call mpi_recv(right_border, 1, MPI_REAL, right_neighbor, 0, mpi_comm_world, stat)
                if (me == 0) call mpi_send(global_field(1), 1, MPI_REAL, left_neighbor, 0, mpi_comm_world)
            end associate
            d_dx(1) = (global_field(2) - left_border) / (2*dx)
            do i = 2, num_local_points-1
                d_dx(i) = (global_field(i+1) - global_field(i-1)) / (2*dx)
            end do
            d_dx(num_local_points) = (right_border - global_field(num_local_points-1)) / (2*dx)
        end associate
    end function

    function d_dx2(global_field, dx, me, nproc)
        real, intent(in) :: global_field(:)
        real, intent(in) :: dx
        integer, intent(in) :: me, nproc
        real :: d_dx2(size(global_field))

        integer :: i
        real :: left_border, right_border
        type(mpi_status) :: stat

        associate(num_local_points => size(global_field))
            associate( &
                    left_neighbor => merge(nproc-1, me-1, me == 0), &
                    right_neighbor => merge(0, me+1, me==nproc-1))
                if (me < nproc-1) call mpi_send(global_field(num_local_points), 1, MPI_REAL, right_neighbor, 0, mpi_comm_world)
                call mpi_recv(left_border, 1, MPI_REAL, left_neighbor, 0, mpi_comm_world, stat)
                if (me == nproc-1) call mpi_send(global_field(num_local_points), 1, MPI_REAL, right_neighbor, 0, mpi_comm_world)
                if (me > 0) call mpi_send(global_field(1), 1, MPI_REAL, left_neighbor, 0, mpi_comm_world)
                call mpi_recv(right_border, 1, MPI_REAL, right_neighbor, 0, mpi_comm_world, stat)
                if (me == 0) call mpi_send(global_field(1), 1, MPI_REAL, left_neighbor, 0, mpi_comm_world)
            end associate
            d_dx2(1) = (global_field(2) &
                    - 2*global_field(1) &
                    + left_border) / dx**2
            do i = 2, num_local_points-1
                d_dx2(i) = (global_field(i+1) - 2*global_field(i) + global_field(i-1)) / dx**2
            end do
            d_dx2(num_local_points) = (right_border &
                    - 2*global_field(num_local_points) &
                    + global_field(num_local_points-1)) / dx**2
        end associate
    end function

    pure function ten_sin(x)
        real, intent(in) :: x
        real :: ten_sin

        ten_sin = 10 * sin(x)
    end function

    pure function diffusion_stability_limit( &
            diffusivity, delta_x, order_of_accuracy) result(stable_time_step)
        real, intent(in) :: diffusivity, delta_x
        integer, intent(in) :: order_of_accuracy
        real :: stable_time_step

        real, parameter :: stability_limit(*)=[2., 2., 2.5, 2.79]

        stable_time_step = stability_limit(order_of_accuracy)*delta_x**2/(4*diffusivity)
    end function
end program