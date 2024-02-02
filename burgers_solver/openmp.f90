program openmp_burgers_solver
    ! compile and run with:
    !   $ ftn -h omp openmp.f90 -o openmp.exe
    !   $ ./openmp.exe
    implicit none

    real, parameter :: nu=1., final_time=0.1, tolerance=1.e-3, safety_factor=0.1
    integer, parameter :: nodes=240
    real, allocatable :: u(:), u_half(:), half_uu(:)
    real :: dx, dt, time

    time = 0
    call init(global_field=u, initial_function=ten_sin, num_points=nodes, dx=dx)
    allocate(u_half(size(u)), half_uu(size(u)))
    dt = safety_factor * diffusion_stability_limit(nu, dx, order_of_accuracy=2)
    do while (time < final_time)
        half_uu(:) = 0.5*u**2
        u_half(:) = u + (dt/2)*(nu*d_dx2(u,dx) - d_dx(half_uu,dx))
        half_uu(:) = 0.5*u_half**2
        u(:) = u + dt*(nu*d_dx2(u_half,dx) - d_dx(half_uu,dx))
        time = time + dt
    end do
    print *, "u = ", u
contains
    subroutine init(global_field, initial_function, num_points, dx)
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

        real, parameter :: two_pi = 2 * 3.14159265
        integer :: i

        allocate(global_field(num_points))
        dx = two_pi / real(num_points)
        associate(grid => [(i - 1, i = 1, num_points)]*dx)
            !$omp parallel do
            do i = 1, num_points
                global_field(i) = initial_function(grid(i))
            end do
        end associate
    end subroutine

    function d_dx(global_field, dx)
        real, intent(in) :: global_field(:)
        real, intent(in) :: dx
        real :: d_dx(size(global_field))

        integer :: i

        associate(num_points => size(global_field))
            d_dx(1) = (global_field(2) - global_field(num_points)) / (2*dx)
            !$omp parallel do
            do i = 2, num_points-1
                d_dx(i) = (global_field(i+1) - global_field(i-1)) / (2*dx)
            end do
            d_dx(num_points) = (global_field(1) - global_field(num_points-1)) / (2*dx)
        end associate
    end function

    function d_dx2(global_field, dx)
        real, intent(in) :: global_field(:)
        real, intent(in) :: dx
        real :: d_dx2(size(global_field))

        integer :: i

        associate(num_points => size(global_field))
            d_dx2(1) = (global_field(2) &
                    - 2*global_field(1) &
                    + global_field(num_points)) / dx**2
            !$omp parallel do
            do i = 2, num_points-1
                d_dx2(i) = (global_field(i+1) - 2*global_field(i) + global_field(i-1)) / dx**2
            end do
            d_dx2(num_points) = (global_field(1) &
                    - 2*global_field(num_points) &
                    + global_field(num_points-1)) / dx**2
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