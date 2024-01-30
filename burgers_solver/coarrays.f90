program coarray_burgers_solver
    implicit none

    real, parameter :: nu=1., final_time=0.1, tolerance=1.e-3, safety_factor=0.1
    integer, parameter :: nodes=240
    real, allocatable :: u(:)[:], u_half(:)[:], half_uu(:)[:]
    real :: dx, dt, time

    time = 0
    call init(global_field=u, initial_function=ten_sin, num_points=nodes, dx=dx)
    allocate(u_half(size(u))[*], half_uu(size(u))[*])
    dt = safety_factor * diffusion_stability_limit(nu, dx, order_of_accuracy=2)
    do while (time < final_time)
        call assign_and_synchronize(lhs=half_uu, rhs=0.5*u**2)
        call assign_and_synchronize(lhs=u_half, rhs=u + (dt/2)*(nu*d_dx2(u,dx) - d_dx(half_uu,dx)))
        call assign_and_synchronize(lhs=half_uu, rhs=0.5*u_half**2)
        call assign_and_synchronize(lhs=u, rhs=u + dt*(nu*d_dx2(u_half,dx) - d_dx(half_uu,dx)))
        time = time + dt
    end do
    print *, "On image ", this_image(), " u = ", u
contains
    subroutine init(global_field, initial_function, num_points, dx)
        real, allocatable, intent(inout) :: global_field(:)[:]
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

        if (mod(num_points, num_images()) /= 0) error stop "num_points not evenly divisible by num_images()"
        associate(num_local_points => num_points / num_images())
            allocate(global_field(num_local_points)[*])
            dx = two_pi / real(num_points)
            associate(local_grid => [((this_image()-1)*num_local_points + i - 1, i = 1, num_local_points)]*dx)
                do i = 1, num_local_points
                    global_field(i) = initial_function(local_grid(i))
                end do
            end associate
        end associate
        call synchronize
    end subroutine

    subroutine synchronize
        if (num_images() > 1) then
            associate(me => this_image(), first_image => 1, last_image => num_images())
                associate( &
                        left_neighbor => merge(last_image, me-1, me == first_image), &
                        right_neighbor => merge(first_image, me+1, me == last_image))
                    if (left_neighbor == right_neighbor) then ! occurs if num_images()==2
                        sync images ([left_neighbor])
                    else
                        sync images([left_neighbor, right_neighbor])
                    end if
                end associate
            end associate
        end if
    end subroutine

    subroutine assign_and_synchronize(lhs, rhs)
        real, intent(out) :: lhs(:)
        real, intent(in) :: rhs(:)

        lhs = rhs
        call synchronize
    end subroutine

    pure function d_dx(global_field, dx)
        real, intent(in) :: global_field(:)[*]
        real, intent(in) :: dx
        real :: d_dx(size(global_field))

        integer :: i

        associate(me => this_image(), num_local_points => size(global_field))
            associate(left_neighbor => merge(num_images(), me-1, me == 1))
                d_dx(1) = (global_field(2) - global_field(num_local_points)[left_neighbor]) / (2*dx)
            end associate
            do i = 2, num_local_points-1
                d_dx(i) = (global_field(i+1) - global_field(i-1)) / (2*dx)
            end do
            associate(right_neighbor => merge(1, me+1, me==num_images()))
                d_dx(num_local_points) = (global_field(1)[right_neighbor] - global_field(num_local_points-1)) / (2*dx)
            end associate
        end associate
    end function

    pure function d_dx2(global_field, dx)
        real, intent(in) :: global_field(:)[*]
        real, intent(in) :: dx
        real :: d_dx2(size(global_field))

        integer :: i

        associate(me => this_image(), num_local_points => size(global_field))
            associate(left_neighbor => merge(num_images(), me-1, me == 1))
                d_dx2(1) = (global_field(2) &
                        - 2*global_field(1) &
                        + global_field(num_local_points)[left_neighbor]) / dx**2
            end associate
            do i = 2, num_local_points-1
                d_dx2(i) = (global_field(i+1) - 2*global_field(i) + global_field(i-1)) / dx**2
            end do
            associate(right_neighbor => merge(1, me+1, me==num_images()))
                d_dx2(num_local_points) = (global_field(1)[right_neighbor] &
                        - 2*global_field(num_local_points) &
                        + global_field(num_local_points-1)) / dx**2
            end associate
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