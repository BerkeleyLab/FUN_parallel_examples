program hello_do_concurrent
    implicit none

    integer, parameter :: num_iterations = 10
    integer :: i

    do concurrent (i = 1 : num_iterations)
        print *, "Hello on iteration ", i
    end do
end program