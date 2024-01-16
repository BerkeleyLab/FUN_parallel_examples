program hello_do_concurrent
    ! compile and run with:
    !   $ ftn -h thread_do_concurrent hello_do_concurrent.f90 -o hello_do_concurrent.exe
    !   $ ./hello_do_concurrent.exe
    implicit none

    integer, parameter :: num_iterations = 10
    integer :: i

    do concurrent (i = 1 : num_iterations)
        print *, "Hello on iteration ", i
    end do
end program